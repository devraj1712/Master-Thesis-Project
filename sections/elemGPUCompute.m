%Update the GPU Arrays after each NR iteration

for i=1:nElemBulk %for bulk elements
    gpuxnBulk(:,:,i) = xn(cpConArrayBulk(i, :), 1:2)';
end

for i=1:nElemXi1 %for contact elements
    gpuxnCont(:,:,i) = xn(cpConArrayCont(i, :), 1:2)'; 
end

R   = zeros(nDof, 1);
ksp = zeros(nksz, 1);
I2 = eye(2,"gpuArray");

startTimeContactPart(ls, ns) = toc;

% initialization of element stiffness matrices and internal load vector
gpuKemat   = zeros(nDofCont_Elem, nDofCont_Elem,nElemXi1, ngpv,"gpuArray"); % material 
gpuKegeo   = zeros(nDofCont_Elem, nDofCont_Elem,nElemXi1, ngpv,"gpuArray"); % geometric
gpuFe     = zeros(nDofCont_Elem,nElemXi1, ngpv,"gpuArray");
jterminate_Elem   = 0;

%% Calculation of stiffness matrices and internal force vectors in elementwise manner for contact region

parfor j=1:ngpve
    
    R_1        = pagefun(@mtimes,gpuBasFunContBottomElem(j,:,:),gpuBasFunContVert(j,1,:)); % tensor product of N1 and M to take care first layer forming contact region, size = [1 x (orderXi+1)] 
    R_2        = pagefun(@mtimes,gpuBasFunContTop(j,:,:),gpuBasFunContVert(j,2,:)); % tensor product of N2 and M to take care rest layers forming bulk region, size = [(orderXi+1) x orderEta]
    Ru         = cat(2,R_1,R_2) ; % all possible tensor products of basis functions in two directions layer-wise (row matrix), size = [1 x ((orderXi+1)+(orderXi+1)*orderEta)] = [1 x (orderXi+1)*(orderEta+1)]
    
    dN2xi_a      = pagefun(@mtimes,gpuBasFunDerContBottomElem(:,j,:),gpuBasFunContVert(j,1,:)); % tensor product of dN2_dXi and M, size = [(orderXi+1) x orderEta]
    dM2eta_a     = pagefun(@mtimes,gpuBasFunDerContVert(1,j,:),gpuBasFunContBottomElem(j,:,:)); % tensor product of dN2 and dM_dEta, size = [(orderXi+1) x orderEta]
    dN2xi_b      = pagefun(@mtimes,gpuBasFunDerContTop(:,j,:),gpuBasFunContVert(j,2,:)); 
    dM2eta_b     = pagefun(@mtimes,gpuBasFunDerContVert(2,j,:),gpuBasFunContTop(j,:,:));
    dpR_a = cat(1,dN2xi_a,dN2xi_b) ; %concatenating the 2 vectors
    dpR_b = cat(2,dM2eta_a,dM2eta_b) ;
    dpR_bT = pagefun(@transpose,dpR_b);
    dpR = cat(2,dpR_a,dpR_bT) ; % array of function derivatives in two directions, size =[(orderXi+1)*(orderEta+1) x 2]
    
    first_row = reshape([Ru;zeros(1,size(Ru,2),size(Ru,3))],1,[],size(Ru,3));
    second_row = reshape([zeros(1,size(Ru,2),size(Ru,3));Ru],1,[],size(Ru,3));
    Rb = cat(1,first_row,second_row) ; %Kronecker tensor product of Ru
    RbT = pagefun(@transpose,Rb);
    
    % Jacobian for parameter to physical (current configuraton) size (2x2)
    J1     = pagefun(@mtimes,gpuxnCont,dpR);
    detJ1 = arrayfun(@(u,v,x,y) u*y-v*x ,J1(1,1,:),J1(1,2,:),J1(2,1,:),J1(2,2,:));
    invJ1 = pagefun(@inv,J1);
    if any(detJ1 < 0) 
        jterminate_Elem = 1;
        fprintf('\n Determinant of Jacobian of element %g at time-step %g is not physical.', find(detJ1<0), ls)
    end
    
    % derivative of shape functions in physical space
    dR = pagefun(@mtimes,dpR,invJ1);
    dRT = pagefun(@transpose,dR);

    % strain-displacement matrix
    first_row = reshape([dRT(1,:,:);zeros(1,size(dRT,2),size(dRT,3))],1,[],size(dRT,3));
    second_row = reshape([zeros(1,size(dRT,2),size(dRT,3));dRT(2,:,:)],1,[],size(dRT,3));
    last_row = reshape([dRT(2,:,:);dRT(1,:,:)],1,[],size(dRT,3));
    B = cat(1,first_row,second_row,zeros(1,nDofCont_Elem,nElemXi1),last_row);
    BT = pagefun(@transpose,B);
    
    % deformation gradient
    invFgr = pagefun(@mtimes,gpuXnCont,dR);% invF = dX/dx = Xj*(dR/dx) = XnElem_Elem * dR
    Fgr = pagefun(@inv,invFgr);
    FgrT = pagefun(@transpose,Fgr);
    detFgr = arrayfun(@(u,v,x,y) u*y-v*x ,Fgr(1,1,:),Fgr(1,2,:),Fgr(2,1,:),Fgr(2,2,:));
    if any(detFgr < 0)
        jterminate_Elem = 1;
        fprintf('\n Deformation gradient of element %g at time-step %g is non-physical.', find(detFgr<0), ls)
    end

    % Cauchy stress (in Voigt notation and plane strain condition)  
    % sigma = (lambda/J)*(log J)*I + (mu/J)*(B-I) % Eq. 3.31
    Upr     = arrayfun(@(x,y) x/y * log(y) ,lambda,detFgr);% (lambda/J)*(log J)
    s1 = pagefun(@mtimes,Fgr,FgrT);
    s2 = pagefun(@minus,s1,I2);
    s3 = pagefun(@rdivide,mu,detFgr);
    sit     = pagefun(@mtimes,s3,s2);% (mu/J)*(FF^T-I)
    sig   = cat(1,Upr + sit(1,1,:), Upr + sit(2,2,:), Upr, sit(1,2,:)); % [sigma_11; sigma_22; .......; sigma_12]
    sigma = reshape(sig,4,1,[]);

    % spatial elasticity tensor
    % c = (lambda*(I tensor product I) + 2*mu*(fourth order identity))/J -
    % (lambda/J)*(log J)*2*(fourth order identity) % Eq. 3.32
    Ds      = arrayfun(@(u,v,x,y) u/v - x*y ,D,detFgr,Upr , Is4);
    
    % body forces
    bh    = arrayfun(@(u,v) u*v ,ldv , gpubCont);
    
    % elemental material tangent stiffness matrix
    % int{(B^T).Ds.B}
    k1  = pagefun(@mtimes,detJ1,xigve(j,3));
    k2  = pagefun(@mtimes,k1,detJ2Cont);
    k3  = pagefun(@mtimes,B,k2);
    k4  = pagefun(@mtimes,Ds,k3);
    gpuKemat(:,:,:,j)  = pagefun(@mtimes,BT,k4);
    
    % elemental geometric tangent stiffness matrix
    % int{((dR/dx)^T).(sigma).(dR/dx)}
    a1 = pagefun(@mtimes,dR(:,1,:),dRT(1,:,:));
    a2 = pagefun(@mtimes,dR(:,1,:),dRT(2,:,:));
    a3 = pagefun(@mtimes,dR(:,2,:),dRT(1,:,:));
    a4 = pagefun(@mtimes,dR(:,2,:),dRT(2,:,:));
    b1 = pagefun(@mtimes,a1,sigma(1,1,:));
    b2 = pagefun(@mtimes,a2,sigma(4,1,:));
    b3 = pagefun(@mtimes,a3,sigma(4,1,:));
    b4 = pagefun(@mtimes,a4,sigma(2,1,:));
    NsN   = b1+b2+b3+b4;
    row_stack=[];
    for jj = 1:size(NsN,1)
        first_row = reshape([NsN(jj,:,:);zeros(1,size(NsN,2),size(NsN,3))],1,[],size(NsN,3));
        second_row = reshape([zeros(1,size(NsN,2),size(NsN,3));NsN(jj,:,:)],1,[],size(NsN,3));
        row_stack = cat(1,row_stack,first_row,second_row) ;     
    end
    gpuKegeo(:,:,:,j)  = pagefun(@mtimes,row_stack,k2);
    %Kegeo = Kegeo + kron(NsN, eye(2)) * detJ1 * J2 * wg;
    
    % element load vector
    % int{(B^T).(sigma).(dSigma)}
    f1  = pagefun(@mtimes,BT,sigma);
    f2  = pagefun(@mrdivide,bh,detFgr);
    f3  = pagefun(@mtimes,RbT,f2);
    gpuFe(:,:,j) = pagefun(@mtimes,f1-f3,k2);
    %Fe= (B' * sigma - Rb' * bh / detFgr) * detJ1 * J2 * wg;
    
end

gpuKe = gpuKemat + gpuKegeo;
Ke_dummy = sum (gpuKe,4); %summation over the gauss points
Fe_dummy = sum (gpuFe,3);
clear gpuKemat gpuKegeo gpuKe gpuFe ;

% Compute the contact contribution (adhesive zone)
time_start_contact  = cputime  ;

    % Initial and Current Surface Nodal Coordinates
    gpuXn = gpuXnCont(:,1:orderXiCont+1,1:nxc);
    XIs = pagefun(@transpose,gpuXn);
    gpuxn = gpuxnCont(:,1:orderXiCont+1,1:nxc);
    xIs = pagefun(@transpose,gpuxn);
    clear gpuXn gpuxn

    % Initialize contact force vector and stiffness matrix
    rACi   = zeros(2*(orderXiCont+1),nxc,ngps,"gpuArray") ;
    kACii  = zeros(2*(orderXiCont+1), 2*(orderXiCont+1),nxc,ngps,"gpuArray") ;
    index=[];

    % Loop over Gauss Points
    parfor ii = 1:ngps
        
        % Shape functions
        first_row = reshape([gpuBasFunContBottomSurf(ii,:,:);zeros(1,size(gpuBasFunContBottomSurf,2),size(gpuBasFunContBottomSurf,3))],1,[],size(gpuBasFunContBottomSurf,3));
        second_row = reshape([zeros(1,size(gpuBasFunContBottomSurf,2),size(gpuBasFunContBottomSurf,3));gpuBasFunContBottomSurf(ii,:,:)],1,[],size(gpuBasFunContBottomSurf,3));
        Ni = cat(1,first_row,second_row) ;
        NiT = pagefun(@transpose,Ni); 
        %Ni  = kron(N,I2);
        
        % Coordinates of Gauss points
        xgpi = pagefun(@mtimes,gpuBasFunContBottomSurf(ii,:,:),xIs(:,1,:));
        ygpi = pagefun(@mtimes,gpuBasFunContBottomSurf(ii,:,:),xIs(:,2,:));
        
        % Jacobian Transformation
        dSdxi = pagefun(@mtimes,gpuBasFunDerContBottomSurf(ii,:,:),XIs(:,1,:));
        
        % For Straight projection, rs is
        rkP  = pagefun(@minus,ygpi,yo);
        
        % Force
        Fo   = arrayfun(@(x,y,z) x/(z.^9) - y/(z.^3) ,c1,c2,rkP);
        %c1/rkP^9    - c2/rkP^3 ;
        dFo = arrayfun(@(x,y,z) (-9*x)/(z.^10) + (3*y)/(z.^4) ,c1,c2,rkP);
        %-9*c1/rkP^10 + 3*c2/rkP^4 ;
        
        % Modified Force,
        index = find(rkP < rm);
        Fo(1,1,index(1):index(end))  = arrayfun(@(u,v,x,y) u + v*(x-y) ,Fm,dFm,rkP(1,1,index(1):index(end)),rm);
        %Fm + dFm*(rkP-rm) ;
        dFo(1,1,index(1):index(end)) = dFm        ;
        
        % Body Force and associated Gradient
        bACi  = pagefun(@mtimes,Fo,gpurbk);
        temp = pagefun(@mtimes,gpurbk,gpurbkT);
        dbACk = pagefun(@mtimes,-dFo,temp);
         
        % Residual and Stiffness
        f1  = pagefun(@mtimes,xigs(ii,2),J2c);
        f2  = pagefun(@mtimes,dSdxi,f1);
        f3  = pagefun(@mtimes,bACi,f2);
        rACi(:,:,ii) = pagefun(@mtimes,NiT,f3);
        %rACi  = Ni' * bACi       * dSdxi * wg * J2 ;
        k1 = pagefun(@mtimes,Ni,f2);
        k2 = pagefun(@mtimes,dbACk,k1);
        kACii(:,:,:,ii) = pagefun(@mtimes,NiT,k2);
        %kACii = Ni' * dbACk * Ni * dSdxi * wg * J2 ;
        
    end
    
    KeCont = sum(kACii,4);%summation over gauss points
    FeCont = sum(rACi,3);
    
    % Updating values of element stiffness matrix and force vector
    iS          = (1:2*(orderXiCont+1))      ;
    Fe_dummy(iS,1:nxc)   = Fe_dummy(iS,1:nxc)  + FeCont  ;
    Ke_dummy(iS,iS,1:nxc)  = Ke_dummy(iS,iS,1:nxc) + KeCont ;

    Ke = gather(Ke_dummy);%transferring GPU array to CPU
    Fe = gather(Fe_dummy);
    KeVec = reshape(Ke,[],size(Ke,3)); % storing all the entries in a single column
    clear kACii rACi KeCont FeCont Ke_dummy Fe_dummy

    % End of contact contribution file

total_time_contact_inc_ele = cputime - time_start_contact ;

if jterminate_Elem > 0
    jterminate = 2;
    fprintf('\n Error in computation inside contact element loop ...\n Exiting the Newton iteration loop ...\n') ;
    error('jterminate > 0')
end

ie = zeros(1,nDofCont_Elem);
%% Assembly of Ke and Fe to their global parts for contact part
for i = 1:nElemXi1 % loop over number of contact elements
    cons          = 2 * cpConArrayCont(i, :);
    ie(1, 1:2:nDofCont_Elem) = cons - 1; % odd entries of DOF array for ith element
    ie(1, 2:2:nDofCont_Elem) = cons; % even entries of DOF array for ith element
    R(ie, 1)                 = R(ie, 1) + Fe(:, i); % storing entries in the global internal force vector for ith element
    ksp( (i-1) * nDofCont_Elem^2 + ...
       1 :  i * nDofCont_Elem^2, 1) ...
        = KeVec(:, i); % storing entries in the global sparse stiffness matrix for ith element
end

clear Ke Fe KeVec ie  

% computation time taken for the contact part
endTimeContactPart(ls,ns) = toc - startTimeContactPart(ls, ns);


%% Calculation of stiffness matrices and internal force vectors in elementwise manner for bulk region

time1         = toc;
startTimeBulkPart(ls, ns) = toc;

% initialization of element stiffness matrices and load vector
gpuKemat   = zeros(nDofBulk_Elem, nDofBulk_Elem,nElemBulk, ngpv,"gpuArray"); % material 
gpuKegeo   = zeros(nDofBulk_Elem, nDofBulk_Elem,nElemBulk, ngpv,"gpuArray"); % geometric
gpuFe      = zeros(nDofBulk_Elem,nElemBulk, ngpv,"gpuArray");
jterminate_Elem   = 0;

parfor j=1:ngpv
    
    % Kronecker tensor product of basis function with I2, size = kron([1 x (orderXi+1)*(orderEta+1)] and [2 x 2]) = [2 x 2*(orderXi+1)*(orderEta+1)]
    first_row = reshape([gpuBasFunBulk(j,:,:);zeros(1,size(gpuBasFunBulk,2),size(gpuBasFunBulk,3))],1,[],size(gpuBasFunBulk,3));
    second_row = reshape([zeros(1,size(gpuBasFunBulk,2),size(gpuBasFunBulk,3));gpuBasFunBulk(j,:,:)],1,[],size(gpuBasFunBulk,3));
    Rb = cat(1,first_row,second_row) ;%concatenating the 2 row vectors
    RbT = pagefun(@transpose,Rb);
    
    dpR = cat(2,gpuBasFunDerXiBulk(:,j,:),gpuBasFunDerEtaBulk(:,j,:)) ; % array of function derivatives in two directions, size =[(orderXi+1)*(orderEta+1) x 2]
    
    % Jacobian for parameter to physical (current configuraton)
    J1     = pagefun(@mtimes,gpuxnBulk,dpR);
    detJ1 = arrayfun(@(u,v,x,y) u*y-v*x ,J1(1,1,:),J1(1,2,:),J1(2,1,:),J1(2,2,:));
    invJ1 = pagefun(@inv,J1);
    if any(detJ1 < 0) 
        jterminate_Elem = 1;
        fprintf('\n Determinant of Jacobian of element %g at time-step %g is not physical.', find(detJ1<0), ls)
    end

    % derivative of shape functions in physical space
    dR = pagefun(@mtimes,dpR,invJ1);
    dRT = pagefun(@transpose,dR);
    
    % strain-displacement matrix
    first_row = reshape([dRT(1,:,:);zeros(1,size(dRT,2),size(dRT,3))],1,[],size(dRT,3));
    second_row = reshape([zeros(1,size(dRT,2),size(dRT,3));dRT(2,:,:)],1,[],size(dRT,3));
    last_row = reshape([dRT(2,:,:);dRT(1,:,:)],1,[],size(dRT,3));    
    B = cat(1,first_row,second_row,zeros(1,nDofBulk_Elem,nElemBulk),last_row);
    BT = pagefun(@transpose,B);
    
    % deformation gradient
    invFgr = pagefun(@mtimes,gpuXnBulk,dR); % invF = dX/dx = Xj*(dR/dx) = XnElem_Elem * dR
    Fgr = pagefun(@inv,invFgr);
    FgrT = pagefun(@transpose,Fgr);
    detFgr = arrayfun(@(u,v,x,y) u*y-v*x ,Fgr(1,1,:),Fgr(1,2,:),Fgr(2,1,:),Fgr(2,2,:));    
    if any(detFgr < 0)
        jterminate_Elem = 1;
        fprintf('\n Deformation gradient of element %g at time-step %g is non-physical.', find(detFgr<0), ls)
    end

    % Cauchy stress (in Voigt notation and plane strain condition)  
    % sigma = (lambda/J)*(log J)*I + (mu/J)*(B-I) % Eq. 3.31
    Upr     = arrayfun(@(x,y) x/y * log(y) ,lambda,detFgr);% (lambda/J)*(log J)
    s1 = pagefun(@mtimes,Fgr,FgrT);
    s2 = pagefun(@minus,s1,I2);
    s3 = pagefun(@rdivide,mu,detFgr);
    sit     = pagefun(@mtimes,s3,s2);% (mu/J)*(FF^T-I)    
    sig   = cat(1,Upr + sit(1,1,:), Upr + sit(2,2,:), Upr, sit(1,2,:)); % [sigma_11; sigma_22; .......; sigma_12]
    sigma = reshape(sig,4,1,[]);

    % spatial elasticity tensor
    % c = (lambda*(I tensor product I) + 2*mu*(fourth order identity))/J -
    % (lambda/J)*(log J)*2*(fourth order identity) % Eq. 3.32
    Ds      = arrayfun(@(u,v,x,y) u/v - x*y ,D,detFgr,Upr , Is4);
    
    % body forces
    bh    = arrayfun(@(u,v) u*v ,ldv , gpubBulk);
    
    % elemental material tangent stiffness matrix
    % int{(B^T).Ds.B}
    k1  = pagefun(@mtimes,detJ1,xigv(j,3));
    k2  = pagefun(@mtimes,k1,detJ2Bulk);
    k3  = pagefun(@mtimes,B,k2);
    k4  = pagefun(@mtimes,Ds,k3);
    gpuKemat(:,:,:,j)  = pagefun(@mtimes,BT,k4);
    
    % elemental geometric tangent stiffness matrix
    % int{((dR/dx)^T).(sigma).(dR/dx)}
    a1 = pagefun(@mtimes,dR(:,1,:),dRT(1,:,:));
    a2 = pagefun(@mtimes,dR(:,1,:),dRT(2,:,:));
    a3 = pagefun(@mtimes,dR(:,2,:),dRT(1,:,:));
    a4 = pagefun(@mtimes,dR(:,2,:),dRT(2,:,:));
    b1 = pagefun(@mtimes,a1,sigma(1,1,:));
    b2 = pagefun(@mtimes,a2,sigma(4,1,:));
    b3 = pagefun(@mtimes,a3,sigma(4,1,:));
    b4 = pagefun(@mtimes,a4,sigma(2,1,:));
    NsN   = b1+b2+b3+b4;
    row_stack=[];
    for jj = 1:size(NsN,1)
        first_row = reshape([NsN(jj,:,:);zeros(1,size(NsN,2),size(NsN,3))],1,[],size(NsN,3));
        second_row = reshape([zeros(1,size(NsN,2),size(NsN,3));NsN(jj,:,:)],1,[],size(NsN,3));
        row_stack = cat(1,row_stack,first_row,second_row) ;     
    end
    gpuKegeo(:,:,:,j)  = pagefun(@mtimes,row_stack,k2);
    %Kegeo = Kegeo + kron(NsN, eye(2)) * detJ1 * J2 * wg;
    
    % element load vector
    % int{(B^T).(sigma).(dSigma)}
    f1  = pagefun(@mtimes,BT,sigma);
    f2  = pagefun(@mrdivide,bh,detFgr);
    f3  = pagefun(@mtimes,RbT,f2);
    gpuFe(:,:,j) = pagefun(@mtimes,f1-f3,k2);
    %Fe= (B' * sigma - Rb' * bh / detFgr) * detJ1 * J2 * wg;
    
end

gpuKe = gpuKemat + gpuKegeo;
Ke_dummy = gather(gpuKe); % transferring GPU array to CPU
Fe_dummy = gather(gpuFe);
Ke = sum (Ke_dummy,4);% Summation over gauss points
Fe = sum (Fe_dummy,3);
KeVec = reshape(Ke,[],size(Ke,3)); % storing all the entries in a single column
clear gpuKemat gpuKegeo gpuKe  gpuFe Ke_dummy Fe_dummy ;

if jterminate_Elem > 0
    jterminate = 1;
    fprintf('\n Error in computation inside bulk element loop ...\n Exiting the Newton iteration loop ...\n') ;
    error('jterminate > 0')
end

%% Assembly of Ke and Fe to their global parts for bulk part
ie = zeros(1,nDofBulk_Elem);
for i = 1:nElemBulk % loop over number of bulk elements
    cons          = 2 * cpConArrayBulk(i, :);
    ie(1, 1:2:nDofBulk_Elem) = cons - 1; % odd entries of DOF array for ith element
    ie(1, 2:2:nDofBulk_Elem) = cons; % even entries of DOF array for ith element
    R(ie, 1)                 = R(ie, 1) + Fe(:, i); % storing entries in the global internal force vector for ith element
    ksp( nElemXi1*nDofCont_Elem^2 + (i-1) * nDofBulk_Elem^2 + 1 : nElemXi1*nDofCont_Elem^2 +  i * nDofBulk_Elem^2, 1) ...
        = KeVec(:, i); % storing entries in the global sparse stiffness matrix for ith element
end

clear  Ke Fe KeVec

% computation time taken for the bulk part
endTimeBulkPart(ls, ns) = toc - startTimeBulkPart(ls, ns);
