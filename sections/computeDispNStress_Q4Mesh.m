%% Compute stresses at the nodes of projected elements

% Number of control points
nCPCont1 = (orderXiCont+1) + (orderXi+1) * orderEta; % contact region
nCPBulk1 = (orderXi+1) * (orderEta+1); % bulk region

% Control point connectivity array (transposed)
cpConArrayCont1 = cpConArrayCont'; % contact region
cpConArrayBulk1 = cpConArrayBulk'; % bulk region

% Coordinate array
U  = zeros(nDof, 1);
U(1:2:nDof, 1) = xn_storeG(1:nCP, 1, zz) - Xn(:, 1);
U(2:2:nDof, 1) = xn_storeG(1:nCP, 2, zz) - Xn(:, 2);
UX = U(1:2:nDof);
UY = U(2:2:nDof);

stress = zeros(nElem_Proj, 4, 4); % stress of all projected elements, (4, 4) denote the number of nodes in Q4 and number of componenets of stress, respectively
disp = zeros(nElem_Proj, 4, 2); % displacement of all projected elements, (4, 2) denote the number of nodes in Q4 and number of componenets of displacement, respectively
stressNode = zeros(max(cpConArray_Q4Mesh(:,3)), 4); % size = [total number of nodes of all projected elements, number of stress components]
dispNode = zeros(max(cpConArray_Q4Mesh(:,3)), 2); % size = [total number of nodes of all projected elements, number of displacement components]

for i = 1:nElem_Proj % loop over total number of projected elements
    
    % Determine the corresponding actual element to ith projected element
    [r2, c2] = find(elemConArray_Q4Mesh == i);
    I        = r2(end,1); % identity of actual element
    
    % Connectivity array for each quadrilateral element
    if i < (nElem_Proj1+1)
        cons = cpConArrayCont1(:, I);
        nCP  = nCPCont1;
    else
        cons = cpConArrayBulk1(:, I-nElemXi);
        nCP  = nCPBulk1;
    end
    
    Xn_ProjElem      = Xn(cons, :); % undeformed coordinate array for ith projected element
    xn_ProjElem      = xn(cons, :); % deformed coordinate array for ith projected element
    kntSpanXi_ProjElem     = kntSpanXi_Proj(kntConArray_Proj(i, 1), :); % knot span along xi for ith projected element
    kntSpanEta_ProjElem    = kntSpanEta_Proj(kntConArray_Proj(i, 2), :); % knot span along eta for ith projected element
    gp      = 1;
    sel     = zeros(4, 4);
    mll     = zeros(4, 1);
    
    for iv = 1:2 % loop over number of Q4 element nodes along eta
        
        if (iv==2) 
            kntSpanXi_ProjElem = sort(kntSpanXi_ProjElem, 'descend'); % if first layer is [0 0.2], then second layer is [0.2, 0] because [1 2 3 4] -> [1 2 4 3] for plot within Q4
        end
        for iu = 1:2 % loop over number of Q4 element nodes along xi
            xi  = kntSpanXi_ProjElem(iu); % coordinate in parameter space
            eta = kntSpanEta_ProjElem(iv); % coordinate in parameter space    
            
            % NURBS basis functions and their derivatives
            if i < (nElem_Proj1+1)
                [N1]             = nrbbasisfun(xi, crvContBottom); % function value at the bottom layer of contact region
                [N2]             = nrbbasisfun(xi, crvContTop); % function value at the top layer of contact region
                [M]              = nrbbasisfun(eta, crvVertical); % function value at the vertical surface of the geometry
                [dN1_dXi]        = nrbbasisfunder(xi, crvContBottom); % function derivatives at the bottom layer of contact region
                [dN2_dXi]        = nrbbasisfunder(xi, crvContTop); % function derivatives at the top layer of contact region
                [dM_dEta]        = nrbbasisfunder(eta, crvVertical); % function derivatives at the vertical surface of the geometry
                R_1              = N1 * M(1,1); % tensor product of N1 and M to take care first layer forming contact region
                R_2              = N2 * M(1,2); % tensor product of N2 and M to take care rest layers forming bulk region
                Ru               = [R_1, R_2]; % all possible tensor products of basis functions in two directions layer-wise (row matrix)
                dN2xi            = dN2_dXi' * M(1, 2); % tensor product of dN2_dXi and M
                dM2eta           =  dM_dEta(1, 2)*N2' ; % tensor product of dN2 and dM_dEta
                dR_dXi           = [dN1_dXi' * M(1,1); dN2xi]';
                dR_dEta          = [ dM_dEta(1,1)*N1' ; dM2eta]';
                dpR              = [dR_dXi', dR_dEta']; 
            else
                Ru                = nrbbasisfun({xi, eta}, nrbObj);
                [dR_dXi, dR_dEta] = nrbbasisfunder({xi, eta}, nrbObj);
                dpR              = [dR_dXi', dR_dEta'];
            end
            
            % Jacobian for parameter to physical (current configuraton)
            J1    = zeros(2, 2) ;
            for j = 1:nCP
                J1 = J1 + xn_ProjElem(j, :)' * dpR(j, :) ;
            end
            detJ1 = det(J1);
            invJ1 = inv(J1);
            
            % derivative of shape functions in physical space
            dR = dpR * (invJ1);
            
            % deformation gradient
            invFgr = zeros(2, 2) ;
            for j = 1:nCP
                invFgr = invFgr + Xn_ProjElem(j, :)' * dR(j, :) ;
            end
            Fgrd  = inv(invFgr);
            detFgr = det(Fgrd);
            if detFgr < 0
                jterminate = 1 ;
                fprintf('\n Deformation gradient of element %g is non-physical.', i)
            end
            
            % Cauchy stress (in Voigt notation and plane strain condition)
            Upr       = lambda / detFgr * log(detFgr);
            sit       = mu / detFgr * (Fgrd * Fgrd' - eye(2));
            sigma     = [Upr + sit(1,1); Upr + sit(2,2); Upr; sit(1,2)] ;
            
            dispNode(cpConArray_Q4Mesh(i, gp), :)   =  Ru * [UX(cons, 1), UY(cons, 1)];
            stressNode(cpConArray_Q4Mesh(i, gp), :) =  sigma';
            gp = gp + 1;
            
        end
    end
end