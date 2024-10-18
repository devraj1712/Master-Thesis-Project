%Create GPU Arrays for NURBS Basis functions and their derivatives

gpuBasFunBulk = zeros(ngpv,prod(nrbObj.order),nElemBulk,"gpuArray");
gpuBasFunDerXiBulk = zeros( prod(nrbObj.order),ngpv,nElemBulk,"gpuArray");
gpuBasFunDerEtaBulk = zeros( prod(nrbObj.order),ngpv,nElemBulk,"gpuArray");
gpuBasFunContTop = zeros(ngpve, crvContTop.order,nElemXi1,"gpuArray");
gpuBasFunContBottomElem = zeros(ngpve, crvContBottom.order,nElemXi1,"gpuArray");
gpuBasFunContBottomSurf = zeros( ngps,crvContBottom.order,nxc,"gpuArray");
gpuBasFunContVert = zeros(ngpve, crvVertical.order,nElemXi1,"gpuArray");
gpuBasFunDerContTop = zeros( crvContTop.order,ngpve,nElemXi1,"gpuArray");
gpuBasFunDerContBottomElem = zeros( crvContBottom.order,ngpve,nElemXi1,"gpuArray");
gpuBasFunDerContBottomSurf = zeros( ngps,crvContBottom.order,nxc,"gpuArray");
gpuBasFunDerContVert = zeros( crvVertical.order,ngpve,nElemXi1,"gpuArray");

%Define GPU Arrays for other parameters for the functions

gpuKntSpanXiBulk = zeros(nElemBulk,2,"gpuArray");
gpuKntSpanEtaBulk = zeros(nElemBulk,2,"gpuArray");
gpuKntSpanXiCont = zeros(nElemXi1,2,"gpuArray");
gpuKntSpanEtaCont = zeros(nElemXi1,2,"gpuArray");
gpuXnBulk = zeros(2,size(cpConArrayBulk,2),nElemBulk,"gpuArray");
gpuxnBulk = zeros(2,size(cpConArrayBulk,2),nElemBulk,"gpuArray");
gpuXnCont = zeros(2,size(cpConArrayCont,2),nElemXi1,"gpuArray");
gpuxnCont = zeros(2,size(cpConArrayCont,2),nElemXi1,"gpuArray");
gpubBulk    = zeros(2,1,nElemBulk,"gpuArray");
gpubCont    = zeros(2,1,nElemXi1,"gpuArray");
gpurbk_1  = zeros(1,1,nxc,"gpuArray");
gpurbk_2  = ones(1,1,nxc,"gpuArray");
gpurbk_3 = -1*gpurbk_2;
gpurbk = cat(1,gpurbk_1,gpurbk_3);
gpurbkT = pagefun(@transpose , gpurbk);
clear gpurbk_1 gpurbk_2 gpurbk_3


%Assign values to the GPU Arrays
parfor i=1:nElemBulk %for bulk elements

    gpuKntSpanXiBulk(i,:) = kntSpanXi(kntConArray(i+nElemXi1, 1), 1:2);
    gpuKntSpanEtaBulk(i,:) = kntSpanEta(kntConArray(i+nElemXi1, 2), 1:2);
    gpuXnBulk(:,:,i) = Xn(cpConArrayBulk(i, :), 1:2)';
    
end

% coordinates in parameter space
xiBulk=zeros(nElemBulk,ngpv);
etaBulk=zeros(nElemBulk,ngpv);
parfor j= 1:ngpv
    % coordinates in parameter space
    xidummy=arrayfun(@(x,y,z) 0.5*((y - x).*z + (y + x)),gpuKntSpanXiBulk(:,1),gpuKntSpanXiBulk(:,2),xigv(j,1));
    etadummy=arrayfun(@(x,y,z) 0.5*((y - x).*z + (y + x)),gpuKntSpanEtaBulk(:,1),gpuKntSpanEtaBulk(:,2),xigv(j,2));
    xiBulk(:,j)=xidummy;
    etaBulk(:,j)=etadummy;
end

% mapping from parameter to master space
dXi_dXiBar    = arrayfun(@(x,y) 0.5 *(y - x),gpuKntSpanXiBulk(:,1),gpuKntSpanXiBulk(:,2) );
dEta_dEtaBar   = arrayfun(@(x,y) 0.5 *(y - x),gpuKntSpanEtaBulk(:,1),gpuKntSpanEtaBulk(:,2) );
J2            = dXi_dXiBar .* dEta_dEtaBar;  % Jacobian for parameter to master space
detJ2Bulk = reshape(J2',1,1,[]);

% NURBS basis functions and their derivatives
parfor i=1:nElemBulk 
    for j= 1:ngpv
        gpuBasFunBulk(j,:,i) =  nrbbasisfun({xiBulk(i,j), etaBulk(i,j)}, nrbObj);
        [gpuBasFunDerXiBulk(:,j,i) , gpuBasFunDerEtaBulk(:,j,i)] = nrbbasisfunder({xiBulk(i,j), etaBulk(i,j)}, nrbObj);
    end
end


parfor i=1:nElemXi1 %for contact elements

    gpuKntSpanXiCont(i,:) = kntSpanXi(kntConArray(i, 1), 1:2);
    gpuKntSpanEtaCont(i,:) = kntSpanEta(kntConArray(i, 2), 1:2);
    gpuXnCont(:,:,i) = Xn(cpConArrayCont(i, :), 1:2)';
    
end

% coordinates in parameter space
xiCont=zeros(nElemXi1,ngpve);
etaCont=zeros(nElemXi1,ngpve);
parfor j= 1:ngpve
    % coordinates in parameter space
    xidummy=arrayfun(@(x,y,z) 0.5*((y - x).*z + (y + x)),gpuKntSpanXiCont(:,1),gpuKntSpanXiCont(:,2),xigve(j,1));
    etadummy=arrayfun(@(x,y,z) 0.5*((y - x).*z + (y + x)),gpuKntSpanEtaCont(:,1),gpuKntSpanEtaCont(:,2),xigve(j,2));
    xiCont(:,j)=xidummy;
    etaCont(:,j)=etadummy;
end

% mapping from parameter to master space
dXi_dXiBar    = arrayfun(@(x,y) 0.5 *(y - x),gpuKntSpanXiCont(:,1),gpuKntSpanXiCont(:,2) );
dEta_dEtaBar   = arrayfun(@(x,y) 0.5 *(y - x),gpuKntSpanEtaCont(:,1),gpuKntSpanEtaCont(:,2) );
J2            = dXi_dXiBar .* dEta_dEtaBar;  % Jacobian for parameter to master space
detJ2Cont = reshape(J2',1,1,[]);

% NURBS basis functions and their derivatives
parfor i=1:nElemXi1 
    for j= 1:ngpve
        % coordinates in parameter space
        gpuBasFunContTop(j,:,i)       = nrbbasisfun(xiCont(i,j), crvContTop); % function values at the bottom layer of contact region, size = [1 x (orderXi+1)]
        gpuBasFunContBottomElem(j,:,i)       = nrbbasisfun(xiCont(i,j), crvContBottom); % function values at the top layer of contact region, size = [1 x (orderXi+1)]
        gpuBasFunContVert(j,:,i)        = nrbbasisfun(etaCont(i,j), crvVertical); % function values at the vertical layer, size = [1 x (orderEta+1)]
        gpuBasFunDerContBottomElem(:,j,i)  = nrbbasisfunder(xiCont(i,j), crvContBottom); % function derivatives at the bottom layer of contact region, size = [1 x (orderXi+1)]
        gpuBasFunDerContTop(:,j,i)  = nrbbasisfunder(xiCont(i,j), crvContTop); % function derivatives at the top layer of contact region, size = [1 x (orderXi+1)]
        gpuBasFunDerContVert(:,j,i)  = nrbbasisfunder(etaCont(i,j), crvVertical); % function derivatives at the vertical layer, size = [1 x (orderEta+1)]
    end
end

% coordinates in parameter space
xic=zeros(nxc,ngps);
gpuKntSpanCont = gpuKntSpanXiCont(1:nxc,:);
parfor j= 1:ngps
    xidummy=arrayfun(@(x,y,z) 0.5*((y - x).*z + (y + x)),gpuKntSpanCont(:,1),gpuKntSpanCont(:,2),xigs(j,1));
    xic(:,j)=xidummy;
end

% Mapping from parent element to Parametric element(for adhesive zone)
dxi_dxitilda = arrayfun(@(x,y) 0.5 *(y - x),gpuKntSpanCont(:,1),gpuKntSpanCont(:,2) );
J2c  = reshape(dxi_dxitilda',1,1,[]);
clear dxi_dxitilda dXi_dXiBar dEta_dEtaBar J2

% NURBS basis functions and their derivatives (for adhesive zone)
parfor i=1:nxc 
    for j= 1:ngps
        gpuBasFunContBottomSurf(j,:,i) = nrbbasisfun(xic(i,j), crvContBottom);
        gpuBasFunDerContBottomSurf(j,:,i) = nrbbasisfunder(xic(i,j), crvContBottom);
    end
end

clear xic xiCont etaCont xiBulk etaBulk
clear gpuKntSpanCont gpuKntSpanXiCont gpuKntSpanEtaCont gpuKntSpanXiBulk gpuKntSpanEtaBulk