
%% Building nodes for projected elements

nElemXi_ProjElem  = 1; % number of projected elements in a single physical element along xi
nElemEta_ProjElem = 1; % number of projected elements in a single physical element along eta
nElem_ProjElem    = nElemXi_ProjElem * nElemEta_ProjElem; % total number of projected elements in a single physical element
nElemXi_Proj      = nElemXi_ProjElem * nElemXi; % total number of projected elements for the body along xi
nElemEta_Proj     = nElemEta_ProjElem * nElemEta; % total number of projected elements for the body along eta
nElem_Proj        = nElemXi_Proj * nElemEta_Proj; % total number of projected elements for the body
nElem_Proj1       = nElem_ProjElem * nElemXi;
nKntXi_Proj       = nElemXi_Proj + 1; % total number of unique knots along xi for NURBS geometry of projected elements
nKntEta_Proj      = nElemEta_Proj + 1; % total number of unique knots along eta for NURBS geometry of projected elements


%% Knot vector in each parameter direction for projected elements

kntVecUnqXi_Proj = unique(nrbObj.knots{1});
kntVecUnqEta_Proj = unique(nrbObj.knots{2});
kntVecUnqXi_Proj = linspace(kntVecUnqXi_Proj(1,1), kntVecUnqXi_Proj(1,end), nKntXi_Proj); % unique knot vector along xi for NURBS geometry of projected elements
kntVecUnqEta_Proj = linspace(kntVecUnqEta_Proj(1,1), kntVecUnqEta_Proj(1,end), nKntEta_Proj); % unique knot vector along eta for NURBS geometry of projected elements


%% Physical coordinates array of all projected nodes

% Contact region
crvBottom = crvContBottom; % bottom or contact layer of deformable body 
crvBottom.coefs(1:2, :) = xn_storeG(1:nCPXiCont, :, zz)'; % updating the coordinates in the contact region only
crvBottom.coefs(1, :)   = crvBottom.coefs(1, :) .* crvBottom.coefs(4, :); % [x] >> [wx]
crvBottom.coefs(2, :)   = crvBottom.coefs(2, :) .* crvBottom.coefs(4, :); % [y] >> [wy]
cpCorArrayCont = nrbeval(crvBottom, {kntVecUnqXi_Proj}); % evaluation of NURBS geometry as cartesian coordinates at parameter points along one direction

% Bulk region
body    = nrbObj;
body.coefs(1:2, nCPXi+1:end) = xn_storeG(nCPXiCont+1:nCP, :, zz)'; % updating the coordinates in the bulk region
body.coefs(1, :)   = body.coefs(1, :) .* body.coefs(4, :); % [x] >> [wx]
body.coefs(2, :)   = body.coefs(2, :) .* body.coefs(4, :); % [x] >> [wx]
cpCorArrayBulk = nrbeval(body, {kntVecUnqXi_Proj, kntVecUnqEta_Proj}); % evaluation of NURBS geometry as cartesian coordinates at parameter points along two directions

% Whole body
cpCorArray_Q4Mesh = [cpCorArrayCont(1:2, :)'; cpCorArrayBulk(1:2, nKntXi_Proj+1:end)']; % coordinates array for the NURBS geometry of four-noded quadrilateral elements 


%% Knot vector connectivity array for projected elements
[kntConArray_Proj, kntSpanXi_Proj, kntSpanEta_Proj] ...
    = conKnt(nElemXi_Proj, nElemEta_Proj, kntVecUnqXi_Proj, kntVecUnqEta_Proj);


%% Control point connectivity array for projected elements
cpConArray_Q4Mesh = zeros(nElemXi_Proj * nElemEta_Proj, (1+1)*(1+1)); % (orderXi + 1) * (orderEta + 1)
for i = 1:nElemEta_Proj % loop over number of projected elements for the body along eta
    for j = 1:nElemXi_Proj % loop over number of projected elements for the body along xi
    cpConArray_Q4Mesh((nElemXi_Proj*(i-1) + j), :) = [nElemXi_Proj*(i-1)+j+(i-1), ...
            nElemXi_Proj*(i-1)+j+1+(i-1), nElemXi_Proj*(i-1)+j+1+(i-1)+nKntXi_Proj, ...
            nElemXi_Proj*(i-1)+j+(i-1)+nKntXi_Proj];
    end
end


%% Element connectivity array for projected elements
elemConArray_Q4Mesh = zeros(nElem, 1+nElemXi_ProjElem * nElemEta_ProjElem);
count     = 1; % count of physical element number

for i = 1:nElemEta
    for j = 1:nElemXi
        
        % first column define the physical element number, say 60 for 60th element
        % for elements along xi
        elemConArray_Q4Mesh(nElemXi*(i-1) + j, 1:nElemXi_ProjElem+1) ...
            = [count, nElemEta_ProjElem*(i-1)*nElemXi_Proj ...
            + (j-1)*nElemXi_ProjElem + (1:nElemXi_ProjElem)];  
        %(2*(j-1)*pnelU_a+(2*i-1)),  (2*(j-1)*pnelU_a+ 2*i),...  
        % ((2*j-1)*pnelU_a+ (2*i-1)),  ((2*j-1)*pnelU_a+ 2*i)
        
        % for elements along eta
        temp = zeros(nElemXi_ProjElem, nElemEta_ProjElem-1);
        for k = 2:nElemEta_ProjElem
            temp(:, k-1) = (nElemEta_ProjElem*(i-1)*nElemXi_Proj ...
                + (j-1)*nElemXi_ProjElem + (k-1)*nElemXi_Proj ...
                + (1:nElemXi_ProjElem))';
        end
        temp = temp(:)';
        elemConArray_Q4Mesh(nElemXi*(i-1) + j, nElemXi_ProjElem ...
            + 2: nElem_ProjElem + 1) = temp;
        count = count + 1;
        
    end
end