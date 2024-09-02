% BFs = Basis functions
% CPs = Control points

%% Geometric and parametric details

Lx      = 200; % lenght of the body
Ly      = 10; % height of the body
nxc      = 150 ;
shiftX  = 0; % x-coordinate shift from the reference
shiftY  = 0; % y-coordinate shift from the reference

% Initial order of BFs
orderXi      = 1; % order of BFs along xi
orderEta     = 1; % order of BFs along eta
orderXiCont  = orderXi; % for later use in contact surface construction
orderEtaCont = orderEta; % for later use in contact surface construction

% Initial number of elements
resXi    = 200; % minimum number of elements along xi
resEta   = 10; % minimum number of elements along eta
nElemXi  = resXi * meshNo; % number of elements along xi
nElemEta = resEta * meshNo; % number of elements along xi
nElemXi1 = nElemXi;
nElemEta1 = nElemEta;

% Initial knot vector
kntVecXi      = [zeros(1, orderXi+1), Lx * ones(1, orderXi+1)]; % knot vector along xi
kntVecEta     = [zeros(1, orderEta+1), Ly * ones(1, orderEta+1)]; % knot vector along eta
kntVecXiCont  = kntVecXi; % for later use in contact surface construction
kntVecEtaCont = kntVecEta; % for later use in contact surface construction

% Initial control points array
cpArray(:, :, 1) = [0, Lx; 0, 0; 0, 0; 1, 1]; % bottom layer CPs
cpArray(:, :, 2) = [0, Lx; Ly, Ly; 0, 0; 1, 1]; % top layer CPs

% Control point shift, if necessary
cpArray(1, :) = cpArray(1, :) + shiftX; % updating x-coordinates
cpArray(2, :) = cpArray(2, :) + shiftY; % updating y-coordinates

% Converting control point coordinates into its weighted form
cpArray(1, :) = cpArray(1, :) .* cpArray(4, :); % [x] --> [wx]
cpArray(2, :) = cpArray(2, :) .* cpArray(4, :); % [y] --> [wy]


%% Geometry construction and discretization

% Geometry construction
nrbObj = nrbmak(cpArray, {kntVecXi, kntVecEta}); % initial geometry construction

% Order elevation
orderElevXi  = str2double(bulkType(1,2)) - orderXi; % elevation of order along xi
orderElevEta = str2double(bulkType(1,4)) - orderEta; % elevation of order along eta
nrbObj       = nrbdegelev(nrbObj, [orderElevXi, orderElevEta]); % reconstruction of geometry

% Knot insertion
kntVecInsXi  = setdiff(linspace(0, Lx, nElemXi+1), [0, Lx]); % inserted knots along xi
kntVecInsEta = setdiff(linspace(0, Ly, nElemEta+1), [0, Ly]); % inserted knots along eta
nrbObj       = nrbkntins(nrbObj, {kntVecInsXi, kntVecInsEta}); % reconstruction of geometry
nrbObjCom    = nrbObj; 
% nrbPlot(nrbObj);

% Update the parametric details
orderXi   = nrbObj.order(1) - 1; % order of basis functions along xi
orderEta  = nrbObj.order(2) - 1; % order of basis functions along eta
kntVecXi  = nrbObj.knots{1}; % knot vector along xi
kntVecEta = nrbObj.knots{2}; % knot vector along eta
nCPXi     = nrbObj.number(1); % number of control points along xi
nCPEta    = nrbObj.number(2); % number of control points along eta
nElemXi   = size(unique(nrbObj.knots{1}), 2) - 1; % number of elements along xi 
nElemEta  = size(unique(nrbObj.knots{2}), 2) - 1; % number of elements along eta
cpArray   = nrbObj.coefs; % control points array


%% Construction and discretization of vertical surface and top and bottom part of the contact elements

% Vertical curve
cpArrayCrvVertical       = [0, 0; 0, Ly; 0, 0; 1, 1]; 
cpArrayCrvVertical(1, :) = cpArrayCrvVertical(1, :) + shiftX; 
cpArrayCrvVertical(2, :) = cpArrayCrvVertical(2, :) + shiftY; 
crvVertical = nrbmak(cpArrayCrvVertical, kntVecEtaCont); 
crvVertical = nrbdegelev(crvVertical, str2double(bulkType(1,4)) - orderEtaCont); 
crvVertical = nrbkntins(crvVertical, kntVecInsEta); 
% nrbPlot(crvVertical); % plotting of vertical curve

% Top part of the contact region
crvContTop = nrbmak(cpArray(:,:,2), kntVecXi);
% nrbPlot(crvBottomXiElem); % plotting of bottom part

% Bottom part of the contact region
cpArrayCrvContBottom(:, :, 1) = [0, Lx; 0, 0; 0, 0; 1, 1]; % top layer
cpArrayCrvContBottom(1, :)    = cpArrayCrvContBottom(1, :) + shiftX; % updating x-coordinates
cpArrayCrvContBottom(2, :)    = cpArrayCrvContBottom(2, :) + shiftY; % updating y-coordinates
cpArrayCrvContBottom(1, :)    = cpArrayCrvContBottom(1, :) .* cpArrayCrvContBottom(4, :); % [x] --> [wx]
cpArrayCrvContBottom(2, :)    = cpArrayCrvContBottom(2, :) .* cpArrayCrvContBottom(4, :); % [y] --> [wy]
crvContBottom = nrbmak(cpArrayCrvContBottom, kntVecXiCont); % construction of Bottom part
crvContBottom = nrbdegelev(crvContBottom, str2double(elemType(1,2)) - orderXiCont); % reconstruction of Bottom part using order elevation
crvContBottom = nrbkntins(crvContBottom, kntVecInsXi); % reconstruction of Bottom part using knot insertion
orderXiContTemp  = crvContBottom.order - 1; % temporary order
kntVecXiContTemp = crvContBottom.knots; % temporary knot vector
nCPXiContTemp    = crvContBottom.number; % temporary number of control points
crvContBottom       = nrbdegelev(crvContBottom, str2double(elemType(1,4))); % reconstruction of Bottom part using order elevation
% nrbPlot(crvTopXiElem); % plotting of top part

% Update the parametric details
orderXiCont       = crvContBottom.order - 1; % order of BFs for Bottom part along xi
kntVecXiCont      = crvContBottom.knots; % knot vector for Bottom part along xi
nCPXiCont        = crvContBottom.number; % number of CPs for Bottom part along xi
cpArrayCrvContBottom = crvContBottom.coefs; % control points array for Bottom part

%Top curve of the top surface
tucrv = nrbmak(cpArray(:,:,end), kntVecXi);

%% Geometrical details for overall domain (standard + VO)

% Converting weighted form of control point coordinates into its original form
nrbObj.coefs(1, :) = nrbObj.coefs(1, :) ./ nrbObj.coefs(4, :); % [wx] --> [x] for bulk part
nrbObj.coefs(2, :) = nrbObj.coefs(2, :) ./ nrbObj.coefs(4, :); % [wy] --> [y] for bulk part
crvContBottom.coefs(1, :) = crvContBottom.coefs(1, :) ./ crvContBottom.coefs(4, :); % [wx] --> [x] for contact part
crvContBottom.coefs(2, :) = crvContBottom.coefs(2, :) ./ crvContBottom.coefs(4, :); % [wy] --> [y] for contact part
crvVertical.coefs(1, :) = crvVertical.coefs(1, :) ./ crvVertical.coefs(4, :); % [wx] --> [x]
crvVertical.coefs(2, :) = crvVertical.coefs(2, :) ./ crvVertical.coefs(4, :); % [wy] --> [y]

% Storing x- and y-coordinates of contact part

XnBulk = nrbObj.coefs(1:2, 1+nCPXi:end)'; % bulk part
XnCont = crvContBottom.coefs(1:2, :)'; % contact part
Xn = [XnCont; XnBulk]; % whole geometry
xn = Xn; % Storing reference into spatial coordinates array for later use

%%  Control points connectivity arrays

nCP           = nCPXi * nCPEta - nCPXi + nCPXiCont; % total number of CPs
nCPBulk       = nCP - nCPXiCont; % number of CPs for bulk part
nDof         = 2 * nCP; % total degrees of freedom
nDofCont_Elem = 2 * ((orderXiCont+1) + (orderXi+1) * orderEta); % degrees of freedom per element for contact part
nDofBulk_Elem = 2 * (orderXi + 1) * (orderEta + 1); % degrees of freedom per element for bulk part
nElem         = nElemXi * nElemEta; % total number of elements
nElemCont     = nElemXi; % number of contact elements
nElemBulk     = nElem - nElemXi1; % number of bulk elements

% Array for the whole domain
[cpConArray, nrbsCorArray] = conCP(orderXi, orderEta, nCPXi, nCPEta); 
cpConArray                   = cpConArray';

% Array for the contact curve in contact region
cpConArrayCont1  = 1:(orderXiCont+1);

if any(strcmp(elemType, {'N1.0', 'N2.0', 'N3.0', 'N4.0', 'N5.0', 'N6.0', 'N8.0'})) == 1
    for i = 1:nElemCont-1
        temp = 1;
        cpConArrayCont1(1+i,:) = cpConArrayCont1(i,:) + temp;
    end
elseif any(strcmp(elemType, {'N1.1', 'N2.1', 'N3.1', 'N4.1', 'N5.1'})) == 1
    for i = 1:nElemCont-1
        temp = 2;
        cpConArrayCont1(1+i,:) = cpConArrayCont1(i,:) + temp;
    end
elseif any(strcmp(elemType, {'N1.2', 'N2.2', 'N3.2', 'N4.2', 'N5.2'})) == 1
    for i = 1:nElemCont-1
        temp = 3;
        cpConArrayCont1(1+i,:) = cpConArrayCont1(i,:) + temp;
    end
elseif any(strcmp(elemType, {'N1.3', 'N2.3', 'N3.3', 'N4.3', 'N5.3'})) == 1
    for i = 1:nElemCont-1
        temp = 4;
        cpConArrayCont1(1+i,:) = cpConArrayCont1(i,:) + temp;
    end
elseif any(strcmp(elemType, {'N1.4', 'N2.4', 'N3.4', 'N4.4', 'N5.4'})) == 1
    for i = 1:nElemCont-1
        temp = 5;
        cpConArrayCont1(1+i,:) = cpConArrayCont1(i,:) + temp;
    end       
elseif any(strcmp(elemType, {'N1.5', 'N2.5', 'N3.5', 'N4.5', 'N5.5'})) == 1
    for i = 1:nElemCont-1
        temp = 6;
        cpConArrayCont1(1+i,:) = cpConArrayCont1(i,:) + temp;
    end 
end

% VO-based updated array
excessCP = cpConArrayCont1(nElemXi, end) - cpConArray(nElemXi, orderXi+1); % excess control points for VO discretization
cpConArray(1:nElemXi, (orderXi+1)+1:end) = cpConArray(1:nElemXi, (orderXi+1)+1:end) + excessCP; % updating global number for bulk CPs participating in contact elements
cpConArray(nElemXi+1:end, :) = cpConArray(nElemXi+1:end, :) + excessCP; % updating global number of CPs in the bulk region

% Array for the contact region
cpConArrayCont = zeros(nElemCont, ((orderXiCont + 1) + (orderXi + 1) * orderEta));
for i = 1:nElemCont
    cpConArrayCont(i, :) = [cpConArrayCont1(i, :), cpConArray(i, (orderXi + 1) + 1:end)];
end
cpConArrayCont2 = cpConArrayCont1;

% Array for the bulk region
cpConArrayBulk = cpConArray(nElemCont+1:end, :);

%% Knot locations connectivity arrays

[kntConArray, kntSpanXi, kntSpanEta] = conKnt(nElemXi1, nElemEta1, kntVecXi, kntVecEta);

