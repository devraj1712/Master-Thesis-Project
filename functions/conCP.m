function [cpConArray, nrbsCorArray] = conCP(orderXi, orderEta, nCPXi,nCPEta)

% Objective of the function: to find the control point connectivity between
% the elements as an array

% INPUT
% integer :: orderXi       = order of basis functions along xi
% integer :: orderEta      = order of basis functions along eta    
% integer :: nCPXi         = number of control points along xi     
% integer :: nCPEta        = number of control points along eta    

% OUTPUT
% integer :: cpConArray    = control point connectivity array      
% integer :: nrbsCorArray  = NURBS coordinates array               

    nElemXi         = nCPXi - orderXi;                % number of elements along xi
    nElemEta        = nCPEta - orderEta;              % number of elements along eta
    nElem           = nElemXi * nElemEta;             % total number of elements
    nCP_Elem        = (orderXi + 1) * (orderEta + 1); % number of control points in an element
    nCP             = nCPXi * nCPEta;                 % total number of control points
    cpConArray      = zeros(nCP_Elem, nElem); 
    nrbsCorArray    = zeros(nCP, 2);
    A               = 1;                        % counter variable
    e               = 1;                        % counter variable
    for i = 1:nCPEta                            % loop over number of control points along xi
        for j = 1:nCPXi                         % loop over number of control points along eta
            nrbsCorArray(A, 1)           = j;
            nrbsCorArray(A, 2)           = i;
            if j >= (orderXi + 1) && i >= (orderEta + 1)
                for ii = 0:orderEta     % loop over order of basis functions along xi
                    for jj = 0:orderXi  % loop over order of basis functions along eta
                        B                = A - ii * nCPXi - jj;          % global function number
                        b                = ii * (orderXi + 1) + jj + 1;  % local function number
                        cpConArray(b, e) = B;
                    end
                end
                e = e + 1; % increment of counter variable by unity
            end
            A = A + 1;     % increment of counter variable by unity
        end
    end
    cpConArray = sort(cpConArray); % sorting in ascending order
    
end