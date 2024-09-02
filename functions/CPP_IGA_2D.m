function [xiP, kntSpanP, xIeP, itr, XmP, nP, aP, gN, a11, k11] = ...
    CPP_IGA_2D(xgp, xiPo, xn_B, cpConArrayCont2_B, ...
    orderXiContTemp_B, nCPXiContTemp_B, kntVecXiContTemp_B, ...
    orderXiCont_B, kntVecXiCont_B, crvContTop_B)

% Objective of the function: to find the closest projection point

% INPUT
% float     :: xgp             = coordinate corresponding to GP
% integer   :: orderXiElem     = order of basis functions representing contact curve
% float     :: xiPo            = old projection point that is used as initial value for N-R iterations
% float     :: kntVecXiElem    = knot vector for the contact curve
% float     :: xn              = coordinate of the deformed configuration
% structure :: crvTopXiElem    = contact curve
% integer   :: cpConArrayNew   = control point connectivity array of contact curve
% integer   :: noCnpnue_b_temp = number of control points along the contact curve
% integer   :: pe_b_temp       = order of basis functions representing contact curve
% float     :: uKnote_b_temp   = knot vector for the contact curve

% OUTPUT
% float   :: xiP      = parametric coordinate of the closest point projected on contact curve
% float   :: kntSpanP = knot span of the contact element where xiP is located
% float   :: xIeP     = coordinates of the contact element where xiP is located
% integer :: itr      = total number of iterations used to solve Fu =
% ||(xgp - C(xiP))|| --> minimum

% Initialize NR itr count
itr        = 0; % N-R iteration count
cppConv    = 0;

while cppConv == 0
    
    itr = itr + 1; % update iteration
    
    if itr == 1
        % knot span corresponding to given projection point
        kntSpanP = findspan(nCPXiContTemp_B - 1, orderXiContTemp_B, xiPo, kntVecXiContTemp_B);
        kntSpanP = kntSpanP - (orderXiContTemp_B - 1);
        
        % find coordinates of CPs corresponding to kntSpanP
        ie_B = cpConArrayCont2_B(kntSpanP, :);
        xIeP = xn_B(ie_B,:);
    end
    
    % Point, tangent, and derivative of tangent at xiPo
    [NP, dNP_dXi, d2NP_dXi2] = NURBS1DBasis2ndDers(xiPo, orderXiCont_B, ...
                crvContTop_B.knots, crvContTop_B.coefs(4, :)');
    Cu      = NP * xIeP;         % x-coordinate of projection point i.e. C(u)
    CuDash  = dNP_dXi * xIeP ;   % first order derivative of C(u)
    CudDash = d2NP_dXi2 * xIeP ; % second order derivative of C(u)
    
    % Update new projection point xiP
    PCu      = xgp - Cu;
    Fud      = PCu * CuDash';
    FudDash  = PCu * CudDash' - (CuDash * CuDash');
    xiP      = xiPo - (Fud / FudDash); % new parametric point on C(u)
    
    % Check if out of current segment
    if xiP <  kntVecXiCont_B(1, 1)
        xiP = kntVecXiCont_B(1, 1);
    elseif xiP > kntVecXiCont_B(1, end)
        xiP = kntVecXiCont_B(1, end);
    end

    % Find the corresponding knot span or element where xiP is located
    kntSpanP      = findspan(nCPXiContTemp_B - 1, orderXiContTemp_B, xiP, kntVecXiContTemp_B);
    kntSpanP      = kntSpanP - (orderXiContTemp_B - 1);
    ie_B          = cpConArrayCont2_B(kntSpanP, :);
    xIeP          = xn_B(ie_B, :);
    [NP, dNP_dXi] = NURBS1DBasisDers(xiP, orderXiCont_B, crvContTop_B.knots, crvContTop_B.coefs(4,:)');
    XmP           = NP * xIeP; % coordinates of xiP in deformed configuration
    CuDash        = dNP_dXi * xIeP; % covariant tangent vector
    
    % Convergence criteria
    RelErr  = abs(xiP - xiPo);         % criteria 1: error b/w new and old projection points
    dFu     = abs((xgp - XmP) * CuDash'); % criteria 2
    
    if itr > 500  % this is just to break the while loop
        disp('CPP: itr > 500')
        break  % to avoid  for non-contacting xgpi points
    end
    
    % Update the computed point
    xiPo   = xiP;
    
    % Declare convergence
    if (RelErr > 1e-14) && dFu > 10 * eps
        cppConv = 0;
    else
        cppConv = 1;
    end
end  

% Check if computed point is outside the knot span
if xiP < kntVecXiCont_B(1, 1)
    xiP = kntVecXiCont_B(1, 1);
elseif xiP > kntVecXiCont_B(1, end)
    xiP = kntVecXiCont_B(1, end);
end

% Compute quantites at xiP
[~, ~, d2NP_dXi2] = NURBS1DBasis2ndDers(xiP, orderXiCont_B, ...
                    crvContTop_B.knots, crvContTop_B.coefs(4,:)');
aP  = CuDash;             % covariant tangent vector at xiP
daP = d2NP_dXi2 * xIeP;   % derivative of aP at xiP
%aP = aP';                % aP as nP
nP  = [-aP(2); aP(1)] / norm(aP); % unit normal vector at xiP
gN  = (xgp - XmP) * nP;            % normal gap
a11 = aP * aP';                   % metric tensor
k11 = daP * nP;                   % Curvature tensor

end 