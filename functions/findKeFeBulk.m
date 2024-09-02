function [Ke, Fe, jterminate_Elem] = findKeFeBulk (lambda, mu, D, Is4, ...
    b, ldv, ls, ngpv, xigv, elemNo, kntSpanXi_Elem, kntSpanEta_Elem, ...
    nDofBulk_Elem, dofConArray_Elem, Xn_Elem, xn_Elem, nrbObj)
    
% Objective of the function: to calculate the element stiffness matrix and
% load vector for bulk part

% INPUT
% real      :: lambda           = Lame's constant
% real      :: mu               = Lame's constant
% real      :: D                = 
% real      :: Is4              = 
% integer   :: I2               = second-order identity matrix
% real      :: b                = body forces
% real      :: ldv              = 
% integer   :: ngpv             = number of Gauss points for the bulk volume
% real      :: xigv             = Gauss points and weights
% integer   :: elemNo           = counter of element number
% real      :: knotSpanXi_Elem  = knot span for each element along xi
% real      :: knotSpanEta_Elem = knot span for each element along eta
% integer   :: nDofBulk_Elem    = number of DOFs in an element for bulk part
% real      :: XnBulk_Elem      = coordinate data in reference
% configuration for bulk part
% real      :: xnBulk_Elem      = coordinate data in current configuration
% for bulk part
% structure :: nrbObj           = NURBS-based object

% OUTPUT
% real      :: Ke               = element stiffness matrix
% real      :: Fe               = element load vector
% real      :: jterminate       = flag for termination

jterminate_Elem   = 0;

% global DOFs
ie(1, 1:2:nDofBulk_Elem) = dofConArray_Elem - 1;
ie(1, 2:2:nDofBulk_Elem) = dofConArray_Elem;

% initialization of element stiffness matrices and load vector
Kemat   = zeros(nDofBulk_Elem, nDofBulk_Elem); % material 
Kegeo   = zeros(nDofBulk_Elem, nDofBulk_Elem); % geometric
Fe      = zeros(nDofBulk_Elem, 1);

for i = 1:ngpv % loop over number of Gauss points
   
    % coordinates and weight in master space
    xiBar    = xigv(i, 1) ;
    etaBar   = xigv(i, 2) ;
    wg       = xigv(i, 3) ;
    
    % coordinates in parameter space
    xi   = 0.5 * ((kntSpanXi_Elem(2) - kntSpanXi_Elem(1)) * xiBar + ...
             (kntSpanXi_Elem(2) + kntSpanXi_Elem(1)));
    eta  = 0.5 * ((kntSpanEta_Elem(2) - kntSpanEta_Elem(1)) * etaBar + ...
             (kntSpanEta_Elem(2) + kntSpanEta_Elem(1)));
    
    % mapping from parameter to master     
    dXi_dXiBar    = 0.5 * (kntSpanXi_Elem(2) - kntSpanXi_Elem(1));
    dEta_dEtaBar  = 0.5 * (kntSpanEta_Elem(2) - kntSpanEta_Elem(1));
    J2            = dXi_dXiBar * dEta_dEtaBar;  % Jacobian for parameter to master space
    
    % NURBS basis functions and their derivatives
    % size of the function/derivative array: [1 x (order+1)]
    [Ru]               = nrbbasisfun({xi, eta}, nrbObj); % all possible tensor products of basis functions in two directions layer-wise (row matrix), size = [1 x (orderXi+1)*(orderEta+1)]
    [dR_dXi, dR_dEta]  = nrbbasisfunder({xi, eta}, nrbObj); % all possible tensor products of function derivatives in two directions layer-wise (row matrix), size = [1 x (orderXi+1)*(orderEta+1)]
    dpR                = [dR_dXi', dR_dEta']; % array of function derivatives in two directions, size =[(orderXi+1)*(orderEta+1) x 2]
    Rb                 = kron(Ru, eye(2)); % Kronecker tensor product of Ru with I2, size = kron([1 x (orderXi+1)*(orderEta+1)] and [2 x 2]) = [2 x 2*(orderXi+1)*(orderEta+1)]
    
    % Jacobian for parameter to physical (current configuraton)
    J1     = xn_Elem' * dpR; % size = [2 x 2]
    detJ1  = det(J1);
    invJ1  = inv(J1);
    if detJ1 < 0
        jterminate_Elem = 1;
        fprintf('\n Determinant of Jacobian of element %g at time-step %g is not physical.', elemNo, ls)
    end
    
%     % Jacobian for parameter to physical (reference configuraton)
%     JRef    = XnBulk_Elem' * dpR;
%     detJRef = det(JRef);
%     invJRef = inv(JRef);
    
    % derivative of shape functions in physical space
    dR = dpR * (invJ1);
    
    % strain-displacement matrix
    B              = zeros(4, size(ie,2));
    B(1, 1:2:end)  = dR(:, 1)'; % first row of B
    B(2, 2:2:end)  = dR(:, 2)'; % second row of B
    B(4, 1:2:end)  = dR(:, 2)'; % fourth row of B
    B(4, 2:2:end)  = dR(:, 1)'; % fourth row of B
    
    % deformation gradient
    invFgr = zeros(2, 2); % invF = dX/dx = Xj*(dR/dx) = XnElem_Elem * dR
    for j = 1:(nDofBulk_Elem/2)
        invFgr = invFgr + Xn_Elem(j, :)' * dR(j, :);
    end
    Fgr    = inv(invFgr);
    detFgr = det(Fgr);
    if detFgr < 0
        jterminate_Elem = 1;
        fprintf('\n Deformation gradient of element %g at time-step %g is non-physical.', elemNo, ls)
    end
    
    % Cauchy stress (in Voigt notation and plane strain condition)  
    % sigma = (lambda/J)*(log J)*I + (mu/J)*(B-I) % Eq. 3.31
    Upr     = lambda / detFgr * log(detFgr); % (lambda/J)*(log J)
    sit     = mu / detFgr * (Fgr * Fgr' - eye(2)); % (mu/J)*(FF^T-I)
    sigma   = [Upr + sit(1,1); Upr + sit(2,2); Upr; sit(1,2)]; % [sigma_11; sigma_22; .......; sigma_12]
    
    % spatial elasticity tensor
    % c = (lambda*(I tensor product I) + 2*mu*(fourth order identity))/J -
    % (lambda/J)*(log J)*2*(fourth order identity) % Eq. 3.32
    Ds      = D / detFgr - Upr * Is4;
    
    % body forces
    bh    = ldv * b;
    
    % elemental material tangent stiffness matrix
    % int{(B^T).Ds.B}
    Kemat  = Kemat + B' * Ds * B * detJ1 * J2 * wg;
    
    % elemental geometric tangent stiffness matrix
    % int{((dR/dx)^T).(sigma).(dR/dx)}
    NsN   = dR(:,1) * dR(:,1)' * sigma(1) + dR(:,1) * dR(:,2)' * sigma(4) ...
            + dR(:,2) * dR(:,1)' * sigma(4) + dR(:,2) * dR(:,2)' * sigma(2);
    Kegeo = Kegeo + kron(NsN, eye(2)) * detJ1 * J2 * wg;
    
    % element load vector
    % int{(B^T).(sigma).(dSigma)}
    Fe = Fe + (B' * sigma - Rb' * bh / detFgr) * detJ1 * J2 * wg;
end

% element tangent stiffness matrix
Ke = Kemat + Kegeo;

end