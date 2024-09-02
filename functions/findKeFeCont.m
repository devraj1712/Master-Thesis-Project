function [Ke, Fe, jterminate_Elem] = findKeFeCont(lambda, mu, D, Is4, ...
    b, ldv, ls,i,yo,c1,c2,rm,Fm,dFm, ngpv,ngps, xigv, xigs,elemNo, kntSpanXi_Elem, kntSpanEta_Elem, ...
    nDofCont_Elem, dofConArray_Elem, Xn_Elem, xn_Elem, crvContTop, ...
    crvContBottom, crvVertical,nxc,orderXiCont)

% Objective of the function: to calculate the element stiffness matrix and
% load vector for contact part

% INPUT
% real      :: lambda              = Lame's constant
% real      :: mu                  = Lame's constant
% real      :: D                   = 
% real      :: Is4                 = 
% integer   :: I2                  = second-order identity matrix
% real      :: b                   = body forces
% real      :: ldv                 = 
% integer   :: ngpv                = number of Gauss points for the bulk volume
% real      :: xigv                = Gauss points and weights
% integer   :: elemNo              = counter of element number
% real      :: knotSpanXi_Elem     = knot span for each element along xi
% real      :: knotSpanEta_Elem    = knot span for each element along eta
% integer   :: nDofElem_Elem       = number of DOFs in an element for contact
% part
% real      :: XnElem_Elem         = coordinate data in reference
% configuration for contact part
% real      :: xnElem_Elem         = coordinate data in current configuration
% for contact part
% structure :: crvBottomXiElem     = bottom part of the contact elements
% structure :: crvTopXiElem        = top part of the contact elements
% structure :: crvVerticalEtaElem  = vertical surface

% OUTPUT
% real      :: Ke               = element stiffness matrix
% real      :: Fe               = element load vector
% real      :: jterminate       = flag for termination

jterminate_Elem   = 0;

% global DOFs
ie(1, 1:2:nDofCont_Elem) = dofConArray_Elem - 1;
ie(1, 2:2:nDofCont_Elem) = dofConArray_Elem;

% initialization of element stiffness matrices and internal load vector
Kemat    = zeros(nDofCont_Elem, nDofCont_Elem); % material 
Kegeo    = zeros(nDofCont_Elem, nDofCont_Elem); % geometric
Fe       = zeros(nDofCont_Elem, 1);

for gp = 1:ngpv % loop over number of Gauss points
    
    % coordinates and weight in master space
    xiBar    = xigv(gp, 1); 
    etaBar   = xigv(gp, 2); 
    wg       = xigv(gp, 3); 
    
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
    [N1]       = nrbbasisfun(xi, crvContBottom); % function values at the bottom layer of contact region, size = [1 x (orderXi+1)]
    [N2]       = nrbbasisfun(xi, crvContTop); % function values at the top layer of contact region, size = [1 x (orderXi+1)]
    [M]        = nrbbasisfun(eta, crvVertical); % function values at the vertical layer, size = [1 x (orderEta+1)]
    [dN1_dXi]  = nrbbasisfunder(xi, crvContBottom); % function derivatives at the bottom layer of contact region, size = [1 x (orderXi+1)]
    [dN2_dXi]  = nrbbasisfunder(xi, crvContTop); % function derivatives at the top layer of contact region, size = [1 x (orderXi+1)]
    [dM_dEta]  = nrbbasisfunder(eta, crvVertical); % function derivatives at the vertical layer, size = [1 x (orderEta+1)]
    R_1        = N1 * M(1, 1); % tensor product of N1 and M to take care first layer forming contact region, size = [1 x (orderXi+1)] 
    R_2        = N2 * M(1, 2); % tensor product of N2 and M to take care rest layers forming bulk region, size = [(orderXi+1) x orderEta]
    Ru         = [R_1, R_2]; % all possible tensor products of basis functions in two directions layer-wise (row matrix), size = [1 x ((orderXi+1)+(orderXi+1)*orderEta)] = [1 x (orderXi+1)*(orderEta+1)]
    dN2xi      = dN1_dXi' * M(1, 1); % tensor product of dN2_dXi and M, size = [(orderXi+1) x orderEta]
    dM2eta     =  dM_dEta(1, 1)*N1'; % tensor product of dN2 and dM_dEta, size = [(orderXi+1) x orderEta]
    dpR(:, 1)  = [dN2xi ; dN2_dXi' * M(1,2)]; % first column of array of function derivatives in two directions, size =[(orderXi+1)*(orderEta+1) x 1]
    dpR(:, 2)  = [dM2eta ;  dM_dEta(1,2)*N2' ]; % second column of array of function derivatives in two directions, size =[(orderXi+1)*(orderEta+1) x 1]
    Rb         = kron(Ru, eye(2)); % Kronecker tensor product of Ru with I2, size = kron([1 x (orderXi+1)*(orderEta+1)] and [2 x 2]) = [2 x 2*(orderXi+1)*(orderEta+1)]
    
    % Jacobian for parameter to physical (current configuraton)
    J1     = xn_Elem' * dpR; % size = [2 x 2]
    detJ1  = det(J1);
    invJ1  = inv(J1);
    if detJ1 < 0
        jterminate_Elem = 1;
        fprintf('\n Determinant of Jacobian of element %g at time-step %g is not physical.', elemNo, ls)
    end
    
%     % Jacobian for parameter to physical (reference configuraton)
%     JRef    = XnElem_Elem' * dpR;
%     detJRef = det(JRef);
%     invJRef = inv(JRef);
    
    % derivative of shape functions in physical space
    dR = dpR * (invJ1); % size (dR) = size(dpR)
    
    % strain-displacement matrix
    B              = zeros(4, size(ie,2));
    B(1, 1:2:end)  = dR(:, 1)'; % first row of B
    B(2, 2:2:end)  = dR(:, 2)'; % second row of B
    B(4, 1:2:end)  = dR(:, 2)'; % fourth row of B
    B(4, 2:2:end)  = dR(:, 1)'; % fourth row of B
    
    % deformation gradient
    invFgr = zeros(2, 2); % invF = dX/dx = Xj*(dR/dx) = XnElem_Elem * dR
    for j = 1:(nDofCont_Elem/2)
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
    bh  = ldv * b;
    
    % elemental material tangent stiffness matrix
    % int{(B^T).Ds.B}
    Kemat  = Kemat +  B' * Ds * B * detJ1 * J2 * wg;
    
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
% Compute the contact contribution
if i <= nxc
    time_start_contact  = cputime  ;
    feacg3x_Q1Npp;
    total_time_contact_inc_ele = cputime - time_start_contact ;
end
end