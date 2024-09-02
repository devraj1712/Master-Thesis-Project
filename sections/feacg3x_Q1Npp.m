I2           = eye(2) ;
% Recallind order of strip in xi-direction
p_1 = orderXiCont;

% Initial and Current Surface Nodal Coordinates
XIs = Xn_Elem(1:p_1+1,:);
xIs = xn_Elem(1:p_1+1,:);

% Initialize contact force vector and stiffness matrix
rACi   = zeros(2*(p_1+1),1) ;
kACii  = zeros(2*(p_1+1), 2*(p_1+1)) ;

% Loop over Gauss Points
for ii = 1:ngps
    
    % Gauss point coordinate, weight in parent element (1D)
    xitilda = xigs(ii,1) ;
    wg      = xigs(ii,2) ;
    
    % Mapping from parent element to Parametric element
    xi               = 0.5*((kntSpanXi_Elem(2)-kntSpanXi_Elem(1))*xitilda + (kntSpanXi_Elem(2)+ kntSpanXi_Elem(1)));
    dxi_dxitilda     = 0.5*(kntSpanXi_Elem(2)-kntSpanXi_Elem(1));
    J2               = dxi_dxitilda;
    
    % Shape functions
    N                  = nrbbasisfun(xi, crvContBottom);
    dNdxi              = nrbbasisfunder(xi, crvContBottom);
    Ni                 = kron(N,I2);
    
    % Coordinates of Gauss points
    xgpi = N*xIs(:,1);
    ygpi = N*xIs(:,2);
    
    % Jacobian Transformation
    dSdxi = dNdxi*XIs(:,1);
    
    % For Straight projection, rs is
    rkP  = ygpi - yo ;
    rbk  = [0 ; -1]  ;
    
    % Force
    Fo  =    c1/rkP^9    - c2/rkP^3 ;
    dFo = -9*c1/rkP^10 + 3*c2/rkP^4 ;
    
    % Modified Force,
    if  rkP < rm
        Fo  = Fm + dFm*(rkP-rm) ;
        dFo = dFm               ;
    end
    
    % Body Force and associated Gradient
    bACi  = Fo * rbk ;
    dbACk = -dFo * (rbk)*(rbk)' ;
    
    % Residual and Stiffness
    rACi  = rACi  + Ni' * bACi       * dSdxi * wg * J2 ;
    kACii = kACii + Ni' * dbACk * Ni * dSdxi * wg * J2 ;
end

% Updating values of element stiffness matrix and force vector
iS          = (1:2*(p_1+1))      ;
Fe(iS,1)   = Fe(iS,1)  + rACi  ;
Ke(iS,iS)  = Ke(iS,iS) + kACii ;
% End of contact contribution file