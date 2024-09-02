
%% Calculation of stiffness matrices and internal force vectors in elementwise manner for contact region

% Initialize dummy stiffness matrices and internal force vector in elementwise manner 
KeDummy    = zeros(nDofCont_Elem*nDofCont_Elem, nElemXi1); % all the entries of a element stiffness matrix in a column with nElemXi_A number of columns
FeDummy    = zeros(nDofCont_Elem, nElemXi1); % all the entries of a element load vector in a column with nElemXi_A number of columns
jterminate_contact = zeros(nElemXi1, 1) ;

startTimeContactPart(ls, ns) = toc;

parfor i = 1:nElemXi1 % loop over number of contact elements

    [Ke, Fe, jterminate_Elem] = findKeFeCont(lambda, mu, D, Is4, ...
        b, ldv, ls,i,yo,c1,c2,rm,Fm,dFm, ngpve,ngps, xigve,xigs, nElemBulk+i, ...
        kntSpanXi(kntConArray(i, 1), 1:2), ...
        kntSpanEta(kntConArray(i, 2), 1:2), ...
        nDofCont_Elem, ...
        2 * cpConArrayCont(i, :), ...
        Xn(cpConArrayCont(i, :), 1:2), ...
        xn(cpConArrayCont(i, :), 1:2), ...
        crvContTop, crvContBottom, crvVertical,nxc,orderXiCont); % calculation of element stiffness matrix and internal load vector
    Ke               = Ke(:); % storing all the entries in a single column
    KeDummy(:, i)    = Ke; % storing in the dummy matrix
    FeDummy(:, i)    = Fe; % storing in the dummy vector
    jterminate_contact(i, 1) = jterminate_Elem;
end

if sum(jterminate_contact(:, 1)) > 0
    jterminate = 2;
    fprintf('\n Error in computation inside contact element loop ...\n Exiting the Newton iteration loop ...\n') ;
    error('jterminate > 0')
end

clear ie

%% Assembly of Ke and Fe to their global parts for contact part
R   = zeros(nDof, 1);
ksp = zeros(nksz, 1);
for i = 1:nElemXi1 % loop over number of contact elements
    cons          = 2 * cpConArrayCont(i, :);
    ie(1, 1:2:nDofCont_Elem) = cons - 1; % odd entries of DOF array for ith element
    ie(1, 2:2:nDofCont_Elem) = cons; % even entries of DOF array for ith element
    R(ie, 1)                 = R(ie, 1) + FeDummy(:, i); % storing entries in the global internal force vector for ith element
    ksp( (i-1) * nDofCont_Elem^2 + ...
        1 :  i * nDofCont_Elem^2, 1) ...
        = KeDummy(:, i); % storing entries in the global sparse stiffness matrix for ith element
end

clear KeDummy FeDummy jterminate_contact ie 

% computation time taken for the contact part
endTimeContactPart(ls,ns) = toc - startTimeContactPart(ls, ns);


%% Calculation of stiffness matrices and internal force vectors in elementwise manner for bulk region
time1         = toc;
% Initialize dummy stiffness matrices and internal force vector in elementwise manner 
KeDummy     = zeros(nDofBulk_Elem * nDofBulk_Elem, nElemBulk);
FeDummy     = zeros(nDofBulk_Elem, nElemBulk);
jterminate_bulk  = zeros(nElemBulk, 1);

startTimeBulkPart(ls, ns) = toc;

parfor i = 1:nElemBulk % loop over number of bulk elements

    [Ke, Fe, jterminate_Elem] = findKeFeBulk(lambda, mu, D, Is4, ...
        b, ldv, ls, ngpv, xigv, i, ...
        kntSpanXi(kntConArray(i+nElemXi1, 1), 1:2), ...
        kntSpanEta(kntConArray(i+nElemXi1, 2), 1:2), ...
        nDofBulk_Elem, ...
        2 * cpConArrayBulk(i, :), ...
        Xn(cpConArrayBulk(i, :), 1:2), ...
        xn(cpConArrayBulk(i, :), 1:2), ...
        nrbObj); % calculation of element stiffness matrix and internal load vector
    Ke               = Ke(:); % storing all the entries in a single column
    KeDummy(:, i)    = Ke; % storing in the dummy matrix
    FeDummy(:, i)    = Fe; % storing in the dummy vector
    jterminate_bulk(i, 1) = jterminate_Elem;
end

if sum(jterminate_bulk(:, 1)) > 0
    jterminate = 1;
    fprintf('\n Error in computation inside bulk element loop ...\n Exiting the Newton iteration loop ...\n') ;
    error('jterminate > 0')
end

%% Assembly of Ke and Fe to their global parts for bulk part

for i = 1:nElemBulk % loop over number of bulk elements
    cons          = 2 * cpConArrayBulk(i, :);
    ie(1, 1:2:nDofBulk_Elem) = cons - 1; % odd entries of DOF array for ith element
    ie(1, 2:2:nDofBulk_Elem) = cons; % even entries of DOF array for ith element
    R(ie, 1)                 = R(ie, 1) + FeDummy(:, i); % storing entries in the global internal force vector for ith element
    ksp( nElemXi1*nDofCont_Elem^2 + (i-1) * nDofBulk_Elem^2 + 1 : nElemXi1*nDofCont_Elem^2 +  i * nDofBulk_Elem^2, 1) ...
        = KeDummy(:, i); % storing entries in the global sparse stiffness matrix for ith element
end

clear KeDummy FeDummy jterminate_bulk ie

% computation time taken for the bulk part
endTimeBulkPart(ls, ns) = toc - startTimeBulkPart(ls, ns);

