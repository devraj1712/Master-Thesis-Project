% Initialization of Newton-Iteraton loop
eN       = 1;  % initialization of energy residual value.
eN_mat   = [];
residual = []; % storing all the values of residual for a time-step

UgOld  = Ug;
% VGold  = VG;
% AGold  = AG;
ns     = 0;
ns_mid = 0;

% if flag == 0
itrcon_flag        = 0;
startTimeNR(ls, 1) = toc;
% end

while eN > TOL
    
    ns = ns + 1; % Newton step
    % Updating displacement vector
    Ug(iux) = ubx + uprx + Ug(ixo)*ones(siux,1) ;
    Ug(iuz) = ubz + uprz + Ug(izo)*ones(siuz,1) ;
    % Current configuration
    Ux = [ Ug(1:2:nDof) Ug(2:2:nDof) ] ;
    xn = Xn + Ux(:,1:2)  ;        % New co-ordinate position of control points
    clear Ux;
    startTimeBulkPart(ls, ns) = toc;
    elemGPUCompute; 
    endTimeBulkPart(ls, ns) = toc - startTimeBulkPart(ls, ns);
    
    % Terminate computation
    if jterminate == 1 || jterminate == 2 || jterminate == 3 || jterminate == 4
        disp('jterm ~= 0, programm has been terminated.');
    end
    
    % Global sparse stiffness matrix
    K   = sparse(isp, jsp, ksp); % storing entries of ksp in a sparse matrix
    %clear ksp
    
    % External load vector (if any)
    Fext = zeros(nDof,1);
    if bcType == 1 % Neumann BC
        ext_force;
    end
    
    % Global residual load vector
    % f(internal) - f(external) + f(contact)
    f = R - ldv * F; % f(internal) - f(external)
    f = f - Fext; % f(internal) - f(external from BC)
    
    % BC modification for Center Node,
        f(ixo)     = f(ixo) + sum(f(iux)) ;
        f(izo)     = f(izo) + sum(f(iuz)) ;
        
        K(ir,ixo)  = K(ir,ixo) + K(ir,iux)*ones(siux,1) ;
        K(ir,izo)  = K(ir,izo) + K(ir,iuz)*ones(siuz,1) ;
        
        K(ixo,ir)  = K(ixo,ir) + ones(1,siux)*K(iux,ir) ;
        K(izo,ir)  = K(izo,ir) + ones(1,siuz)*K(iuz,ir) ;
        
        K(ixo,ixo) = K(ixo,ixo) + ones(1,siux)*K(iux,iux)*ones(siux,1) ;
        K(ixo,izo) = K(ixo,izo) + ones(1,siux)*K(iux,iuz)*ones(siuz,1) ;
        K(izo,ixo) = K(izo,ixo) + ones(1,siuz)*K(iuz,iux)*ones(siux,1) ;
        K(izo,izo) = K(izo,izo) + ones(1,siuz)*K(iuz,iuz)*ones(siuz,1) ;
        
    
    startTimeSolver(ls, ns) = toc;
    % Compute displacement field
    Kgr       = K(ir, ir); % reduced global stiffness matrix
    Fgr       = f(ir, 1); % reduced global load vector
    Ugr       = -Kgr\Fgr; % reduced displacement vector
    DU        = zeros(nDof, 1); % displacement vector
    DU(ir,1)  = Ugr; % storing reduced displacement into main displacement vector 
    Ug        = Ug + DU; % update the displacement field
    endTimeSolver(ls, ns) = toc - startTimeSolver(ls, ns);    
    
    % Current coordinates of body 
    
    Ux = [Ug(1:2:nDof) Ug(2:2:nDof) ]; % storing U in required format
    xn = Xn + Ux(:, 1:2); % current coordinate vector
%     clear Ux_A Ux_B;
    
    % Convergence check (energy-residual method)
    eN        = abs(Fgr' * Ugr) / nir; % energy residual per DOF at a N-R iteration
    residual  = [residual eN] ;
    eN_mat    = [eN_mat, log10(eN) ];
    fprintf('\nEnergy Residual at N-R iteration %i: %.6f \n', ns, log10(eN));
    
    % Continue only if ns > 40 and log10(eN) >= -18
    nlimc = nlimc_max; % nlimc_max is set at that start of program
    if (ns - ns_mid > nlimc) && (log10(eN) < eN_slowConvergence)
        eN = 0;
    end
    
    % Terminate if ns > 300 or log10(eN) > 3, i.e., NR didn't converve
    nlim = nlim_max; % nlim_max is set at that start of program
    if (ns - ns_mid > nlim) || log10(eN) > 3
        error('ns > nlim_max');
    end
    
    % Update control point coordinates of NURBS geometry for next N-R iteration
    crvContBottom.coefs(1:2, :)  = xn(1: nCPXiCont , :)'; 
    crvContTop.coefs(1:2, :)     = xn(nCPXiCont + 1 : nCPXiCont+nCPXi , :)'; 
    tucrv.coefs(1:2, :) = xn(end - nCPXi + 1:end, :)';
    clear Kg Fg DU UX fr  % clear KGr FGr;
end                   

% plotting deformed configuration
jplot_body  = 0;           % to plot deformed/reference configuration
jplot_curve = 0;           % to plot contact curve
    if jplot_body == 1  % for 2D geometries
        nrbObj1                = nrbObj;
        nrbObj1.coefs(1,1:nCPXiCont, 1)  = xn(1:nCPXiCont,1)';
        nrbObj1.coefs(2,1:nCPXiCont, 1)  = xn(1:nCPXiCont,2)';
        nrbObj1.coefs(1:2, 2:end)  = xn(1+nCPXiCont:end)';
        nrbObj1.coefs(1, :)    = nrbObj1.coefs(1, :) .* nrbObj1.coefs(4, :);
        nrbObj1.coefs(2, :)    = nrbObj1.coefs(2, :) .* nrbObj1.coefs(4, :);
        
    elseif jplot_body == 2  % for reference configuration
        nrbObj1          = nrbObj;
        nrbObj1.coefs(1:2, :)  = Xn';
        nrbObj1.coefs(1, :)    = nrbObj1.coefs(1, :) .* nrbObj1.coefs(4, :);
        nrbObj1.coefs(2, :)    = nrbObj1.coefs(2, :) .* nrbObj1.coefs(4, :);
        
%         figure('units', 'normalized', 'position', [0 0 1 1]); 
%         nrbkntplot(nrbObj1_B); hold on;
%         nrbkntplot(nrbObj1_A);
%         axis equal; axis off; view(2);
%         set(gcf, 'color', 'w');
    end

    if jplot_curve == 1 % for contact curves
        crvContBottom1                = crvContBottom;
        crvContBottom1.coefs(1, :)    = crvContBottom1.coefs(1, :) ...
                                            .* crvContBottom1.coefs(4,:);
        crvContBottom1.coefs(2, :)    = crvContBottom1.coefs(2, :) ...
                                            .* crvContBottom1.coefs(4,:);
        crvContTop1                   = crvContTop;
        crvContTop1.coefs(1, :)       = crvContTop1.coefs(1, :) ...
                                            .* crvContTop1.coefs(4,:);
        crvContTop1.coefs(2, :)       = crvContTop1.coefs(2, :) ...
                                            .* crvContTop1.coefs(4,:);
        crvVertical1             = crvVertical;
        
        crvVertical1.coefs(1, :) = crvVertical1.coefs(1, :) ...
                                            .* crvVertical1.coefs(4,:);
        crvVertical1.coefs(2, :) = crvVertical1.coefs(2, :) ...
                                            .* crvVertical1.coefs(4,:);
%         figure('units', 'normalized', 'position', [0 0 1 1], 'visible', 'on'); 
%         nrbkntplot(crvContTop1_B); hold on;
%         nrbkntplot(crvContBottom1_A); hold on;
%         nrbkntplot(crvVertical1_A); hold on;
%         axis auto; grid on; viewset(gcf, 'color', 'w'); view(2); hold on;
    elseif jplot_curve == 2 % Contact curve in the deformed configuration
        crvContBottom1  = crvContBottom;
        crvContTop1     = crvContTop;
        
        crvContBottom1.coefs(1:2, :) = xn(1:nCPXiCont, :)';
        crvContBottom1.coefs(1, :)   = crvContBottom1.coefs(1, :) ...
                                        .* crvContBottom1.coefs(4, :);
        crvContBottom1.coefs(2, :)   = crvContBottom1.coefs(2, :)...
                                        .* crvContBottom1.coefs(4, :);
        crvContTop1.coefs(1:2, :) = xn((cpConArrayCont(nElemXi, ...
                    (orderXiCont + 1) + (orderXi + 1) * orderEta) ...
                    - nCPXi + 1:cpConArrayCont(nElemXi, (orderXiCont + 1) ...
                    + (orderXi + 1) *  orderEta)), :)';
        crvContTop1.coefs(1, :)   = crvContTop1.coefs(1, :) ...
                                        .* crvContTop1.coefs(4, :);
        crvContTop1.coefs(2, :)   = crvContTop1.coefs(2, :) ...
                                        .* crvContTop1.coefs(4, :);
        crvContTop1.coefs(1:2, :) = xn(end - nCPXiCont + 1 : end, :)';
        
%         figure('units', 'normalized', 'position', [0 0 1 1], 'visible', 'on'); 
%         nrbkntplot(crvContTop1_B); hold on;
%         nrbkntplot(crvContBottom1_A); hold on;
%         axis auto; grid on; viewset(gcf, 'color', 'w'); view(2); hold on;
    end


% Terminate computation
if jterminate == 1 || jterminate == 2 || jterminate == 3 || jterminate == 4
    fprintf('\nProgram gets terminated as jterminate is not equal to 0 \n');
end

% % Newton steps
% if ns - ns_mid > nsm && jterminate == 0
%     nsm = ns - ns_mid; % maximum number of iteration steps
% end

% Declaring N-R converged
if eN < TOL  && jterminate == 0 % && flag == 0
    itrcon_flag = 1;
    endTimeNR(ls, 1) = toc - startTimeNR(ls, 1);
%     fprintf('\n Time taken by N-R iterations is %g \n', endTimeNR(ls, 1));
end

% Load value as N-R iterations converge
if itrcon_flag == 1 &&  jterminate == 0
    fprintf('\nTotal N-R iterations for Load-step %i: %i \n', ls - 1, ns);
    % fprintf(fipf,'\n Load step %2i used %2i Newton-steps \n\n',l-1+e0,ns-ns_mid);
%     if flag_num_scheme == 2 && flag == 0
%         ns_mid = ns;
%     end
    % Load under Displacement Drive
    if ls == 1
        Reas    = 0 ;
        Reas3   = 0 ;
        Moment1 = 0 ;             % Initial values of load and moment
    end
    
    ReasOld   = Reas ;
    MomentOld = Moment1 ;
    Reac      = zeros(siux,1); Reac2 = zeros(siux,1);
    
    FN = zeros(size(iuz,1),1);
    FT = zeros(size(iux,1),1);
        for i = 1:siux
            distancex = xn(crvRight(i),1)-xn(centerNode,1) ;
            distancez = xn(crvRight(i),2)-xn(centerNode,2) ;
            Reac(i)   = R(iuz(i))*distancex - R(iux(i))*distancez ;
            Reac2(i)  = Reac(i) ;
            FN(i)     = R(iuz(i));
            FT(i)     = R(iux(i));
            if i == 1
                Reas2 = 0 ;
            end
            Reas2 = Reas2 + sqrt((Reac(i))^2 + (Reac2(i))^2) ;
        end
    Reas  = sum(Reac) ;
    Reas3 = Reas ;
    Moment1 = Reas ;
end 
