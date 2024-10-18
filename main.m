close all; clear; clc;

path1 = genpath('nurbs-1.4.3');
path2 = genpath('igafem-v2');
path3 = genpath('functions');
path4 = genpath('sections');
addpath(path1,path2,path3,path4)  

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %  
%                           PRE-PROCESSING                                %  
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars 

clc; tic;

bcType     = 2; % for boundary conditions
% jsave      = 1; % for saving the data
jterminate = 0; 
% jstore     = 0; 
jreset     = 0; 
% jplot      = 0; % for deformed plotting
flag       = 0; % for time integration schemes
% dynamic = 0; % for dynamic analysis
jsave = 1;
fric=0;
% flag_num_scheme = 0; % flag for time-integration scheme
gpuDevice();

%% Mesh generation

meshNo      = 1;
bulkType  = 'P1Q1';
elemType  = 'N2.1';

meshGeneration; % mesh generation code for body

figure;
nrbctrlplot(nrbObjCom); hold on;
axis equal; axis off; view(2);
set(gcf, 'color', 'w');
set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);
saveas(gcf, 'Original_configuration_cp'); 

figure;
nrbkntplot(nrbObjCom); hold on;
axis equal; axis off; view(2);
set(gcf, 'color', 'w');
set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);
saveas(gcf, 'Original_configuration_knt');

%close all

%% Boundary conditions

% Identification of boundaries
crvLeft   = find(Xn(:,1) == min(Xn(:,1))); % left edge
crvLeft = crvLeft';
crvRightOrg  = find(Xn(:,1) == max(Xn(:,1))); % right edge
cent_nde_dist       = (length(crvRightOrg)+1)/2;        % It holds for od no of control points along y-axis
centerNode          = crvRightOrg(cent_nde_dist);       % Stores control point No at x=Lx, y= Ly/2
crvRight           = setdiff(crvRightOrg,centerNode);  % Stores no of control point at x=Lx, y=0:Ly except y=Ly/2
crvRight           =crvRight';
crvBottom = find(Xn(:,2) == min(Xn(:,2))); % bottom edge

crvTop    = find(Xn(:,2) > max(Xn(:,2)) - 0.002 ); % top edge

% Global numbers of boundaries
crvLeftIndex   = [2 ,      2* crvLeft - 1]; % left edge (except the first node, left end is constrained along x)
crvRightIndex  = [2 * crvRight - 1;     2 * crvRight]; % right edge
centerNodeIndex    = [2*centerNode-1, 2*centerNode];     % Stores dof of control point at x=Lx, y=Ly/2
crvBottomIndex = [2 * crvBottom - 1,    2 * crvBottom]; % bottom edge
crvTopIndex    = [2 * crvTop - 1,       2 * crvTop]; % top edge

% Reduced global number due to given boundary conditions
crvRightIndex  = crvRightIndex(:)';
if fric == 0 && bcType == 2
    % 1. Remove dofs of fixed bcs
    ir = setdiff(1:nDof,crvLeftIndex) ; 
    % 2. Subtracte non-homogenenous BC.
    ir  = setdiff(ir,crvRightIndex) ;
    % 3. size of dofs to be analysed.
    nir = size(ir,2) ;
end

ixo   = centerNodeIndex(1,1:2:end);                    % x-dof of center nodes
izo   = centerNodeIndex(1,2:2:end);                    % y-dof of center nodes
iux   = crvRightIndex(1,1:2:end);                     % x-dof of right edge except for center node
iuz   = crvRightIndex(1,2:2:end);                     % y-dof of right edge except for center node
% Determine size of right boundary's dofs
siux = size(iux,2);
siuz = size(iuz,2);

%% Material parameters (Neo-Hookean model)

E      = 1e0; % Young's modulus
nu     = 0.2; % poisson's ratio
mu     = E/2/(1+nu); % shear modulus (Lame constant)
lambda = 2*mu*nu/(1 - 2*nu); % bulk modulus (Lame constant)
kappa  = lambda + 2*mu/3;
Ip4      = [ones(3,3), zeros(3,1); zeros(1,3), 0]; % (I tensor product I) % [1 1 1 0; 1 1 1 0; 1 1 1 0; 0 0 0 0]
Is4      = [2*eye(3), zeros(3,1); zeros(1,3), 1]; % 2*(fourth order identity) % [2 0 0 0; 0 2 0 0; 0 0 2 0; 0 0 0 1]
D      = lambda*Ip4 + mu*Is4; % one term of elasticity tensor expression

%% Computational parameters

UB    = 120 ;
UBo   = 0 ;                       % Initial displacement/rotation of the strip end
e0    = 0 ;                       %
es    = 90;                     % Number of steps in which the UB has to be acheieved
jpull = 0 ;    
ubdis = (Xn(crvRight,2)-Ly/2);

%--------------------
% Factor after which the after which the application of displcement is
% stoped and the rammed down. This happens when difference in cooradinate
% of node 1 at start and at time t_n+1 is more than this factor times
% length of strip in Y direction i.e. Ly
factor_until_which_disp_applied   = 1 ;   % Used in feangle.m and feangle_composite_schemes.m
factor_until_which_moment_applied = 0.3 ; % Used in feangle_for_moment.m and feangle_for_moment_composite_schemes.m

% Factor after which the after which the code is terminated as the strip
% has completely peeled off from the surface. This happens when difference in cooradinate
% of node 1 at start and at time t_n+1 is more than this factor times
% length of strip in Y direction i.e. Ly
flag_terminate_code_at_complete_peeloff = 1 ; % 0 ---> dont terminate code after strip completely peels off from the ssurface
% 1 ---> terminate
factor_after_which_code_terminates      = 1 ; % Used in IGA_main.m in the time loop as soon as time loop starts.
%--------------------

restart_counter     = 0 ;             % counter which counts the number of times a increment has been cut for convergence
max_restart_counter = 1 ;

%% Initialization of quantities

% Displacement, velocity, and acceleration
% UGOld = zeros(nDof_A + nDof_B, 1);
% VGOld = zeros(nDof_A + nDof_B, 1);
% AGOld = zeros(nDof_A + nDof_B, 1);
Ug    = zeros(nDof, 1);
DU    = zeros(nDof, 1);
% VG    = zeros(nDof_A + nDof_B, 1);
% AG    = zeros(nDof_A + nDof_B, 1);

% External force vectors
F   = zeros(nDof, 1);
ldv = 1;

% Body forces
b = [0; 0];

 % Load value
Reas= 0;

%  Initial Boundary displacement due to prerotation (except center node)
uprx  = zeros(size(crvRight,2),1) ;
uprz  = zeros(size(crvRight,2),1) ;

%  Initial value of rotation converted from degree to radians
phico = pi/180*UBo ;

% Other values
ttime = 0; % time
ttimeOld = 0;  % time from last restart
ls    = 0; % increment number
lsOld    = 0;  % increment number from last restart
ubOld    = 0;  % displacement from last restart
ub        = 0 ;         % Initialize the total angle

%  Initial x- and y-displacement values of right edge's nodes (except the center node)
ubdot = 0 ;
ubx   = 0 ;
ubz   = 0 ;
flag_coming = 0 ;
% lflag   = 0;  % This flag is used while writing the history data to file
% totalSuppliedEnergy = 0;  % total energy supplied

% % Total energy
% if jreset == 0
%     totalSuppliedEnergyVector(1, :) = [0 0 0];
% end

% Details about Gauss quadrature scheme

multi = 1;
jgaus = multi * (orderXiCont + 1); fequad1; xigs = xig1; ngps = jgaus; % GPs for the evaluations of contact surface integral 
jgaus = multi * (orderXiCont + 1); fequad2; xigve = xig2; ngpve = jgaus^2; % GPs for the evaluations of contact element integral 
jgaus = 1 * (orderXi + 1); fequad2; xigv  = xig2; ngpv = jgaus^2; % GPs for the evaluations of bulk element integral 

%---------------------------------------------------------------------------------
% Adhesion Parameters
%---------------------------------------------------------------------------------
% 1. Interaction parameters
yL = 2.5;
yW = 25.3;
yE = yW*yL^3 ;
c1 = pi/45/yE/yL^6 ;
c2 = pi/3/yE ;
% Steinmann discrete gradient method
yo = -(c1/c2)^(1/6) ;       % r_eq

% 3. Regularized adhesive force function
cFmod = 1;
rkPs  = (1/15)^(1/6)/yL ;   % r_eq = (1/15)^(1/6). ro
rm    = cFmod*rkPs ;
Fm    =    c1/rm^9    - c2/rm^3 ;
dFm   = -9*c1/rm^10 + 3*c2/rm^4 ;
%-------
% Parameters for adhesive frictional contact
%-------
actset = zeros(nxc,ngps);
gts    = zeros(nxc,ngps);
gtsold = gts;
xiP    = zeros(nxc,ngps);
xiPo   = xiP;
aslide = zeros(nxc,ngps);

% 4. Specify the rate of loading that is to be achieved
thetaa = 0;                  % No ramming is required for static analysis

% 5. Specify the rate of loading that is to be achieved. (it is in length per
% unit non dimensionalzed time)

ubdota = UB/es;

% 6. Ta is the time till which the ramming is to be applied
Ta = 2*thetaa/ubdota ;

% For restarting
if jreset == 0
    starttime = 0 ;
else
    starttime = ttime ;
end

% 7. Defining endtime of analysis
if jreset == 0
    % starttime = 0 ;
    % Specify the end time of the analysis
    endTime = (UB - thetaa)/ubdota + Ta ;     
end
 
phic = pi/180 ;                      % quantity transforming from degree to radians
uz   = 0 ;

% Parameters for Newton-Raphson
TOL       = 1e-25; % tolerance limit
epsN      = 1e2; % normal penalty parameter
epsT      = 100; % tangential penalty parameter
actsetTOL = 1.0e-6; % active-set tolerance for contact identification
nsm       = 0; % maximum no. of N-R steps, which changes during analysis.
nlimc_max = 40; % If N-R iter. are more than nlimc_max & eN < -18 then fine, otherwise terminate
nlim_max  = 300; % If N-R iteration is more than nlim_max terminate the code
eN_slowConvergence = -20;

% Time steps for transient analysis
dtA         = 1; % fixed time step of 1
dt          = dtA;
dtOriginal = dtA; % save the time period with which we started
dtOld      = dtOriginal;

%---------------------------------------------------------------------------------
% Display Messages in Command Window
%---------------------------------------------------------------------------------

fprintf('\nStarttime = %g, Endtime = %g, Ramming Angle = %g \n',starttime,endTime,thetaa);
fprintf('Time step = %g, Ramming Time = %g, Rate of angle change = %g\n',dt,Ta,ubdota);   

%% Parallel computation
nksz   = nDofCont_Elem * nDofCont_Elem * nElemCont + ...
            nDofBulk_Elem * nDofBulk_Elem * nElemBulk; % total number of entries of all the element stiffness matrices
isp    = zeros(nksz, 1); % column vector of this size
jsp    = zeros(nksz, 1); % column vector of this size
ksp    = zeros(nksz, 1); % column vector of this size
% Contact part
for i = 1:nElemCont % loop over all the contact elements
    con       = 2 * cpConArrayCont(i, :); 
    ie(1, 1:2:nDofCont_Elem) = con - 1; % x-DOF array
    ie(1, 2:2:nDofCont_Elem) = con; % y-DOF array
    isp((i-1) * nDofCont_Elem^2 + 1 : i * nDofCont_Elem^2, 1) ...
          = repmat(ie, [1, nDofCont_Elem])'; % storing all the DOF numbers in a column vector
    jsptemp = repmat(ie, [nDofCont_Elem, 1]);
    jsp((i-1) * nDofCont_Elem^2 + 1 : i * nDofCont_Elem^2, 1) ...
          = jsptemp(:); % storing all the DOF numbers in a column vector
end
clear con ie jsptemp
% Bulk part
for i = 1:nElemBulk % loop over all the bulk elements
    con          = 2 * cpConArrayBulk(i, :);
    ie(1, 1:2:nDofBulk_Elem) = con - 1; % x-DOF array
    ie(1, 2:2:nDofBulk_Elem) = con; % y-DOF array
    isp(nDofCont_Elem * nDofCont_Elem * nElemCont + ...
        (i-1) * nDofBulk_Elem^2 + 1 : nDofCont_Elem * nDofCont_Elem ...
        * nElemCont + i * nDofBulk_Elem^2, 1) ...
        = repmat(ie, [1, nDofBulk_Elem])'; % storing all the DOF numbers in a column vector
    jsptemp = repmat(ie, [nDofBulk_Elem, 1]);
    jsp(nDofCont_Elem * nDofCont_Elem * nElemCont + ...
        (i-1) * nDofBulk_Elem^2 + 1 : nDofCont_Elem * nDofCont_Elem ...
        * nElemCont + i * nDofBulk_Elem^2, 1) = jsptemp(:); % storing all the DOF numbers in a column vector
end
clear con ie jsptemp

%% In case iteration do not converge (for restart)

itrcon_flag = 1;

% coordinate storing array for animation making
Uh             = zeros(endTime, 3); % Stores: [increment_no ttime Moment1] ;
xn_storeG        = zeros(nCP,2, endTime ); % current coordinate vector (endTime+1)

% Other arrays
time_inc_start = zeros(endTime, 1);
if jreset == 0
    R_contact = zeros(nDof, 1);
end

% Contact parameters
nElemCont = nElemXi;
nCPCont   = nCPXiCont; % CCP                = noCnpnue;
xipm_tn   = zeros(nElemCont, ngps);
xiPS_old  = zeros(nElemCont, ngps);
xiPS_Tn   = zeros(nElemCont, ngps);


%% Load step time and element operation time

% Load step computation time
startTimeLoad   = zeros(endTime + 1, 1);
endTimeLoad     = zeros(endTime + 1, 1);

% Newton-Rapshon iteration convergence time
startTimeNR     = zeros(endTime + 1, 1);
endTimeNR       = zeros(endTime + 1, 1);

% Contact part computation time
startTimeContPart   = zeros(endTime + 1, 27);
endTimeContPart     = zeros(endTime + 1, 27);

% Bulk part computation time
startTimeBulkPart = zeros(endTime + 1, 27);
endTimeBulkPart   = zeros(endTime + 1, 27);

% Contact computation time
startTimeCont = zeros(endTime + 1, 27);
endTimeCont   = zeros(endTime + 1, 27);

% Contact computation time (same as above, check)
cont_computation_time_start = zeros(1, endTime + 1);
cont_computation_time_end   = zeros(1, endTime + 1);

% Solver time
startTimeSolver = zeros(endTime + 1, 27);
endTimeSolver   = zeros(endTime + 1, 27);


%% End of pre-processing

endTimePreprocessor = toc;
fprintf('\nPre-processing time: %.6f \n', endTimePreprocessor);
fprintf('---------------------------------------------------------------');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %  
%                               PROCESSING                                %  
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

initGPUArray;
wcounter = 1;
while ttime < endTime
    counter_20   = 0 ;
        if jreset == 0
            flag_no_more_disp = 0;
        end
        UB            = zeros(endTime,4);
    % For Bathe's scheme
    if flag == 0
        time_inc_start = toc;
    end
    
    % Increment identification loop for Bathe's scheme.
    if (ls == 0 || itrcon_flag == 1) && flag == 0
        xnOld  = xn;    
        UgOld    = Ug;      
%         VG_old    = VG;      
%         AG_old    = AG;
        lsOld    = ls;       
        ttimeOld = ttime;
        
    elseif itrcon_flag == 0
        error('\nN-R iterations did not converge at Load-step %g', ls);
    end % end of increment identification loop
    
    % Update load/time increment number
    if flag == 0
        ls = ls + 1; % load step
        startTimeLoad(ls, 1) = toc;
    end
    
    if (itrcon_flag == 1  ||  ls == 1) && flag == 0
            fprintf('\n ------------------------- \n') ;
            fprintf('\n\nStarting load step no: %g\n\n',ls-1) ;
            fprintf('Rate of angle change desired for the load step: %g is: %g \n',ls-1,ubdota) ;
            fprintf('\n\nLoad-step size %g, Load-step %g of %g \n', ubdota * dt, ls, es);
    end
    
    
     %--------------------------------------------------------------------------
     % Angle Change triggered displacement values
     %--------------------------------------------------------------------------
     % Exact Value of displacement
     
     % Angle Change
     % Applied moment
     feangle_for_moment_IGA    % Static analysis
        
     % Newton-Raphson loop
     newtonLoop;
     
     if itrcon_flag == 1 && jterminate == 0        
        if flag == 0
            % load history, time, and other parameters
            Uh(ls, :)                = [ls, ttime, Reas];
            endTimeLoad(ls, 1)         = toc - startTimeLoad(ls, 1);
            xn_storeG(:,1:2, ls)      = xn;
            
            % Increment in time, based on time-integration scheme
%             if flag_num_scheme == 0 || flag_num_scheme == 1 % Static or Newmark
            if bcType == 2
               ttime = ttime + dt;
            elseif bcType == 1
               ttime = ttime + 1;
            end
            
            fprintf('\nTime taken by Load-step %g: %.6f', ls, endTimeLoad(ls, 1));
            fprintf('\n---------------------------------------------------------------');
            ttime_when_error = ttime ;
        end
     else % elseif itrcon_flag == 0 && jterminate == 0
        fprintf('\nN-R iterations did not converge.');    
     end
    wcounter = wcounter + 1;
end
    
close all

%% End of processing
endTimeProcessor = toc - endTimePreprocessor;
fprintf('\n\nProcessing time: %.6f \n', endTimeProcessor);
fprintf('---------------------------------------------------------------');


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                            POST-PROCESSING                              %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Plot displacement and stress contour
zz = ls;
% Buliding projected mesh using four-noded quadrilateral elements (using current coordinates)
build_Q4Mesh;

% Compute stresses at the nodes of projected mesh
computeDispNStress_Q4Mesh;

% Draw displacement and stress at the nodes of projected mesh
jfield = 1;
jsmooth = 1; % node type flag for displacement and stress plot
jdisp = 2; 
jstress = 2;

drawDisp_Q4Mesh;
for j=[1,2,3,4]
    jstress=j;
    drawStress_Q4Mesh;
end
%close all

%% End of post-processing
endTimePostprocessor = toc - endTimeProcessor;
fprintf('\n\nPost-processing time: %.6f \n', endTimePostprocessor);
fprintf('---------------------------------------------------------------');

totalTime = toc;
fprintf('\n\nTotal time: %.6f \n', totalTime);
