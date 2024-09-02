% In this file, current angle is evaluated for the case of static analysis.

    flag_ramming              = 0 ;
    flag_no_more_disp         = 1 ;
    ub_before_load_decreases1 = 0 ;
    ub_before_load_decreases2 = 0 ;
    Tb1                       = 0 ;
    Tb2                       = 0 ;
    ubdot                     = 0 ;
    ubddot                    = 0 ;
    ubdotOld                  = 0 ;
    if jreset == 0
        ub = ttime*ubdota;               % Current angle
    else
        ub = (ls - ls_old ) * dt * ubdota; % Since ub has to be a non zero value
    end
    ub_data(ls,:)= [ttime ub ubdota dt];
    fprintf('---> (Static Analysis) Current Time: %g \t Angle: %g \t Rate of change of angle: %g\n\n',ttime,ub,ubdota);
    % fprintf(fipf,'\n---> (Static Analysis)  Current Time: %g \t Angle: %g \t Rate of change of angle: %g\n\n',ttime,ub,ubdota);


% Before applying the new total displacement, save the old one. This is
% required in fenewtond.m for computing the energy supplied.
% if l == 1
%     ubzOld = 0 ;                    % For the first increment the previous value is 0
%     ubxOld = 0 ;
% else
    ubzOld = ub ;                   % ub is = ttime*ubdota(==ttime*0.02)
    ubxOld = 0  ;
% end


ubz =  ( cos(ub*phic) - cos(phico) ) * ubdis + ub*uz*ones(siuz,1) ;   
ubx = -( sin(ub*phic) - sin(phico) ) * ubdis ;
% Here, ub is initial value of angle which is zero, phic=pi/180, phico = phico = pi/180*UBo,
% ubdis = [ ubdis ; Xn(i,2)-Ly/2 ], ub & uz = 0 is defined initially.