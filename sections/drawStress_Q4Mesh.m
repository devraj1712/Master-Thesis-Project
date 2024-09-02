
jsmooth = 1; % node type flag for stress plot

PX      = cpCorArray_Q4Mesh; % coordinates array for the NURBS geometry of four-noded quadrilateral elements 
CON     = cpConArray_Q4Mesh; % control point connectivity array for four-noded quadrilateral elements
%4-----3
%|     |
%|     |
%1-----2
% [1 2 3 4]

% sigma(I1) = sigma_11 + sigma_22 + sigma_33;
stress_I1 = stressNode(:, 1) + stressNode(:, 2) + stressNode(:, 3);
    figure;
    %jstress = 1; % stress plot flag
    hold on; % Plot in an element-wise manner
    for i = 1:nElem_Proj % loop over number of projected elements for the body

        if jstress == 0 % Plot mesh wired/colored
            % [1 2 3 4 1]
            X = [PX(CON(i,1), 1) PX(CON(i,2), 1) PX(CON(i,3), 1) PX(CON(i,4), 1) PX(CON(i,1), 1)] ;
            Y = [PX(CON(i,1), 2) PX(CON(i,2), 2) PX(CON(i,3), 2) PX(CON(i,4), 2) PX(CON(i,1), 2)] ;        
            plot(X, Y, 'k'); hold on;

        else % Plot stress field onto mesh
            if jsmooth == 0 % four-noded, incorrect interpolation
                %4-----3
                %|     |
                %1-----2
                % [1 2 4 3]
                Xp = [PX(CON(i,1), 1) PX(CON(i,2), 1); PX(CON(i,4), 1) PX(CON(i,3), 1)]; % x-coordinates of nodes of ith projected element
                Yp = [PX(CON(i,1), 2) PX(CON(i,2), 2); PX(CON(i,4), 2) PX(CON(i,3), 2)]; % y-coordinates of nodes of ith projected element
                Sp = [stressNode(CON(i,1), jstress) stressNode(CON(i,2), jstress); 
                        stressNode(CON(i,4), jstress) stressNode(CON(i,3), jstress) ]; % (jstress)th component of stress of nodes of ith projected element
                pcolor(Xp, Yp, Sp); hold on;

            elseif jsmooth == 1 % three-noded, smoother Interpolation
                %4-----3
                %|  5  |
                %1-----2
                Xn1 = PX(CON(i,1), :); % coordinate of 1st point for ith projected element
                Xn2 = PX(CON(i,2), :); % coordinate of 2nd point for ith projected element
                Xn3 = PX(CON(i,3), :); % coordinate of 3rd point for ith projected element
                Xn4 = PX(CON(i,4), :); % coordinate of 4th point for ith projected element
                Xn5 = (Xn1 + Xn2 + Xn3 + Xn4) / 4; % coordinate of 5th point for ith projected element
                
                    Sp1 = stressNode(CON(i,1), jstress); % (jstress)th component of stress at 1st point for ith projected element
                    Sp2 = stressNode(CON(i,2), jstress); % (jstress)th component of stress at 2nd point for ith projected element
                    Sp3 = stressNode(CON(i,3), jstress); % (jstress)th component of stress at 3rd point for ith projected element
                    Sp4 = stressNode(CON(i,4), jstress); % (jstress)th component of stress at 4th point for ith projected element
                    Sp5 = (Sp1 + Sp2 + Sp3 + Sp4) / 4; % (jstress)th component of stress at 5th point for ith projected element

                    % [1 2 5 5], 1st sub-element
                    Xp  = [Xn1(1) Xn2(1); Xn5(1) Xn5(1)]; % x-coordinates of nodes of 1st sub-element of ith projected element
                    Yp  = [Xn1(2) Xn2(2); Xn5(2) Xn5(2)]; % y-coordinates of nodes of 1st sub-element of ith projected element
                    Sp  = [Sp1 Sp2; Sp5 Sp5]; % (jstress)th component of stress of nodes of 1st sub-element of ith projected element
                    pcolor(Xp, Yp, Sp); hold on;

                    % [2 3 5 5], 2nd sub-element
                    Xp  = [Xn2(1)  Xn3(1) ; Xn5(1)  Xn5(1)] ;
                    Yp  = [Xn2(2)  Xn3(2) ; Xn5(2)  Xn5(2)] ;
                    Sp  = [Sp2 Sp3; Sp5 Sp5] ;
                    pcolor(Xp, Yp, Sp); hold on ;

                    % [3 4 5 5], 3rd sub-element
                    Xp  = [Xn3(1) Xn4(1); Xn5(1) Xn5(1)];
                    Yp  = [Xn3(2) Xn4(2); Xn5(2) Xn5(2)];
                    Sp  = [Sp3 Sp4; Sp5 Sp5];
                    pcolor(Xp, Yp, Sp); hold on;

                    % [4 1 5 5], 4th sub-element
                    Xp  = [Xn4(1) Xn1(1); Xn5(1) Xn5(1)];
                    Yp  = [Xn4(2) Xn1(2); Xn5(2) Xn5(2)];
                    Sp  = [Sp4 Sp1; Sp5 Sp5];
                    pcolor(Xp, Yp, Sp); hold on;
               
            else  % First invariant of stress            
                Xp = [PX(CON(i,1), 1) PX(CON(i,2), 1); PX(CON(i,4), 1) PX(CON(i,3), 1)];
                Yp = [PX(CON(i,1), 2) PX(CON(i,2), 2); PX(CON(i,4), 2) PX(CON(i,3), 2)];
                Sp = [stress_I1(CON(i,1), 1) stress_I1(CON(i,2), 1);
                    stress_I1(CON(i,4), 1)  stress_I1(CON(i,3), 1)];
                pcolor(Xp, Yp, Sp); hold on;

            end
    end
end
axis equal; axis off; view(2);
set(gca,'fontsize', 25); set(gcf, 'color', 'w');
set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);
%saveas(gcf, 'Output_Stress.png');
shading interp; % set color shading properties, interp = interpolated setting
colormap jet; % view and set current colormap, jet = colormap name
colorbar; % colorbar showing color scale
%{
clear body;
body                = nrbObj;
body.coefs(1:2, :)  = xn_storeG(1:end, :, zz)';
body.coefs(1, :)    = body.coefs(1, :) .* body.coefs(4, :); % [x] >> [wx]
body.coefs(2, :)    = body.coefs(2, :) .* body.coefs(4, :); % [y] >> [wy]
drawMesh(body, nElemXi_Proj, nElemEta_Proj);
%}