
PX      = cpCorArray_Q4Mesh; % coordinates array for the NURBS geometry of four-noded quadrilateral elements 
CON     = cpConArray_Q4Mesh; % control point connectivity array for four-noded quadrilateral elements
%4-----3
%|     |
%|     |
%1-----2
% [1 2 3 4]

% Plot in element-wise manner
figure;
for i = 1:nElem_Proj % loop over number of projected elements for the body
    
    if jdisp == 0
        % [1 2 3 4 1]
        X = [PX(CON(i,1), 1) PX(CON(i,2), 1) PX(CON(i,3), 1) PX(CON(i,4), 1) PX(CON(i,1), 1)];
        Y = [PX(CON(i,1), 2) PX(CON(i,2), 2) PX(CON(i,3), 2) PX(CON(i,4), 2) PX(CON(i,1), 2)];
        plot(X, Y, 'k'); hold on;
        
    else  % Plot stresses onto mesh
        if jsmooth == 0 % four-noded
            %4-----3
            %|     |
            %1-----2
            % [1 2 4 3]
            Xp = [PX(CON(i,1),1) PX(CON(i,2), 1) ; PX(CON(i,4), 1) PX(CON(i,3), 1)]; % x-coordinates of nodes of ith projected element
            Yp = [PX(CON(i,1),2) PX(CON(i,2), 2) ; PX(CON(i,4), 2) PX(CON(i,3), 2)]; % y-coordinates of nodes of ith projected element
            Dp = [dispNode(CON(i,1), jdisp) dispNode(CON(i,2), jdisp);
                dispNode(CON(i,4), jdisp) dispNode(CON(i,3), jdisp)]; % (jstress)th component of stress of nodes of ith projected element
            pcolor(Xp, Yp, Dp); hold on;
            
        elseif jsmooth == 1 % three-noded
            %4-----3
            %|  5  |
            %1-----2
            Xn1 = PX(CON(i,1), :); % coordinate of 1st point for ith projected element
            Xn2 = PX(CON(i,2), :); % coordinate of 2nd point for ith projected element
            Xn3 = PX(CON(i,3), :); % coordinate of 3rd point for ith projected element
            Xn4 = PX(CON(i,4), :); % coordinate of 4th point for ith projected element
            Xn5 = (Xn1 + Xn2 + Xn3 + Xn4) / 4; % coordinate of 5th point for ith projected element
            
            Dp1 = dispNode(CON(i,1), jdisp) ; % (jdisp)th component of displacement at 1st point for ith projected element
            Dp2 = dispNode(CON(i,2), jdisp) ; % (jdisp)th component of displacement at 2nd point for ith projected element
            Dp3 = dispNode(CON(i,3), jdisp) ; % (jdisp)th component of displacement at 3rd point for ith projected element
            Dp4 = dispNode(CON(i,4), jdisp) ; % (jdisp)th component of displacement at 4th point for ith projected element
            Dp5 = (Dp1 + Dp2 + Dp3 + Dp4) / 4 ; % (jdisp)th component of displacement at 5th point for ith projected element
            
            % [1 2 5 5], 1st sub-element
            Xp = [Xn1(1) Xn2(1); Xn5(1) Xn5(1)]; % x-coordinates of nodes of 1st sub-element of ith projected element
            Yp = [Xn1(2) Xn2(2); Xn5(2) Xn5(2)]; % y-coordinates of nodes of 1st sub-element of ith projected element
            Dp = [Dp1 Dp2; Dp5 Dp5]; % (jdisp)th component of displacement of nodes of 1st sub-element of ith projected element
            s=pcolor(Xp, Yp, Dp); hold on;
            set(s,'EdgeColor','none');
            
            % [2 3 5 5], 2nd sub-element
            Xp  = [Xn2(1)  Xn3(1) ; Xn5(1)  Xn5(1)] ;
            Yp  = [Xn2(2)  Xn3(2) ; Xn5(2)  Xn5(2)] ;
            Dp  = [Dp2 Dp3; Dp5 Dp5] ;
            s=pcolor(Xp, Yp, Dp); hold on ;
            set(s,'EdgeColor','none');
            
            % [3 4 5 5], 3rd sub-element
            Xp  = [Xn3(1) Xn4(1); Xn5(1) Xn5(1)];
            Yp  = [Xn3(2) Xn4(2); Xn5(2) Xn5(2)];
            Dp  = [Dp3 Dp4; Dp5 Dp5];
            s=pcolor(Xp, Yp, Dp); hold on;
            set(s,'EdgeColor','none');
            
            % [4 1 5 5], 4th sub-element
            Xp  = [Xn4(1) Xn1(1); Xn5(1) Xn5(1)];
            Yp  = [Xn4(2) Xn1(2); Xn5(2) Xn5(2)];
            Dp  = [Dp4 Dp1; Dp5 Dp5];
            s=pcolor(Xp, Yp, Dp); hold on;
            set(s,'EdgeColor','none');
            
        end
    end
end
axis auto; colorbar
 axis equal; axis off; view(2);
 set(gca,'fontsize', 25); set(gcf, 'color', 'w');
 set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);
 %saveas(gcf, 'Output_Disp.png');
 colormap jet
%{
 clear body;
body                = nrbObj;
body.coefs(1:2, :)  = xn_storeG(1:end, :, zz)';
body.coefs(1, :)    = body.coefs(1, :) .* body.coefs(4, :); % [x] >> [wx]
body.coefs(2, :)    = body.coefs(2, :) .* body.coefs(4, :); % [y] >> [wy]
drawMesh(body, nElemXi_Proj, nElemEta_Proj);
%}