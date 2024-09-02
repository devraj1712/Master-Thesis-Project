function drawMesh(nrbObj, nElemXi, nElemEta)

% Draw mesh wire of a given NURBS geometry
hold on;


% Multiply the wight points to the Cps (for drawing exact geometry)
%nurbs.coefs(1,:) = nurbs.coefs(1,:).*nurbs.coefs(4,:);
%nurbs.coefs(2,:) = nurbs.coefs(2,:).*nurbs.coefs(4,:);
%nurbs.coefs(3,:) = nurbs.coefs(3,:).*nurbs.coefs(4,:);


% And plot the knots
nKntXi = nElemXi+1; % number of knots along xi
nKntEta = nElemEta+1; % number of knots along eta
kntVecUnqXi = unique(nrbObj.knots{1}(nrbObj.order(1):end-nrbObj.order(1)+1)); % unique knot vector along xi  (p to n-p+1)
kntVecUnqEta = unique(nrbObj.knots{2}(nrbObj.order(2):end-nrbObj.order(2)+1)); % unique knot vector along eta (from q to m-q+1)

% Coordinates along x for nx_a+1 or knt1 elements
p1   = nrbeval(nrbObj, {linspace(kntVecUnqXi(1), kntVecUnqXi(end), nKntXi),  kntVecUnqEta} );
% Coordinates along y for nx_a+1 or knt1 elements
p2   = nrbeval (nrbObj, {kntVecUnqXi,  linspace(kntVecUnqEta(1), kntVecUnqEta(end), nKntEta) });

if (any(nrbObj.coefs(3,:))) % surface in a 3D space
    for i = 1:numel(kntVecUnqXi)
        plot3(squeeze(p1(1,i,:)), squeeze(p1(2,i,:)), squeeze(p1(3,i,:)), 'k');
    end
    for i = 1:numel(kntVecUnqEta)
        plot3(squeeze(p2(1,:,i)), squeeze(p2(2,:,i)), squeeze(p2(3,:,i)), 'k');
    end
else % plain surface
    for i = 1:numel(kntVecUnqXi)
        plot(squeeze(p2(1,i,:)), squeeze (p2(2,i,:)), '-k');
    end
    for i = 1:numel(kntVecUnqEta)
        plot(p1(1,:,i), p1(2,:,i),'-k');
    end
end

end