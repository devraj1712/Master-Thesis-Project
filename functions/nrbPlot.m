function nrbPlot(nrbObject)

% Objective of the function: to plot the NURBS geometry according to
% the control points and the knot vectors

% INPUT
% structure :: nrbsObject  = NURBS-based object

% OUTPUT
% figure    :: nrbctrlplot = plotting geometry showing control points
% figure    :: nrbkntplot  = plotting geometry showing knot locations
    
    figure;
    nrbctrlplot(nrbObject); view(2);
    figure;
    nrbkntplot(nrbObject); view(2);
    
end