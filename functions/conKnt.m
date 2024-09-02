function [kntConArray, kntSpanXi, kntSpanEta] = conKnt(nElemXi, nElemEta, kntVecXi, kntVecEta)

% Objective of the function: to find the knot span and knot vector
% connectivity between the elements as an array

% INPUT          
% integer :: nElemXi       = number of elements along xi   
% integer :: nElemEta      = number of elements along eta 
% real    :: kntVecXi      = knot vector along xi          
% real    :: kntVecEta     = knot vector along eta

% OUTPUT
% integer :: kntConArray   = knot vector connectivity array 
% real    :: kntSpanXi     = knot span along xi             
% real    :: kntSpanEta    = knot span along eta            

    kntVecUnqXi = unique(kntVecXi);   % vector along xi with unique knot points
    kntVecUnqEta = unique(kntVecEta); % vector along eta with unique knot points
    kntSpanXi = [kntVecUnqXi(1:end-1); kntVecUnqXi(2:end)]'; 
    kntSpanEta = [kntVecUnqEta(1:end-1); kntVecUnqEta(2:end)]';
    [array1, array2] = meshgrid(linspace(1, nElemXi, nElemXi), linspace(1, nElemEta, nElemEta));
    array1 = array1';
    array2 = array2';
    kntConArray(:, 1) = array1(:);
    kntConArray(:, 2) = array2(:);

end