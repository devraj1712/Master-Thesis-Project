% Buliding projected mesh using four-noded quadrilateral elements (using current coordinates)
build_Q4Mesh;

% Compute stresses at the nodes of projected mesh
computeDispNStress_Q4Mesh;

% Draw stress at the nodes of projected mesh
drawStress_Q4Mesh;

% Draw displacement at the nodes of projected mesh
drawDisp_Q4Mesh;