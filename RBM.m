% 1. Define points [x, y, z]
P = [0, 0, 0;   % Point 1 (Origin)
     10, 0, 0]; % Point 2 (10 units along X)

% 2. Define edges
E = [1,2];

% 3. Define constraints [PointID, DOF (1-6), Value]
C = [
    1, 1, 0;  % Pt 1: Ux = 0
    1, 2, 0;  % Pt 1: Uy = 0
    1, 3, 0;  % Pt 1: Uz = 0
    1, 4, 0;  % Pt 1: Rx = 0
    1, 5, 0;  % Pt 1: Ry = 0
    1, 6, pi/2; % Pt 1: Rz = 90 degrees
];

% 4. Run Solver and Plot
[~, U_all] = solve_rbm(P, C);

disp('Calculated Displacements [Ux, Uy, Uz, Rx, Ry, Rz]:');
disp(U_all);

plot_rbm(P, U_all, E);