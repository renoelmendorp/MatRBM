function [global_dof, U_all] = solve_rbm(P, C)
    % Inputs:
    % P: [M x 3] matrix of point coordinates [x, y, z]
    % C: [N x 3] matrix of constraints [PointID, DOF_Index, Value]
    %    DOF_Index: 1=Ux, 2=Uy, 3=Uz, 4=RotX, 5=RotY, 6=RotZ
    
    num_constraints = size(C, 1);
    
    %% 1. Check Constraint Status (Instantaneous Kinematics)
    % We build the linear Jacobian evaluated at zero displacement
    A = zeros(num_constraints, 6);
    b_lin = zeros(num_constraints, 1);
    
    for k = 1:num_constraints
        pt_idx = C(k, 1);
        dof = C(k, 2);
        b_lin(k) = C(k, 3);
        
        x = P(pt_idx, 1); y = P(pt_idx, 2); z = P(pt_idx, 3);
        
        % Populate the linear kinematic matrix rows based on DOF
        if dof == 1; A(k,:) = [1, 0, 0, 0, z, -y]; end
        if dof == 2; A(k,:) = [0, 1, 0, -z, 0, x]; end
        if dof == 3; A(k,:) = [0, 0, 1, y, -x, 0]; end
        if dof == 4; A(k,:) = [0, 0, 0, 1, 0, 0]; end
        if dof == 5; A(k,:) = [0, 0, 0, 0, 1, 0]; end
        if dof == 6; A(k,:) = [0, 0, 0, 0, 0, 1]; end
    end
    
    rA = rank(A);
    rAug = rank([A, b_lin]);
    
    if rA < 6
        error('Structure is UNDER-CONSTRAINED. Rigid body modes exist. Rank: %d/6.', rA);
    elseif rA == 6 && num_constraints == 6
        fprintf('Structure is EXACTLY CONSTRAINED. Proceeding to solve.\n');
    elseif rA == 6 && num_constraints > 6
        if rA == rAug
            fprintf('Structure is OVER-CONSTRAINED but CONSISTENT. Proceeding.\n');
        else
            fprintf('Structure is OVER-CONSTRAINED and CONTRADICTORY. Solving via Least Squares fit.\n');
        end
    end

    %% 2. Non-Linear Solver
    % Initial guess for the 6 global DOFs: [tx, ty, tz, rx, ry, rz]
    q0 = zeros(6, 1); 
    
    % Optimization options
    options = optimoptions('lsqnonlin', 'Display', 'off', ...
                           'FunctionTolerance', 1e-9, ...
                           'StepTolerance', 1e-9);
    
    % Solve for global parameters driving residuals to zero (or min squared error)
    global_dof = lsqnonlin(@(q) rb_residual(q, P, C), q0, [], [], options);
    
    %% 3. Calculate All Displacements
    num_points = size(P, 1);
    U_all = zeros(num_points, 6); % Store [Ux, Uy, Uz, Rx, Ry, Rz] for each point
    
    % Global rotation matrix
    R = get_rotation_matrix(global_dof(4), global_dof(5), global_dof(6));
    t = global_dof(1:3);
    
    for i = 1:num_points
        X0 = P(i, :)';
        X_new = R * X0 + t;
        U_trans = X_new - X0;
        
        % Record translations and global rotations
        U_all(i, 1:3) = U_trans';
        U_all(i, 4:6) = global_dof(4:6)'; 
    end
    
    fprintf('Solution Found.\n');
end

%% Helper Functions
function res = rb_residual(q, P, C)
    % Calculates the error between prescribed constraints and current state
    t = q(1:3);
    R = get_rotation_matrix(q(4), q(5), q(6));
    
    res = zeros(size(C, 1), 1);
    for k = 1:size(C, 1)
        pt_idx = C(k, 1);
        dof = C(k, 2);
        val = C(k, 3);
        
        if dof <= 3 % Translational constraint
            X0 = P(pt_idx, :)';
            X_new = R * X0 + t;
            U = X_new - X0;
            res(k) = U(dof) - val;
        else % Rotational constraint
            res(k) = q(dof) - val;
        end
    end
end

function R = get_rotation_matrix(rx, ry, rz)
    % Standard XYZ Euler sequence rotation matrix
    Rx = [1, 0, 0; 0, cos(rx), -sin(rx); 0, sin(rx), cos(rx)];
    Ry = [cos(ry), 0, sin(ry); 0, 1, 0; -sin(ry), 0, cos(ry)];
    Rz = [cos(rz), -sin(rz), 0; sin(rz), cos(rz), 0; 0, 0, 1];
    R = Rz * Ry * Rx;
end