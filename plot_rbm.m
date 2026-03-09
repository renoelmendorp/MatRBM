function plot_rbm(P, U_all, Edges, scale, rigid_mode)
    % Inputs:
    % P: [M x 3] original coordinates
    % U_all: [M x 6] results from the solver
    % Edges: [E x 2] connectivity
    % scale: Magnification factor (e.g., 10, 100)
    % rigid_mode: Boolean. If true, maintains shape during large rotations.

    if nargin < 3; Edges = []; end
    if nargin < 4; scale = 1; end
    if nargin < 5; rigid_mode = false; end

    if ~rigid_mode
        % Standard Linear Scaling (Standard in FEA/Structural tools)
        % Warning: Can distort shape if rotations are very large
        P_displaced = P + (U_all(:, 1:3) * scale);
    else
        % Rigid Scaling (Maintains geometry perfectly)
        % We take the global rotation/translation and scale THEM.
        % We assume the rotation was the same for all points (Rigid Body)
        t_scaled = U_all(1, 1:3) * scale; % Scaled translation
        r_scaled = U_all(1, 4:6) * scale; % Scaled rotation angles
        
        % Re-calculate rotation matrix with scaled angles
        Rx = [1, 0, 0; 0, cos(r_scaled(1)), -sin(r_scaled(1)); 0, sin(r_scaled(1)), cos(r_scaled(1))];
        Ry = [cos(r_scaled(2)), 0, sin(r_scaled(2)); 0, 1, 0; -sin(r_scaled(2)), 0, cos(r_scaled(2))];
        Rz = [cos(r_scaled(3)), -sin(r_scaled(3)), 0; sin(r_scaled(3)), cos(r_scaled(3)), 0; 0, 0, 1];
        R_scaled = Rz * Ry * Rx;
        
        P_displaced = (R_scaled * P')' + t_scaled;
    end
    
    % --- Plotting Logic ---
    figure('Color', 'w'); hold on; grid on; axis equal; view(3);
    
    % Plot Original (Ghosted/Gray)
    if ~isempty(Edges)
        patch('Vertices', P, 'Faces', Edges, 'EdgeColor', [0.8 0.8 0.8], ...
              'LineWidth', 1, 'FaceColor', 'none', 'DisplayName', 'Original');
    else
        scatter3(P(:,1), P(:,2), P(:,3), 'k', 'filled');
    end
    
    % Plot Displaced (Magnified)
    if ~isempty(Edges)
        patch('Vertices', P_displaced, 'Faces', Edges, 'EdgeColor', 'r', ...
              'LineWidth', 2, 'FaceColor', 'none', ...
              'DisplayName', ['Displaced (x', num2str(scale), ')']);
    else
        scatter3(P_displaced(:,1), P_displaced(:,2), P_displaced(:,3), 'r', 'filled');
    end
    
    title(['Rigid Body Motion (Magnification: ', num2str(scale), 'x)']);
    legend('Location', 'northeastoutside');
end