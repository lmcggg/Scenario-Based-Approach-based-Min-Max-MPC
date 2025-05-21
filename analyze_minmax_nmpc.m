%% Analysis of Min-Max NMPC for Nonlinear Tracking Control
% This script compares the standard NMPC and min-max NMPC approaches
% for tracking control of a nonlinear double integrator system with disturbances

clear all; close all; clc;

% Run min-max NMPC simulation first
run_minmax_nmpc = true;     % Set to true to run min-max NMPC
run_standard_nmpc = true;   % Set to true to run standard NMPC
save_results = true;        % Set to true to save results

%% Common parameters for both approaches
Ts = 0.05;              % Sampling time (s)
T_sim_total = 20;       % Total simulation time (s)
N = 5;                 % NMPC prediction horizon
N_d = 5;                % Number of disturbance scenarios for min-max NMPC

% System limits
d_min = -1.5;           % Lower bound on continuous disturbance
d_max = 1.5;            % Upper bound on continuous disturbance
w_min = Ts * d_min;     % Lower bound on discrete disturbance
w_max = Ts * d_max;     % Upper bound on discrete disturbance
u_min = -5;             % Lower bound on control input
u_max = 5;              % Upper bound on control input
x1_min = -5;            % Lower bound on position (optional)
x1_max = 5;             % Upper bound on position (optional)
x2_min = -3;            % Lower bound on velocity (optional)
x2_max = 3;             % Upper bound on velocity (optional)

% Initial state
x_initial = [0; 0];     % Initial state [position; velocity]

% Cost function weights
Q = diag([10, 1]);      % State tracking error weight
R = 0.1;                % Control input weight
P = diag([20, 2]);      % Terminal state weight

% Reference trajectory parameters
A = 1.0;                % Amplitude of the sine reference
omega_ref = 0.5;        % Angular frequency of the reference (rad/s)

% Disturbance parameters for actual system
dist_freq = 1.0;        % Frequency for sinusoidal disturbance

% Function handle for nonlinear dynamics
f_nonlinear = @(x1, x2) -0.5 * x2 - 0.2 * x2^3 - 0.1 * sin(x1);

%% Generate the same disturbance sequence for fair comparison
num_steps = floor(T_sim_total / Ts);
time_vec = (0:num_steps) * Ts;

% Generate reference trajectory
x_ref_history = zeros(2, num_steps+1);
for k = 0:num_steps
    t = k * Ts;
    x_ref_history(1, k+1) = A * sin(omega_ref * t);
    x_ref_history(2, k+1) = A * omega_ref * cos(omega_ref * t);
end

% Generate actual disturbance sequence
d_actual_sequence = zeros(1, num_steps);
for k = 0:num_steps-1
    t = k * Ts;
    d_actual_sequence(k+1) = (d_max - d_min)/2 * sin(dist_freq * t) + (d_max + d_min)/2;
end

%% Run Min-Max NMPC
if run_minmax_nmpc
    disp('Running Min-Max NMPC Simulation...');
    
    % Setup optimization options
    options_minmax = optimoptions('fmincon', ...
        'Algorithm', 'sqp', ...
        'Display', 'iter', ...
        'MaxIterations', 100, ...
        'MaxFunctionEvaluations', 3000, ...
        'ConstraintTolerance', 1e-3, ...
        'OptimalityTolerance', 1e-3);
    
    % Storage for results
    x_minmax = zeros(2, num_steps+1);
    u_minmax = zeros(1, num_steps);
    cpu_time_minmax = zeros(1, num_steps);
    
    % Initialize state
    x_current = x_initial;
    x_minmax(:, 1) = x_current;
    
    % Initial guess for control sequence
    u_prev = zeros(N, 1);
    
    % Main simulation loop
    for k = 0:num_steps-1
        fprintf('Min-Max NMPC: Step %d / %d\n', k+1, num_steps);
        current_time = k * Ts;
        
        % Get reference for prediction horizon
        if k+N+1 <= num_steps+1
            x_ref_horizon = x_ref_history(:, k+1:k+N+1);
        else
            % For steps near the end of simulation, extend the reference by repeating the last value
            last_valid_idx = num_steps+1;
            valid_horizon = x_ref_history(:, k+1:last_valid_idx);
            x_ref_horizon = [valid_horizon, repmat(x_ref_history(:,last_valid_idx), 1, k+N+1-last_valid_idx)];
        end
        u_ref_horizon = zeros(1, N); % Reference control is zero
        
        % Setup NMPC optimization
        tic;
        
        % Decision variables: [u_1, u_2, ..., u_N, z]
        z_guess = 1000; % A high initial value for z
        decision_vars_guess = [u_prev(2:end); 0; z_guess];
        
        % Lower and upper bounds for decision variables
        lb = [u_min * ones(N, 1); -inf];
        ub = [u_max * ones(N, 1); inf];
        
        % Generate disturbance scenarios
        W_scenarios = generate_disturbance_scenarios(N, N_d, w_min, w_max);
        
        % NMPC optimization
        objective = @(decision_vars) decision_vars(end); % Minimize z
        
        % Define constraints function
        constraint_fn = @(decision_vars) min_max_constraints(decision_vars, x_current, ...
            x_ref_horizon, u_ref_horizon, W_scenarios, N, N_d, Ts, ...
            Q, R, P, f_nonlinear, x1_min, x1_max, x2_min, x2_max);
        
        % Solve optimization problem
        [optimal_decision_vars, ~, exitflag, ~] = fmincon(objective, ...
            decision_vars_guess, [], [], [], [], lb, ub, constraint_fn, options_minmax);
        
        cpu_time = toc;
        cpu_time_minmax(k+1) = cpu_time;
        
        % Extract optimal control sequence and apply first control
        U_optimal = optimal_decision_vars(1:N);
        u_prev = U_optimal; % Store for next iteration
        u_applied = U_optimal(1);
        u_minmax(k+1) = u_applied;
        
        % Get actual disturbance
        d_actual = d_actual_sequence(k+1);
        
        % Simulate plant with actual disturbance
        x_next = discrete_model_actual_plant(x_current, u_applied, d_actual, Ts, f_nonlinear);
        
        % Update current state
        x_current = x_next;
        x_minmax(:, k+2) = x_current;
    end
    
    % Save results if needed
    if save_results
        save('minmax_nmpc_results.mat', 'x_minmax', 'u_minmax', 'cpu_time_minmax');
    end
end

%% Run Standard NMPC (without considering worst-case disturbance)
if run_standard_nmpc
    disp('Running Standard NMPC Simulation...');
    
    % Setup optimization options
    options_standard = optimoptions('fmincon', ...
        'Algorithm', 'sqp', ...
        'Display', 'iter', ...
        'MaxIterations', 100, ...
        'MaxFunctionEvaluations', 3000, ...
        'ConstraintTolerance', 1e-3, ...
        'OptimalityTolerance', 1e-3);
    
    % Storage for results
    x_standard = zeros(2, num_steps+1);
    u_standard = zeros(1, num_steps);
    cpu_time_standard = zeros(1, num_steps);
    
    % Initialize state
    x_current = x_initial;
    x_standard(:, 1) = x_current;
    
    % Initial guess for control sequence
    u_prev = zeros(N, 1);
    
    % Main simulation loop
    for k = 0:num_steps-1
        fprintf('Standard NMPC: Step %d / %d\n', k+1, num_steps);
        current_time = k * Ts;
        
        % Get reference for prediction horizon
        if k+N+1 <= num_steps+1
            x_ref_horizon = x_ref_history(:, k+1:k+N+1);
        else
            % For steps near the end of simulation, extend the reference by repeating the last value
            last_valid_idx = num_steps+1;
            valid_horizon = x_ref_history(:, k+1:last_valid_idx);
            x_ref_horizon = [valid_horizon, repmat(x_ref_history(:,last_valid_idx), 1, k+N+1-last_valid_idx)];
        end
        u_ref_horizon = zeros(1, N); % Reference control is zero
        
        % Setup NMPC optimization
        tic;
        
        % Decision variables: U = [u_1, u_2, ..., u_N]
        decision_vars_guess = [u_prev(2:end); 0];
        
        % Lower and upper bounds for decision variables
        lb = u_min * ones(N, 1);
        ub = u_max * ones(N, 1);
        
        % NMPC optimization (assuming zero disturbance in prediction)
        objective = @(U) standard_nmpc_cost(U, x_current, x_ref_horizon, u_ref_horizon, ...
            N, Ts, Q, R, P, f_nonlinear);
        
        % Define constraints function
        constraint_fn = @(U) standard_nmpc_constraints(U, x_current, N, Ts, ...
            f_nonlinear, x1_min, x1_max, x2_min, x2_max);
        
        % Solve optimization problem
        [U_optimal, ~, exitflag, ~] = fmincon(objective, ...
            decision_vars_guess, [], [], [], [], lb, ub, constraint_fn, options_standard);
        
        cpu_time = toc;
        cpu_time_standard(k+1) = cpu_time;
        
        % Store for next iteration
        u_prev = U_optimal;
        u_applied = U_optimal(1);
        u_standard(k+1) = u_applied;
        
        % Get actual disturbance (same as for min-max NMPC for fair comparison)
        d_actual = d_actual_sequence(k+1);
        
        % Simulate plant with actual disturbance
        x_next = discrete_model_actual_plant(x_current, u_applied, d_actual, Ts, f_nonlinear);
        
        % Update current state
        x_current = x_next;
        x_standard(:, k+2) = x_current;
    end
    
    % Save results if needed
    if save_results
        save('standard_nmpc_results.mat', 'x_standard', 'u_standard', 'cpu_time_standard');
    end
end

%% Compare Results
if run_minmax_nmpc && run_standard_nmpc
    % Load results if not running the simulations
    if ~run_minmax_nmpc && exist('minmax_nmpc_results.mat', 'file')
        load('minmax_nmpc_results.mat');
    end
    if ~run_standard_nmpc && exist('standard_nmpc_results.mat', 'file')
        load('standard_nmpc_results.mat');
    end
    
    % Compare tracking performance
    figure(1);
    subplot(2,1,1);
    plot(time_vec, x_minmax(1,:), 'b-', 'LineWidth', 2);
    hold on;
    plot(time_vec, x_standard(1,:), 'g-', 'LineWidth', 2);
    plot(time_vec, x_ref_history(1,:), 'r--', 'LineWidth', 2);
    grid on;
    xlabel('Time (s)');
    ylabel('Position (x_1)');
    legend('Min-Max NMPC', 'Standard NMPC', 'Reference');
    title('Position Tracking Comparison');
    
    subplot(2,1,2);
    plot(time_vec, x_minmax(2,:), 'b-', 'LineWidth', 2);
    hold on;
    plot(time_vec, x_standard(2,:), 'g-', 'LineWidth', 2);
    plot(time_vec, x_ref_history(2,:), 'r--', 'LineWidth', 2);
    grid on;
    xlabel('Time (s)');
    ylabel('Velocity (x_2)');
    legend('Min-Max NMPC', 'Standard NMPC', 'Reference');
    title('Velocity Tracking Comparison');
    
    % Compare control inputs
    figure(2);
    stairs(time_vec(1:end-1), u_minmax, 'b-', 'LineWidth', 2);
    hold on;
    stairs(time_vec(1:end-1), u_standard, 'g-', 'LineWidth', 2);
    grid on;
    xlabel('Time (s)');
    ylabel('Control Input (u)');
    legend('Min-Max NMPC', 'Standard NMPC');
    title('Control Input Comparison');
    
    % Compare tracking errors
    figure(3);
    tracking_error_minmax = x_minmax(1,:) - x_ref_history(1,:);
    tracking_error_standard = x_standard(1,:) - x_ref_history(1,:);
    
    subplot(2,1,1);
    plot(time_vec, tracking_error_minmax, 'b-', 'LineWidth', 2);
    hold on;
    plot(time_vec, tracking_error_standard, 'g-', 'LineWidth', 2);
    grid on;
    xlabel('Time (s)');
    ylabel('Position Error');
    legend('Min-Max NMPC', 'Standard NMPC');
    title('Position Tracking Error Comparison');
    
    subplot(2,1,2);
    % Calculate cumulative absolute error
    cumulative_error_minmax = cumsum(abs(tracking_error_minmax));
    cumulative_error_standard = cumsum(abs(tracking_error_standard));
    
    plot(time_vec, cumulative_error_minmax, 'b-', 'LineWidth', 2);
    hold on;
    plot(time_vec, cumulative_error_standard, 'g-', 'LineWidth', 2);
    grid on;
    xlabel('Time (s)');
    ylabel('Cumulative Absolute Error');
    legend('Min-Max NMPC', 'Standard NMPC');
    title('Cumulative Absolute Error Comparison');
    
    % Compare computation time
    figure(4);
    stairs(time_vec(1:end-1), cpu_time_minmax, 'b-', 'LineWidth', 2);
    hold on;
    stairs(time_vec(1:end-1), cpu_time_standard, 'g-', 'LineWidth', 2);
    grid on;
    xlabel('Time (s)');
    ylabel('Computation Time (s)');
    legend('Min-Max NMPC', 'Standard NMPC');
    title('Computation Time Comparison');
    
    % Phase portrait comparison
    figure(5);
    plot(x_minmax(1,:), x_minmax(2,:), 'b-', 'LineWidth', 2);
    hold on;
    plot(x_standard(1,:), x_standard(2,:), 'g-', 'LineWidth', 2);
    plot(x_ref_history(1,:), x_ref_history(2,:), 'r--', 'LineWidth', 2);
    grid on;
    xlabel('Position (x_1)');
    ylabel('Velocity (x_2)');
    legend('Min-Max NMPC', 'Standard NMPC', 'Reference');
    title('Phase Portrait Comparison');
    
    % Print statistical comparison
    disp('Performance Metrics:');
    fprintf('Min-Max NMPC: Mean Absolute Error = %.4f, Max Absolute Error = %.4f\n', ...
        mean(abs(tracking_error_minmax)), max(abs(tracking_error_minmax)));
    fprintf('Standard NMPC: Mean Absolute Error = %.4f, Max Absolute Error = %.4f\n', ...
        mean(abs(tracking_error_standard)), max(abs(tracking_error_standard)));
    fprintf('Min-Max NMPC: Mean Computation Time = %.4f s\n', mean(cpu_time_minmax));
    fprintf('Standard NMPC: Mean Computation Time = %.4f s\n', mean(cpu_time_standard));
end

%% Helper Functions
function W_scenarios = generate_disturbance_scenarios(N, N_d, w_min, w_max)
    % Generates N_d disturbance scenario sequences of length N
    W_scenarios = zeros(N, N_d);
    
    if N_d >= 2
        % First scenario: all w_min (worst-case negative disturbance)
        W_scenarios(:, 1) = w_min * ones(N, 1);
        
        % Second scenario: all w_max (worst-case positive disturbance)
        W_scenarios(:, 2) = w_max * ones(N, 1);
        
        % If N_d >= 3, add a zero disturbance scenario
        if N_d >= 3
            W_scenarios(:, 3) = zeros(N, 1);
        end
        
        % Remaining scenarios: random combinations
        for s = 4:N_d
            W_scenarios(:, s) = w_min + (w_max - w_min) * rand(N, 1);
        end
    else
        % If only one scenario, use random disturbances
        W_scenarios(:, 1) = w_min + (w_max - w_min) * rand(N, 1);
    end
end

function x_next = discrete_model_actual_plant(x_current, u, d_actual, Ts, f_nonlinear)
    % Simulates the actual plant for one time step with disturbance
    x1 = x_current(1);
    x2 = x_current(2);
    
    % Euler discretization of the nonlinear dynamics
    x1_next = x1 + Ts * x2;
    x2_next = x2 + Ts * (f_nonlinear(x1, x2) + u + d_actual);
    
    x_next = [x1_next; x2_next];
end

function [c, ceq] = min_max_constraints(decision_vars, x_current, x_ref_horizon, ...
    u_ref_horizon, W_scenarios, N, N_d, Ts, Q, R, P, f_nonlinear, ...
    x1_min, x1_max, x2_min, x2_max)
    % Constraint function for min-max NMPC
    % Extract control sequence and auxiliary variable z
    U_sequence = decision_vars(1:N);
    z = decision_vars(end);
    
    % Initialize constraints
    c = [];    % Inequality constraints: c <= 0
    ceq = [];  % Equality constraints: ceq = 0
    
    % State constraints for all scenarios
    x1_constraints_lower = [];
    x1_constraints_upper = [];
    x2_constraints_lower = [];
    x2_constraints_upper = [];
    
    % Cost constraint for each scenario: J_s - z <= 0
    cost_constraints = zeros(N_d, 1);
    
    % For each disturbance scenario
    for s = 1:N_d
        W_s = W_scenarios(:, s);  % Disturbance sequence for scenario s
        
        % Initialize state prediction for scenario s
        x_pred_s = x_current;
        cost_s = 0;
        
        % Prediction along the horizon
        for i = 1:N
            u_i = U_sequence(i);
            w_i_s = W_s(i);
            
            % Current state
            x1_i = x_pred_s(1);
            x2_i = x_pred_s(2);
            
            % Reference state and control for this step
            x_ref_i = x_ref_horizon(:, i);
            u_ref_i = u_ref_horizon(i);
            
            % Compute stage cost
            state_error = x_pred_s - x_ref_i;
            control_error = u_i - u_ref_i;
            stage_cost = state_error' * Q * state_error + R * control_error^2;
            cost_s = cost_s + stage_cost;
            
            % Predict next state for scenario s
            x1_next = x1_i + Ts * x2_i;
            x2_next = x2_i + Ts * (f_nonlinear(x1_i, x2_i) + u_i) + w_i_s;
            x_pred_s = [x1_next; x2_next];
            
            % Add state constraints for this prediction step
            x1_constraints_lower = [x1_constraints_lower; x1_next - x1_min];
            x1_constraints_upper = [x1_constraints_upper; x1_max - x1_next];
            x2_constraints_lower = [x2_constraints_lower; x2_next - x2_min];
            x2_constraints_upper = [x2_constraints_upper; x2_max - x2_next];
        end
        
        % Add terminal cost
        x_ref_terminal = x_ref_horizon(:, N+1);
        terminal_error = x_pred_s - x_ref_terminal;
        terminal_cost = terminal_error' * P * terminal_error;
        cost_s = cost_s + terminal_cost;
        
        % Cost constraint for scenario s: J_s - z <= 0
        cost_constraints(s) = cost_s - z;
    end
    
    % Combine all inequality constraints
    c = [cost_constraints; 
         -x1_constraints_lower; 
         -x1_constraints_upper; 
         -x2_constraints_lower; 
         -x2_constraints_upper];
    
    % No equality constraints
    ceq = [];
end

function J = standard_nmpc_cost(U, x_current, x_ref_horizon, u_ref_horizon, ...
    N, Ts, Q, R, P, f_nonlinear)
    % Cost function for standard NMPC (assumes zero disturbance in prediction)
    
    % Initialize
    x_pred = x_current;
    J = 0;
    
    % Loop over prediction horizon
    for i = 1:N
        % Current control
        u_i = U(i);
        
        % Reference values
        x_ref_i = x_ref_horizon(:, i);
        u_ref_i = u_ref_horizon(i);
        
        % Stage cost
        state_error = x_pred - x_ref_i;
        control_error = u_i - u_ref_i;
        stage_cost = state_error' * Q * state_error + R * control_error^2;
        J = J + stage_cost;
        
        % Predict next state (assuming zero disturbance)
        x1_next = x_pred(1) + Ts * x_pred(2);
        x2_next = x_pred(2) + Ts * (f_nonlinear(x_pred(1), x_pred(2)) + u_i);
        x_pred = [x1_next; x2_next];
    end
    
    % Terminal cost
    terminal_error = x_pred - x_ref_horizon(:, N+1);
    terminal_cost = terminal_error' * P * terminal_error;
    J = J + terminal_cost;
end

function [c, ceq] = standard_nmpc_constraints(U, x_current, N, Ts, ...
    f_nonlinear, x1_min, x1_max, x2_min, x2_max)
    % Constraints for standard NMPC
    
    % Initialize constraints
    x1_constraints_lower = [];
    x1_constraints_upper = [];
    x2_constraints_lower = [];
    x2_constraints_upper = [];
    
    % Initialize state prediction
    x_pred = x_current;
    
    % Loop over prediction horizon
    for i = 1:N
        % Current control
        u_i = U(i);
        
        % Predict next state (assuming zero disturbance)
        x1_next = x_pred(1) + Ts * x_pred(2);
        x2_next = x_pred(2) + Ts * (f_nonlinear(x_pred(1), x_pred(2)) + u_i);
        x_pred = [x1_next; x2_next];
        
        % Add state constraints
        x1_constraints_lower = [x1_constraints_lower; x1_next - x1_min];
        x1_constraints_upper = [x1_constraints_upper; x1_max - x1_next];
        x2_constraints_lower = [x2_constraints_lower; x2_next - x2_min];
        x2_constraints_upper = [x2_constraints_upper; x2_max - x2_next];
    end
    
    % Combine all inequality constraints
    c = [-x1_constraints_lower; 
         -x1_constraints_upper; 
         -x2_constraints_lower; 
         -x2_constraints_upper];
    
    % No equality constraints
    ceq = [];
end 