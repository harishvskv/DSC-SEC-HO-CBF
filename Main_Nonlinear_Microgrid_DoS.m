%% Main_Nonlinear_Microgrid_DoS.m
% Simulates the Cyber-Physical Microgrid under DoS Attack
% Demonstrates the transition from vulnerable WAN-PI control to DSC-sEC.

clear; clc;

% 1. Load Nonlinear System Parameters
params.M = 8.0; params.D = 1.2; params.gamma = 1.5; % Nonlinear load exponent
params.Tt = 0.3; params.Tg = 0.2; params.R = 0.04;
params.TB = 0.1; params.Tac = 0.2; 
params.C_bat = 500; params.R_th = 2; params.C_th = 8; params.COP = 2.5; 

% 2. Simulation Setup
Tsim = 200; Ts = 0.01; t = 0:Ts:Tsim; N = length(t);
X_baseline = zeros(8, N); % [dw, Pm, Pv, PB, SoC, Tin, PAC, int_dw]
X_proposed = zeros(8, N);

% 3. DoS Attack Definition (Breaching SC-ATD bound)
t_attack_start = 115.0; 
t_attack_end = 118.0; 

% 4. Main Integration Loop
for k = 1:N-1
    time = t(k);
    
    % External Exogenous Disturbances
    d_k = [0.1 * sin(0.5*time); 1.0]; % [P_res, P_load0]
    
    % Check WAN Communication Status
    if time >= t_attack_start && time <= t_attack_end
        is_WAN_compromised = true;
    else
        is_WAN_compromised = false;
    end
    
    %% Baseline: Fails to detect DoS, uses stale ZOH signals
    U_base = Baseline_PI_AGC(X_baseline(:,k), is_WAN_compromised);
    dX_base = Plant_Dynamics_Nonlinear(X_baseline(:,k), d_k, U_base, params);
    X_baseline(:, k+1) = X_baseline(:,k) + Ts * dX_base;
    
    %% Proposed: DSC-sEC Framework
    if is_WAN_compromised
        % Triggers Local Edge Control (HO-CBF)
        U_prop = Controller_DSC_sEC_HOCBF(X_proposed(:,k), d_k, params);
    else
        % Nominal Operation
        U_prop = Baseline_PI_AGC(X_proposed(:,k), false); 
    end
    dX_prop = Plant_Dynamics_Nonlinear(X_proposed(:,k), d_k, U_prop, params);
    X_proposed(:, k+1) = X_proposed(:,k) + Ts * dX_prop;
end

disp('Simulation Complete. Run Evaluate_Performance_Metrics.m to view results.');
