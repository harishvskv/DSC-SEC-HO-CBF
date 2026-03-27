%% Main_Nonlinear_Microgrid_DoS.m
% Simulates the Cyber-Physical Microgrid under DoS Attack
% Demonstrates the transition from vulnerable Conventional-PI control to DSC-sEC.

% 1. Load Nonlinear System Parameters
params.M = input ('Enter Value, M:'); params.D = input ('Enter Value, D:'); params.gamma = input ('Enter Value, gamma:'); % Nonlinear load exponent
params.Tt = input ('Enter Value, Tt:'); params.Tg = input ('Enter Value, Tg:'); params.R = input ('Enter Value, R in p.u.:');
params.TB = input ('Enter Value, TB:'); params.Tac = input ('Enter Value, Tac:'); 
params.C_bat = input ('Enter Value, C_BESS:'); params.R_th = input ('Enter Value, RTh:'); params.C_th = input ('Enter Value, CTh:'); params.COP = input ('Enter Value, CoP:'); 

% 2. Simulation Setup
Tsim = input ('Enter Value:'); Ts = 0.1e-3; t = 0:Ts:Tsim; N = length(t);
X_baseline = zeros(8, N); % [df, Pm, Pv, PB, SoC, Tin, PAC, int_df]
X_proposed = zeros(8, N);

% 3. DoS Attack Definition (Breaching SC-ATD bound)
%t_attack_start = 115.0; %Deterministic case 
%t_attack_end = 118.0; 

% Stochastic Case
alpha = ones(1,N);
% Number of attacks (Poisson-like)
num_attacks = randi([1,3]);  
for k_attack = 1:num_attacks
    t_start = 80 + 120*rand();        % anywhere in simulation
    %duration = 0.3 + 3*rand();        % random duration
    duration = 2.0 + 1.4 * rand();
    alpha(t >= t_start & t <= t_start + duration) = 0;
end

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
