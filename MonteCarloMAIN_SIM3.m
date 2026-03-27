%% MonteCarlo_MAIN_SIM3.m
% 3-Case Comparative Study for Cyber-Resilient Microgrid

%% 1. Parameters 
params.gamma = input ('Enter Value, gamma:');
params.Tt = input ('Enter Value, Tt:');
params.Tg = input ('Enter Value, Tg:');
params.R = input ('Enter Value, R in p.u.:');
params.C_bat = input ('Enter Value, C_BESS:');
params.R_th = input ('Enter Value, RTh:'); 
params.C_th = input ('Enter Value, CTh:');
params.COP = input ('Enter Value, CoP:');
params.T_out = input ('Enter Value, TOut:');
params.Tac = input ('Enter Value, Tac:'); 
params.P_ref = input ('Enter Value, PRef:'); 
params.Tin_max = input ('Enter Value, Max. Indoor Temp:');  %ASHRAE Guideliness  

% Sensitivity Variations
params.M = input ('Enter Value, M:');
params.M = params.M * (1 + 0.2*randn());      % ±20% inertia
params.D = input ('Enter Value, D:');
params.D = params.D * (1 + 0.1*randn());    % ±10% damping
params.TB = input ('Enter Value, TB:');
params.TB = params.TB * (1 + 0.15*randn()); 

%% 2. Simulation Setup
Tsim = input ('Enter Value, TSimulation:'); Ts = 1e-3; t = 0:Ts:Tsim; N = length(t);

% Multi-Step Disturbances for Load, WT, and PV
P_Load_raw = ones(1, N); P_WT_raw = zeros(1, N); P_PV_raw = zeros(1, N);
for i = 1:N
    time = t(i);
    % Load Steps
    if time < 40, P_Load_raw(i) = 0.8;
    elseif time < 80, P_Load_raw(i) = 1.0;
    elseif time < 120, P_Load_raw(i) = 0.9;
    elseif time < 160, P_Load_raw(i) = 1.1;  %Overload condition %%Demonstartes the strngth of SEC-HO-CBF Overload and cyber attack at the same instant for deterministic case
    else, P_Load_raw(i) = 0.85; 
    end
    
    % Wind Turbine Steps
    if time < 60, P_WT_raw(i) = 0.8 * 0.15;
    elseif time < 100, P_WT_raw(i) = 0.6 * 0.15;
    elseif time < 140, P_WT_raw(i) = 0.9 * 0.15;
    elseif time < 180, P_WT_raw(i) = 0.7 * 0.15;
    else, P_WT_raw(i) = 1.0 * 0.15; 
    end
    
    % PV Steps
    if time < 20, P_PV_raw(i) = 0.6 * 0.15;
    elseif time < 100, P_PV_raw(i) = 0.8 * 0.15;
    elseif time < 150, P_PV_raw(i) = 0.5 * 0.15;
    else, P_PV_raw(i) = 0.9 * 0.15; 
    end
end

% Add multiplicative noise
    load_noise = 1 + 0.1*randn(1,N);   % ±10%
    pv_noise   = 1 + 0.15*randn(1,N);  % ±15%
    wind_noise = 1 + 0.2*randn(1,N);   % ±20%
    P_Load = P_Load_raw .* load_noise;
    P_PV   = P_PV_raw   .* pv_noise;
    P_WT = P_WT_raw .* wind_noise;

for i = 2:N
    P_load(i) = P_load(i-1) + (Ts/tau_f) * (P_Load_raw(i) - P_load(i-1));
    P_pv(i)   = P_pv(i-1) + (Ts/tau_f) * (P_PV_raw(i) - P_pv(i-1));
    P_wind(i) = P_wind(i-1) + (Ts/tau_f) * (P_WT_raw(i) - P_wind(i-1));
end
P_res = P_pv + P_wind; % Total Renewables 

alpha = ones(1,N);
% Number of attacks (Poisson-like)
num_attacks = randi([1,3]);  
for k_attack = 1:num_attacks
    t_start = 80 + 120*rand();        % anywhere in simulation
    %duration = 0.3 + 3*rand();        % random duration
    duration = 2.0 + 1.4 * rand();
    alpha(t >= t_start & t <= t_start + duration) = 0;
end

%% --- STOCHASTIC DoS ATTACK ---
alpha = ones(1,N);
num_attacks = randi([1,3]);
for k_attack = 1:num_attacks
    t_start = 115 + (125-115)*rand();   % your requirement
    duration = 2.0 + 1.4 * rand();
    alpha(t >= t_start & t <= t_start + duration) = 0;
end

%% --- STOCHASTIC COMMUNICATION DELAY ---
tau = 0.02 + 0.08*rand(1,N);   % base delay
% Add DoS-induced spikes (NOW alpha exists)
tau(alpha==0) = tau(alpha==0) + 0.5 + 1.5*rand(size(tau(alpha==0)));

%% 3. Initializations
%Initial Condition Uncertainty
   x0 = [0.01*randn(); 0.85; 0.85; 0; 0.5; 22 + randn(); 0.2; 0];
   h_M = 0.1; 
buf_len = round(h_M/Ts) + 10;

X_C1 = repmat(x0, 1, buf_len); % Case 1: PI (No Attack)
X_C2 = repmat(x0, 1, buf_len); % Case 2: PI (With Attack)
X_C3 = repmat(x0, 1, buf_len); % Case 3: Proposed (With Attack)

log_C1 = zeros(8,N); log_C2 = zeros(8,N); log_C3 = zeros(8,N);
log_C1(:,1)=x0; log_C2(:,1)=x0; log_C3(:,1)=x0;

UFLS_C1 = false; UFLS_C2 = false; UFLS_C3 = false;
UFLS_time_C2 = NaN; 

% Baseline PI Gains (System is stable under delay due to M=8) IEEE 39-bus three area and Kundur Case
Kp_base = 1.5; 
Ki_base = 0.8; 

disp('Running Corrected 3-Case Simulation...');

%% 4. Main Loop
for k = 1:N-1
    d_k = [P_res(k); P_load(k)];
    
    % Fetch delayed states based on tau(t)
    delay_steps = round(tau(k)/Ts);
    idx = max(1, size(X_C1,2) - delay_steps);
    x1_d = X_C1(:, idx); 
    x2_d = X_C2(:, idx); 
    x3_d = X_C3(:, idx); 
    
    %% --- CASE 1: Baseline PI (NO ATTACK) ---
    u1_c1 = -(Kp_base*x1_d(1) + Ki_base*x1_d(8));
    u2_c1 = 0.2 + (Kp_base*x1_d(1) + Ki_base*x1_d(8));
    U_C1 = [max(-0.4, min(0.4, u1_c1)); max(0, min(0.5, u2_c1))];
    
    %% --- CASE 2: Baseline PI (WITH DoS ATTACK) ---
    if alpha(k) == 1
        u1_c2 = -(Kp_base*x2_d(1) + Ki_base*x2_d(8));
        u2_c2 = 0.2 + (Kp_base*x2_d(1) + Ki_base*x2_d(8));
        U_C2 = [max(-0.4, min(0.4, u1_c2)); max(0, min(0.5, u2_c2))];
    else
        U_C2 = [0; 0]; 
    end
    
    %% --- CASE 3: Proposed CLF-HOCBF (WITH DoS ATTACK) ---
    if alpha(k) == 1
        [u1_c3, u2_c3] = CLF_HOCBF_CONTROLLER(x3_d, d_k, params);
        U_C3 = [u1_c3; u2_c3];
    else
        U_C3 = [0; 0]; 
    end
    
    %% Check Safety & UFLS Triggers
    if X_C1(1,end) <= -0.02 && ~UFLS_C1, UFLS_C1 = true; end
    if X_C2(1,end) <= -0.02 && ~UFLS_C2, UFLS_C2 = true; UFLS_time_C2 = t(k); end
    if X_C3(1,end) <= -0.02 && ~UFLS_C3, UFLS_C3 = true; end
    
    %% Integration Solver
    % C1
    [f1, g1] = microgrid_fg_dde(X_C1(:,end), d_k, params, UFLS_C1); k1 = f1 + g1*U_C1;
    [f2, g2] = microgrid_fg_dde(X_C1(:,end)+0.5*Ts*k1, d_k, params, UFLS_C1); k2 = f2 + g2*U_C1;
    [f3, g3] = microgrid_fg_dde(X_C1(:,end)+0.5*Ts*k2, d_k, params, UFLS_C1); k3 = f3 + g3*U_C1;
    [f4, g4] = microgrid_fg_dde(X_C1(:,end)+Ts*k3, d_k, params, UFLS_C1); k4 = f4 + g4*U_C1;
    X_C1 = [X_C1(:,2:end), X_C1(:,end) + (Ts/6)*(k1 + 2*k2 + 2*k3 + k4)];
    
    % C2
    [f1, g1] = microgrid_fg_dde(X_C2(:,end), d_k, params, UFLS_C2); k1 = f1 + g1*U_C2;
    [f2, g2] = microgrid_fg_dde(X_C2(:,end)+0.5*Ts*k1, d_k, params, UFLS_C2); k2 = f2 + g2*U_C2;
    [f3, g3] = microgrid_fg_dde(X_C2(:,end)+0.5*Ts*k2, d_k, params, UFLS_C2); k3 = f3 + g3*U_C2;
    [f4, g4] = microgrid_fg_dde(X_C2(:,end)+Ts*k3, d_k, params, UFLS_C2); k4 = f4 + g4*U_C2;
    X_C2 = [X_C2(:,2:end), X_C2(:,end) + (Ts/6)*(k1 + 2*k2 + 2*k3 + k4)];
    
    % C3
    [f1, g1] = microgrid_fg_dde(X_C3(:,end), d_k, params, UFLS_C3); k1 = f1 + g1*U_C3;
    [f2, g2] = microgrid_fg_dde(X_C3(:,end)+0.5*Ts*k1, d_k, params, UFLS_C3); k2 = f2 + g2*U_C3;
    [f3, g3] = microgrid_fg_dde(X_C3(:,end)+0.5*Ts*k2, d_k, params, UFLS_C3); k3 = f3 + g3*U_C3;
    [f4, g4] = microgrid_fg_dde(X_C3(:,end)+Ts*k3, d_k, params, UFLS_C3); k4 = f4 + g4*U_C3;
    X_C3 = [X_C3(:,2:end), X_C3(:,end) + (Ts/6)*(k1 + 2*k2 + 2*k3 + k4)];
    
    log_C1(:,k+1) = X_C1(:,end); log_C2(:,k+1) = X_C2(:,end); log_C3(:,k+1) = X_C3(:,end);
end

disp('Simulation Complete.');
