function U_edge = Controller_DSC_sEC_HOCBF(x, d, params)
    % Solves the High-Order CBF Quadratic Program for Edge Control
    
    dw = x(1); Pm = x(2); Pv = x(3); PB = x(4); PAC = x(7);
    P_res = d(1); P_L0 = d(2);

    %% 1. Nonlinear Relative Degree Formulation
    dw_bound = max(dw, -0.5);
    P_load = P_L0 * (1 + params.D * dw_bound)^params.gamma;
    f1 = (1/params.M) * (Pm + P_res + PB - PAC - P_load); % dw_dot

    Pm_dot = (-Pm + Pv) / params.Tt;
    PB_dot = -(1/params.TB) * PB;
    PAC_dot = -(1/params.Tac) * PAC;
    Pload_dot = params.D * params.gamma * P_L0 * (1 + params.D * dw_bound)^(params.gamma - 1) * f1;

    % Lie Derivatives
    Lf2h = (1/params.M) * (Pm_dot + PB_dot - PAC_dot - Pload_dot);
    LgLfh = [1/(params.M * params.TB), -1/(params.M * params.Tac)];

    %% 2. Barrier Definition
    dw_min = -0.02; % Safety Limit (pu)
    h = dw - dw_min;
    
    % Class-K parameter tuning
    cb1 = 25; cb2 = 25;
    
    A_cbf = -LgLfh;
    B_cbf = Lf2h + (cb1 + cb2) * f1 + cb1 * cb2 * h;

    %% 3. Quadratic Program Formulation
    % Optimize BESS and AIAC dispatch while penalizing thermodynamic deviation
    H = diag([10, 500]); % Heavily penalize AIAC usage to preserve comfort
    f_qp = [0; 0];
    
    % Physical Actuator Limits
    lb = [-0.6; -0.6];
    ub = [0.6; 0.6];
    
    % Solve QP
    options = optimoptions('quadprog', 'Display', 'off');
    [U_edge, ~, exitflag] = quadprog(H, f_qp, A_cbf, B_cbf, [], [], lb, ub, [], options);
    
    if exitflag ~= 1
        U_edge = [0;0]; % Fallback if infeasible
    end
end
