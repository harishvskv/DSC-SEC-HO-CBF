function [u1, u2] = CLF_HOCBF_CONTROLLER(x, d, params)
    dw = x(1); Pm = x(2); Pv = x(3); PB = x(4); S = x(5); Tin = x(6); PAC = x(7); int_dw = x(8);
    Pres = d(1); P_L0 = d(2);

    %% 1. Grid Dynamics & Relative Degree 2 Variables
    dw_bound = max(dw, -0.5);
    P_load = P_L0 * (1 + params.D * dw_bound)^params.gamma;

    f1 = (1/params.M) * (Pm + Pres + PB - PAC - P_load); % dw_dot
    Pm_dot = (-Pm + Pv) / params.Tt;
    PB_dot = -(1/params.TB) * PB;
    PAC_dot = -(1/params.Tac) * PAC;
    Pload_dot = params.D * params.gamma * P_L0 * (1 + params.D * dw_bound)^(params.gamma - 1) * f1;

    Lf2h = (1/params.M) * (Pm_dot + PB_dot - PAC_dot - Pload_dot);
    LgLfh = [1/(params.M * params.TB), -1/(params.M * params.Tac)];

    %% 2. Reference Controller (PI + Virtual Inertia)
    Kp = 1.5; Ki = 0.8; K_vi = 0.6;
    u1_ref = -(Kp * dw + Ki * int_dw + K_vi * f1); 
    u2_ref = 0.2 + (Kp * dw + Ki * int_dw); 

    %% 3. HOCBF (Safety)
    dw_min = -0.02;
    h = dw - dw_min;
    cb1 = 8; cb2 = 8;
    A_cbf = -LgLfh;
    B_cbf = Lf2h + (cb1 + cb2) * f1 + cb1 * cb2 * h;

    %% 4. Backstepping CLF (Stability)
    c1 = 4; c2 = 4;
    z1 = dw; z2 = f1 + c1 * z1;
    V = 0.5 * z1^2 + 0.5 * z2^2;
    A_clf = [z2 * LgLfh(1), z2 * LgLfh(2), -1];
    B_clf = -c2 * V + c1 * z1^2 - z1 * z2 - z2 * (Lf2h + c1 * f1);

    %% 5. QP Solver
    A = [A_cbf, 0; A_clf]; b = [B_cbf; B_clf];
    H = diag([2, 2, 50000]); fqp = [-2 * u1_ref; -2 * u2_ref; 0];
    lb = [-0.4; 0; 0]; ub = [0.4; 0.5; Inf];

    options = optimoptions('quadprog', 'Display', 'off');
    [U, ~, flag] = quadprog(H, fqp, A, b, [], [], lb, ub, [], options);

    if flag ~= 1
        U = [max(min(u1_ref, ub(1)), lb(1)); max(min(u2_ref, ub(2)), lb(2)); 0]; 
    end
    u1 = U(1); u2 = U(2);
end