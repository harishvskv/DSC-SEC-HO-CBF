function [f,g] = microgrid_fg_dde(x, d, params, UFLS_flag)
    % STATES
    dw = x(1); Pm = x(2); Pv = x(3); PB = x(4); S = x(5); Tin = x(6); PAC = x(7); int_dw = x(8);
    
    % DISTURBANCES
    Pres = d(1); 
    P_L0 = d(2);
    
    % Apply UFLS if triggered (shed 18% of base load)
    if UFLS_flag
        P_L0 = P_L0 - 0.18;
    end
    
    f = zeros(8,1);
    g = zeros(8,2);
    
    dw_bound = max(dw, -0.5);
    P_load = P_L0 * (1 + params.D * dw_bound)^params.gamma;
    
    %% Swing Equation
    f(1) = (1/params.M) * (Pm + Pres + PB - PAC - P_load);
    
    %% Turbine
    f(2) = (-Pm + Pv) / params.Tt;
    
    %% Governor (Primary + AGC interaction)
    Pref = params.P_ref; 
    f(3) = (Pref - Pv - (1/params.R)*dw) / params.Tg;
    
    %% Battery Power
    f(4) = -(1/params.TB) * PB;
    g(4,1) = 1/params.TB;
    
    %% Battery SOC
    f(5) = -(1/params.C_bat) * PB;
    
    %% Indoor Temperature
    f(6) = -(1/(params.R_th*params.C_th)) * (Tin - params.T_out) - (params.COP/params.C_th) * PAC;
    
    %% AC Power
    f(7) = -(1/params.Tac) * PAC;
    g(7,2) = 1/params.Tac;
    
    %% Area Control Error (Integral for AGC)
    f(8) = dw; 
end