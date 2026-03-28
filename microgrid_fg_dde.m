function [f,g] = microgrid_fg_dde(x, d, params, UFLS_flag)
    % STATES
    df = x(1); Pm = x(2); Pv = x(3); PB = x(4); S = x(5); Tin = x(6); PAC = x(7); int_df = x(8);
    
    % DISTURBANCES
    Pres = d(1); 
    P_L0 = d(2);
    
    % Apply UFLS if triggered (shed 10% of base load)
    if UFLS_flag
        P_L0 = P_L0 - 0.10;
    end
    
    f = zeros(8,1);
    g = zeros(8,2);
    
    % Fixed typo: changed dw to df
    df_bound = max(df, -0.5);
    P_load = P_L0 * (1 + params.D * df_bound)^params.gamma;
    
    %% Swing Equation
    f(1) = (1/params.M) * (Pm + Pres + PB - PAC - P_load);
    
       
    %% Governor Dead-Band (GDB) Filter
    GDB_range = 0.0006;           % ±0.0006 Hz dead band (NERC)
    
    if abs(df) <= GDB_range
        df_gdb = 0;               % Inside dead-band: no governor action
    elseif df > GDB_range
        df_gdb = df - GDB_range;  % Positive deviation beyond dead-band
    else
        df_gdb = df + GDB_range;  % Negative deviation beyond dead-band
    end

    %% Governor (Primary + AGC interaction)
    Pref = params.P_ref; 
    % Author Tried using df_gdb instead of df so the governor ignores tiny 
variations
    f(3) = (Pref - Pv - (1/params.R)*df_gdb) / params.Tg;
    
    %% Generation Rate Constraint (GRC) & Turbine
    GRC_limit = 0.0017;           % pu/s generation rate constraint
    
    % 1. unbounded rate of change of mechanical power
    dPm_dt_unconstrained = (-Pm + Pv) / params.Tt;
    
    % 2. Clamp the rate of change using max/min limits
    f(2) = max(-GRC_limit, min(GRC_limit, dPm_dt_unconstrained));

    % =========================================================
    
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
    f(8) = df; 
end
