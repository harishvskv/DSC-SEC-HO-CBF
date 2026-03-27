function dX = Plant_Dynamics_Nonlinear(x, d, U, params)
    % Nonlinear CPMS State-Space Dynamics
    % x = [df, Pm, Pv, PB, SoC, Tin, PAC, int_dw]
    
    df = x(1); Pm = x(2); Pv = x(3); PB = x(4); 
    SoC = x(5); Tin = x(6); PAC = x(7); int_df = x(8);
    
    P_res = d(1); P_L0 = d(2);
    u_BESS = U(1); u_AIAC = U(2);
    
    dX = zeros(8,1);
    
    % 1. Nonlinear Swing Equation with Voltage/Freq Dependent Load
    df_bound = max(df, -0.5); % singularity
    P_load = P_L0 * (1 + params.D * df_bound)^params.gamma;
    dX(1) = (1/params.M) * (Pm + P_res + PB - PAC - P_load);
    
    % 2. Governor & Turbine Dynamics
    dX(2) = (-Pm + Pv) / params.Tt;
    dX(3) = (0.85 - Pv - (1/params.R)*df) / params.Tg;
    
    % 3. BESS Dynamics
    dX(4) = -(1/params.TB) * PB + (1/params.TB) * u_BESS;
    dX(5) = -PB / params.C_bat;
    
    % 4. AIAC Smart Building Thermodynamics
    dX(6) = -params.COP/params.C_th * PAC - (Tin - 32)/(params.R_th * params.C_th);
    dX(7) = -(1/params.Tac) * PAC + (1/params.Tac) * u_AIAC;
    
    % 5. Integral Tracking Error
    dX(8) = df;
end
