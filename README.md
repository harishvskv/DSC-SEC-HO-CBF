# DSC-SEC-HO-CBF
Safety-Critical Edge-Control for Cyber-Physical Microgrids Under DoS Attacks: A Time-Delay High-Order CBF Approach
This repository contains the core MATLAB implementation of the DSC-sEC framework for Cyber-Physical Microgrid Systems (CPMS) under Denial-of-Service (DoS) attacks, as presented in Manuscript ID	TSG-00739-2026.
To protect proprietary research extensions and multi-area grid topologies, this repository provides a **representative nonlinear isolated microgrid model**. It demonstrates the fundamental implementation of:
1. The nonlinear delay-aware plant dynamics.
2. The High-Order Control Barrier Function (HO-CBF) Quadratic Program.
3. The stochastic Monte Carlo robustness validation framework.

In order to demonstrate the scalabality and universatality, you may enter the paper's values or the values of your own interest. please note that the SEC-HO-CBF Frequency/ACE regulatore works optimally well for low inertia values aswell.

### Files Included:
* `Main_Nonlinear_Microgrid_DoS.m`: The primary simulation script executing the switching logic between WAN (baseline) and Edge (DSC-sEC) under DoS.
* `Plant_Dynamics_Nonlinear.m`: The highly nonlinear state-space representation of the CPMS (Swing equations, BESS, AIAC thermal dynamics).
* `Controller_DSC_sEC_HOCBF.m`: The decentralized edge-controller solving the HO-CBF QP in real-time.
* `MonteCarlo_Robustness_Test.m`: The 1000-run stochastic uncertainty wrapper.
* `Evaluate_Performance_Metrics.m`: Computes IAE, ISE, ITAE, and thermal comfort indices.
