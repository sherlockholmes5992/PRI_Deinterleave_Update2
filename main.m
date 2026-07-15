% =========================================================================
% PROJECT: Enhanced Radar Deinterleaving via Modular SDIF & Fusing PDWs
% Based on the paper: Hasani & Khosravi (2021)
% Architecture: Modular Design for Optimization and Easy Debugging
% =========================================================================

clc; clear; close all;

%% --- STEP 1: EMITTER CONFIGURATION ---
t_sim = 0.15; % Total simulation time window (120 ms)

clear emCfg; 
emCfg(1).type = 'Fixed';     emCfg(1).PRI = 157e-5;             emCfg(1).RF = 5e9;   emCfg(1).PW = 2e-6;   emCfg(1).jitter = 0;   emCfg(1).t_start = 0;      emCfg(1).p_missing = 0.1;          emCfg(1).toa_error = 2e-6;
emCfg(2).type = 'Fixed';     emCfg(2).PRI = 157e-5;             emCfg(2).RF = 5e9;   emCfg(2).PW = 2e-6;   emCfg(2).jitter = 0;   emCfg(2).t_start = 0.0015; emCfg(2).p_missing = 0.1;          emCfg(2).toa_error = 2e-6;
emCfg(3).type = 'Staggered'; emCfg(3).PRI = [106e-5, 112e-5, 90.8e-5]; emCfg(3).RF = 5e9;  emCfg(3).PW = 6e-6;   emCfg(3).jitter = 0;   emCfg(3).t_start = 0.0005; emCfg(3).p_missing = 0.1;          emCfg(3).toa_error = 2e-6;
emCfg(4).type = 'Jittered';  emCfg(4).PRI = 10e-4;              emCfg(4).RF = 5e9;  emCfg(4).PW = 9e-6;   emCfg(4).jitter = 0.1; emCfg(4).t_start = 0.0030; emCfg(4).p_missing = 0.1;          emCfg(4).toa_error = 2e-6;
emCfg(5).type = 'Fixed';     emCfg(5).PRI = 157e-5;             emCfg(5).RF = 5e9;   emCfg(5).PW = 2.4e-6; emCfg(5).jitter = 0;   emCfg(5).t_start = 0.0035; emCfg(5).p_missing = 0.1;          emCfg(5).toa_error = 2e-6;
emCfg(6).type = 'Fixed';     emCfg(6).PRI = 200e-5;             emCfg(6).RF = 5e9;   emCfg(6).PW = 4.2e-6;   emCfg(6).jitter = 0;   emCfg(6).t_start = 0.0045; emCfg(6).p_missing = 0.1;          emCfg(6).toa_error = 2e-6;
emCfg(7).type = 'Fixed';     emCfg(7).PRI = 201e-5;             emCfg(7).RF = 5e9;   emCfg(7).PW = 4e-6;   emCfg(7).jitter = 0;   emCfg(7).t_start = 0.0005; emCfg(7).p_missing = 0.1;          emCfg(7).toa_error = 2e-6;
emCfg(8).type = 'Staggered';     emCfg(8).PRI = [50e-5, 75e-5, 90e-5];             emCfg(8).RF = 5e9;   emCfg(8).PW = 4.5e-6;   emCfg(8).jitter = 0;   emCfg(8).t_start = 0.0022; emCfg(8).p_missing = 0.1;          emCfg(8).toa_error = 2e-6;
%% --- STEP 2: ALGORITHM HYPERPARAMETERS ---
algoParams.C_max = 25;       % Maximum difference level (C-level) to prevent infinite loops
algoParams.M0 = 1;         % Distance threshold for physical radar fingerprints matching (RF, PW)
algoParams.RF_0 = 3e9;       % Radio Frequency (RF) normalization radius
algoParams.PW_0 = 1e-6;      % Pulse Width (PW) normalization radius
algoParams.t_Bin = 1e-5;     % Time bin width for the Sequential Difference Histogram (SDIF)
algoParams.x_emp = 0.08;      % Empirical scaling factor for the optimal threshold curve bounds
algoParams.k_emp = 1;        % Decay coefficient controlling the dynamic threshold slope curvature

%% --- STEP 3: EXECUTE SIMULATION ENGINE & DEINTERLEAVING ---

% 1. Generate interleaved radar pulse stream with realistic channel degradation and dropouts
[pdwData, initialSDIF] = generate_pulses(emCfg, t_sim);

% 2. Execute the core SDIF difference loops and multi-parameter sequence search engine
[separated_sequences, extracted_flags] = sdif_deinterleave(pdwData, algoParams);

% % 3. Advanced Post-Processing: Recover harmonic fragments and safeguard co-channel paths
% mergeParams.rf_merge_tol = 0.5e9;   % Tight RF fingerprint matching boundary tolerance (0.5 GHz)
% mergeParams.pw_merge_tol = 0.5e-6;  % Tight PW fingerprint matching boundary tolerance (0.5 us)
% mergeParams.pri_ratio_tol = 0.08;   % Allowed fundamental/harmonic ratio deviation tolerance (8%)
% final_sequences = merge_sequences(pdwData, separated_sequences, mergeParams);
% 
% %% --- STEP 4: GENERATE AUDIT REPORT & VISUAL PLOTS ---
% export_results(pdwData, final_sequences, extracted_flags);
% plot_results(pdwData, final_sequences, initialSDIF, length(emCfg));

%% --- STEP 4: GENERATE AUDIT REPORT & VISUAL PLOTS ---
export_results(pdwData, separated_sequences, extracted_flags);
plot_results(pdwData, separated_sequences, initialSDIF, length(emCfg));