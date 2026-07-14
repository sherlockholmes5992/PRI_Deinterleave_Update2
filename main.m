% =========================================================================
% PROJECT: Enhanced Radar Deinterleaving via Modular SDIF & Fusing PDWs
% Based on the paper: Hasani & Khosravi (2021)
% Architecture: Modular Design for Optimization and Easy Debugging
% =========================================================================

clc; clear; close all;

%% --- B??C 1: C?U HÌNH THAM S? ?ÀI PHÁT (EMITTER CONFIGURATION) ---
t_sim = 0.12; % Th?i gian mô ph?ng t?ng th? (120 ms)

clear emCfg; 
emCfg(1).type = 'Fixed';     emCfg(1).PRI = 157e-5;             emCfg(1).RF = 3e9;   emCfg(1).PW = 2e-6;   emCfg(1).jitter = 0;   emCfg(1).t_start = 0;      emCfg(1).p_missing = 0.15;          emCfg(1).toa_error = 2e-6;
emCfg(2).type = 'Fixed';     emCfg(2).PRI = 157e-5;             emCfg(2).RF = 3e9;   emCfg(2).PW = 2e-6;   emCfg(2).jitter = 0;   emCfg(2).t_start = 0.0015; emCfg(2).p_missing = 0.15;          emCfg(2).toa_error = 2e-6;
emCfg(3).type = 'Staggered'; emCfg(3).PRI = [106e-5, 112e-5, 90.8e-5]; emCfg(3).RF = 10e9;  emCfg(3).PW = 6e-6;   emCfg(3).jitter = 0;   emCfg(3).t_start = 0.0005; emCfg(3).p_missing = 0.15;   emCfg(3).toa_error = 2e-6;
emCfg(4).type = 'Jittered';  emCfg(4).PRI = 10e-4;              emCfg(4).RF = 17e9;  emCfg(4).PW = 9e-6;   emCfg(4).jitter = 0.1; emCfg(4).t_start = 0.0030; emCfg(4).p_missing = 0.15;          emCfg(4).toa_error = 2e-6;
emCfg(5).type = 'Fixed';     emCfg(5).PRI = 157e-5;             emCfg(5).RF = 3e9;   emCfg(5).PW = 2.4e-6; emCfg(5).jitter = 0;   emCfg(5).t_start = 0.0035; emCfg(5).p_missing = 0.15;          emCfg(5).toa_error = 2e-6;
emCfg(6).type = 'Fixed';     emCfg(6).PRI = 200e-5;             emCfg(6).RF = 5e9;   emCfg(6).PW = 4e-6;   emCfg(6).jitter = 0;   emCfg(6).t_start = 0.0045; emCfg(6).p_missing = 0.15;          emCfg(6).toa_error = 2e-6;
emCfg(6).type = 'Fixed';     emCfg(6).PRI = 200e-5;             emCfg(6).RF = 5e9;   emCfg(6).PW = 4e-6;   emCfg(6).jitter = 0;   emCfg(6).t_start = 0.0045; emCfg(6).p_missing = 0.15;          emCfg(7).toa_error = 2e-6;
%% --- B??C 2: C?U HÌNH THAM S? THU?T TOÁN (ALGORITHM HYPERPARAMETERS) ---
algoParams.C_max = 30;       % Gi?i h?n m?c sai phân t?i ?a ?? tránh l?p vô h?n
algoParams.M0 = 1.0;         % Ng??ng t??ng ??ng kho?ng cách vân tay v?t lý (RF, PW)
algoParams.RF_0 = 3e9;       % Bán kính chu?n hóa t?n s? vô tuy?n
algoParams.PW_0 = 1e-6;      % Bán kính chu?n hóa ?? r?ng xung
algoParams.t_Bin = 1e-5;     % ?? r?ng c?a m?t ô ??m (Bin) trên Histogram
algoParams.x_emp = 0.5;      % H? s? th?c nghi?m tính hàm ng??ng t?i ?u
algoParams.k_emp = 5;        % H? s? suy gi?m c?a ???ng cong ng??ng ??ng

%% --- B??C 3: KÍCH HO?T ??NG C? MÔ PH?NG & TÁCH CHU?I ---

% 1. Sinh lu?ng d? li?u xung ?an xen b? nhi?u suy hao th?c t?
[pdwData, initialSDIF] = generate_pulses(emCfg, t_sim);

% 2. Ch?y thu?t toán vòng l?p sai phân SDIF & Tìm ki?m ?a tham s? g?c
[separated_sequences, extracted_flags] = sdif_deinterleave(pdwData, algoParams);

% 3. Kh?i h?u x? lý cao c?p: G?p các m?nh sóng hài và b?o v? ?ài ??ng kênh
mergeParams.rf_merge_tol = 0.5e9;   % Dung sai th?t ch?t vân tay RF (0.5 GHz)
mergeParams.pw_merge_tol = 0.5e-6;  % Dung sai th?t ch?t vân tay PW (0.5 us)
mergeParams.pri_ratio_tol = 0.08;   % Sai s? t? l? sóng hài cho phép (8%)
final_sequences = merge_sequences(pdwData, separated_sequences, mergeParams);

%% --- B??C 4: XU?T BÁO CÁO KI?M TOÁN & ?? TH? TR?C QUAN ---
export_results(pdwData, final_sequences, extracted_flags);
plot_results(pdwData, final_sequences, initialSDIF, length(emCfg));