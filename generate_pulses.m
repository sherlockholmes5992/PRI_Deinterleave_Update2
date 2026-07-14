function [pdwData, initialSDIF] = generate_pulses(emCfg, t_sim)
    num_emitters = length(emCfg);
    TOA_all = []; RF_all = []; PW_all = []; ID_all = [];
    total_theoretical_pulses = 0;

    for i = 1:num_emitters
        t_current = emCfg(i).t_start;
        TOA_temp = [];
        
        % T?o chu?i TOA g?c lý thuy?t theo t?ng lo?i hình phát
        if strcmp(emCfg(i).type, 'Fixed')
            TOA_temp = emCfg(i).t_start : emCfg(i).PRI : t_sim;
            
        elseif strcmp(emCfg(i).type, 'Staggered')
            idx = 1; levels = emCfg(i).PRI; num_levels = length(levels);
            while t_current <= t_sim
                TOA_temp = [TOA_temp, t_current];
                t_current = t_current + levels(mod(idx-1, num_levels) + 1);
                idx = idx + 1;
            end
            
        elseif strcmp(emCfg(i).type, 'Jittered')
            while t_current <= t_sim
                TOA_temp = [TOA_temp, t_current];
                jitter_val = emCfg(i).PRI * (1 + emCfg(i).jitter * (2*rand - 1));
                t_current = t_current + jitter_val;
            end
        end
        
        N_theoretical = length(TOA_temp);
        total_theoretical_pulses = total_theoretical_pulses + N_theoretical;
        
        % Áp d?ng c? ch? Pulse Dropout Mask (M?t xung th?c t?)
        keep_mask = rand(1, N_theoretical) > emCfg(i).p_missing;
        TOA_temp = TOA_temp(keep_mask);
        N_actual = length(TOA_temp);
        
        if N_actual > 0
            TOA_noise = emCfg(i).toa_error * randn(1, N_actual);
            TOA_temp = TOA_temp + TOA_noise;
            
            ID_temp = ones(1, N_actual) * i;
            % Cài ??t nhi?u ?o ??c v?t lý cho RF và PW
            RF_temp = emCfg(i).RF * (1 + 0.1 * (2*rand(1, N_actual) - 1)); 
            PW_temp = emCfg(i).PW * (1 + 0.05 * randn(1, N_actual));      
            
            TOA_all = [TOA_all, TOA_temp];
            RF_all  = [RF_all, RF_temp];
            PW_all  = [PW_all, PW_temp];
            ID_all  = [ID_all, ID_temp];
        end
    end

    % S?p x?p tr?n lu?ng theo trình t? th?i gian TOA ??n máy thu
    [TOA_sorted, sort_idx] = sort(TOA_all);
    N_pulses  = length(TOA_sorted);
    RF_sorted = RF_all(sort_idx);
    PW_sorted = PW_all(sort_idx);
    True_ID   = ID_all(sort_idx); 

    % ?óng gói d? li?u ra
    pdwData.TOA_sorted = TOA_sorted;
    pdwData.RF_sorted = RF_sorted;
    pdwData.PW_sorted = PW_sorted;
    pdwData.True_ID = True_ID;
    pdwData.N_pulses = N_pulses;
    pdwData.total_theoretical_pulses = total_theoretical_pulses;

    % L?u tr? tr?ng thái Histogram t?i C=1 ph?c v? v? ?? th? ki?m tra s? 1
    diff_TOA_C1 = abs(TOA_sorted(1 : end-1) - TOA_sorted(2 : end));
    t_Bin_C1 = 1e-5; 
    max_diff_C1 = max(diff_TOA_C1);
    edges_C1 = 0 : t_Bin_C1 : (max_diff_C1 + t_Bin_C1);
    [N_counts_C1, ~] = histcounts(diff_TOA_C1, edges_C1);
    bin_centers_C1 = edges_C1(1:end-1) + t_Bin_C1/2; 

    x_emp_C1 = 0.5; k_emp_C1 = 0.1 / length(bin_centers_C1); 
    tau_C1 = 1:length(bin_centers_C1); 
    Threshold_C1 = x_emp_C1 * (N_pulses - 1) * exp(-tau_C1 * k_emp_C1 * length(bin_centers_C1));
    pot_idx_C1 = find(N_counts_C1 > Threshold_C1);

    initialSDIF.bin_centers_C1 = bin_centers_C1;
    initialSDIF.N_counts_C1 = N_counts_C1;
    initialSDIF.Threshold_C1 = Threshold_C1;
    initialSDIF.pot_idx_C1 = pot_idx_C1;
    initialSDIF.max_diff_C1 = max_diff_C1;
end