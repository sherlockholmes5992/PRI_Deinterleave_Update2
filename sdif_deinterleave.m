function [separated_sequences, extracted_flags] = sdif_deinterleave(pdwData, algoParams)
    fprintf('\n--- STARTING DEINTERLEAVING LOOP ---\n');
    
    TOA_sorted = pdwData.TOA_sorted;
    RF_sorted = pdwData.RF_sorted;
    PW_sorted = pdwData.PW_sorted;
    N_pulses = pdwData.N_pulses;

    C_max = algoParams.C_max;
    M0 = algoParams.M0;
    RF_0 = algoParams.RF_0;
    PW_0 = algoParams.PW_0;
    t_Bin = algoParams.t_Bin;
    x_emp = algoParams.x_emp;
    k_emp = algoParams.k_emp;

    extracted_flags = false(1, N_pulses);
    separated_sequences = {};
    C = 1;

    while C <= C_max
        active_idx = find(~extracted_flags);
        if length(active_idx) < 10
            break; 
        end 
        
        active_TOA = active_idx_subset(TOA_sorted, active_idx);
        active_RF  = active_idx_subset(RF_sorted, active_idx);
        active_PW  = active_idx_subset(PW_sorted, active_idx);
        
        if length(active_TOA) <= C
            break;
        end
        
        % Tính toán Histogram m?c sai phân C
        diff_TOA = abs(active_TOA(1 : end-C) - active_TOA(1+C : end));
        diff_TOA = round(diff_TOA, 7); 
        max_diff = max(diff_TOA);
        edges = 0 : t_Bin : (max_diff + t_Bin);
        [N_counts, ~] = histcounts(diff_TOA, edges);
        bin_centers = edges(1:end-1) + t_Bin/2; 
        
        E_active = length(active_TOA);
        N_bins = length(bin_centers);
        tau = 1:N_bins; 
        Threshold = x_emp * (E_active - C) * exp(-tau * k_emp * N_bins);
        
        pot_idx = find(N_counts > Threshold);
        pot_PRIs_raw = bin_centers(pot_idx);
        
        if isempty(pot_PRIs_raw)
            C = C + 1; 
            continue;
        end
        
        % Lo?i b? sóng hài v?t lý trên ?? th?
        is_harmonic = false(1, length(pot_PRIs_raw));
        for i = 1:length(pot_PRIs_raw)
            for j = 1:(i-1)
                if ~is_harmonic(j)
                    ratio = pot_PRIs_raw(i) / pot_PRIs_raw(j);
                    if abs(ratio - round(ratio)) < (t_Bin / pot_PRIs_raw(j))
                        is_harmonic(i) = true;
                        break;
                    end
                end
            end
        end
        pot_PRIs = pot_PRIs_raw(~is_harmonic);
        pot_idx_filtered = pot_idx(~is_harmonic);
        
        % Phân bi?t Jittered và Non-Jittered ?ng viên d?a trên ?? li?n k? vách ng?n
        is_jittered = false(1, length(pot_idx_filtered));
        for i = 1:length(pot_idx_filtered)
            if i < length(pot_idx_filtered) && (pot_idx_filtered(i+1) - pot_idx_filtered(i) <= 2)
                is_jittered(i) = true;
                is_jittered(i+1) = true;
            end
        end
        
        jittered_PRIs = pot_PRIs(is_jittered);
        non_jittered_PRIs = pot_PRIs(~is_jittered);
        all_pot_PRIs = [non_jittered_PRIs, jittered_PRIs];
        
        % In nhanh tr?ng thái ?? ph?c v? g? l?i tr?c quan
        fprintf('\n[M?c sai phân C=%d] Các PRI ti?m n?ng (PRI_p) v??t ng??ng:\n', C);
        if ~isempty(non_jittered_PRIs), fprintf('  + Non-Jittered PRIs (ms): %s\n', num2str(non_jittered_PRIs * 1e3, '%.4f  ')); end
        if ~isempty(jittered_PRIs),     fprintf('  + Jittered PRIs (ms):     %s\n', num2str(jittered_PRIs * 1e3, '%.4f  ')); end
        
        success_in_this_C = false;
        
        % Quét bám ?uôi chu?i xung ?a tham s? (Multi-Parameter Sequence Search)
        for i = 1:length(all_pot_PRIs)
            current_PRI = all_pot_PRIs(i);
            W = iif(ismember(current_PRI, jittered_PRIs), current_PRI * 0.35, 2e-5);
            
            p = 1; 
            while p <= length(active_TOA)
                k = p;
                temp_seq_idx_local = k; 
                pulses_found = 1;
                
                % Giai ?o?n 1: Tìm ki?m móng xích (4 xung ??u)
                while pulses_found < 4
                    found_next = false;
                    for j = (k+1):length(active_TOA)
                        TOA_diff = active_TOA(j) - active_TOA(k);
                        Z_k = abs(current_PRI - TOA_diff);
                        M_k = ((active_RF(j) - active_RF(k))^2) / (RF_0^2) + ...
                              ((active_PW(j) - active_PW(k))^2) / (PW_0^2);
                        
                        if Z_k < W && M_k < M0
                            temp_seq_idx_local = [temp_seq_idx_local, j];
                            k = j; 
                            pulses_found = pulses_found + 1;
                            found_next = true;
                            break; 
                        end
                    end
                    if ~found_next, break; end 
                end
                
                % Giai ?o?n 2: Bám ?uôi ch?p nh?n b? qua xung r?i r?ng (Missing Pulses Engine)
                if pulses_found >= 4
                    M_multiplier = 1; 
                    while k < length(active_TOA)
                        found_next = false;
                        for j = (k+1):length(active_TOA)
                            TOA_diff = active_TOA(j) - active_TOA(k);
                            Z_k = abs(current_PRI * M_multiplier - TOA_diff);
                            M_k = ((active_RF(j) - active_RF(k))^2) / (RF_0^2) + ...
                                  ((active_PW(j) - active_PW(k))^2) / (PW_0^2);
                            
                            if Z_k < (W * M_multiplier) && M_k < M0
                                temp_seq_idx_local = [temp_seq_idx_local, j];
                                k = j;
                                M_multiplier = 1; 
                                found_next = true;
                                break;
                            end
                        end
                        
                        if ~found_next
                            M_multiplier = M_multiplier + 1;
                            if M_multiplier > 15, break; end % C?t chu?i n?u m?t quá 15 xung liên ti?p
                        end
                    end
                    
                    % Ki?m tra tiêu chu?n ?? dài chu?i xung t?i thi?u
                    if length(temp_seq_idx_local) >= 10
                        global_indices = active_idx(temp_seq_idx_local);
                        extracted_flags(global_indices) = true;
                        separated_sequences{end+1} = global_indices;
                        
                        fprintf('>> Successfully extracted 1 Emitter! (C=%d, PRI=%.4f ms, Pulses=%d)\n', ...
                                C, current_PRI * 1e3, length(global_indices));
                        
                        success_in_this_C = true;
                        break; 
                    end
                end
                p = p + 1;
            end
            if success_in_this_C, break; end
        end
        C = iif(success_in_this_C, 1, C + 1); % Quy lu?t Hasani: Thành công thì reset C v? 1, th?t b?i t?ng b?c C
    end
end

% --- Hàm b? tr? c?c b? ?? vi?t code ng?n g?n ---
function val = active_idx_subset(arr, idx), val = arr(idx); end
function res = iif(cond, true_val, false_val), if cond, res = true_val; else, res = false_val; end; end