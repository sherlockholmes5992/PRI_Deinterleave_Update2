function [separated_sequences, extracted_flags] = sdif_deinterleave(pdwData, algoParams)
    fprintf('\n--- STARTING DEINTERLEAVING LOOP ---\n');
    
    % --- Extract packed Pulse Description Word (PDW) parameters ---
    TOA_sorted = pdwData.TOA_sorted;
    RF_sorted = pdwData.RF_sorted;
    PW_sorted = pdwData.PW_sorted;
    N_pulses = pdwData.N_pulses;

    % --- Map algorithm hyperparameters ---
    C_max = algoParams.C_max;
    M0 = algoParams.M0;
    RF_0 = algoParams.RF_0;
    PW_0 = algoParams.PW_0;
    t_Bin = algoParams.t_Bin;
    x_emp = algoParams.x_emp;
    k_emp = algoParams.k_emp;

    % --- Initialize tracking arrays ---
    extracted_flags = false(1, N_pulses); 
    separated_sequences = {};             
    C = 1;                                

    while C <= C_max
        % L?c ra các xung ch?a ???c trích xu?t (Active pulses)
        active_idx = find(~extracted_flags);
        if length(active_idx) < 10 
            break; 
        end 
        
        % ??ng b? hóa lu?ng xung ??u vào Xin theo ?úng ??nh ngh?a bài báo (M?c 3)
        Xin.TOA = TOA_sorted(active_idx);
        Xin.RF  = RF_sorted(active_idx);
        Xin.PW  = PW_sorted(active_idx);
        
        if length(Xin.TOA) <= C
            break; 
        end
        
        % --- COMPUTE SDIF HISTOGRAM AT DIFFERENCE LEVEL C ---
        diff_TOA = abs(Xin.TOA(1 : end-C) - Xin.TOA(1+C : end));
        diff_TOA = round(diff_TOA, 7); 
        max_diff = max(diff_TOA);
        edges = 0 : t_Bin : (max_diff + t_Bin);
        [N_counts, ~] = histcounts(diff_TOA, edges);
        bin_centers = edges(1:end-1) + t_Bin/2; 
        
        % --- CALCULATE DYNAMIC OPTIMAL THRESHOLD CURVE ---
        E_active = length(Xin.TOA);
        N_bins = length(bin_centers);
        tau = 1:N_bins; 
        Threshold = x_emp * (E_active - C) * exp(-tau ./ (k_emp * N_bins));
        
        % --- EXTRACT POTENTIAL RADAR PRI CANDIDATES {PRI_p} ---
        pot_idx = find(N_counts > Threshold);
        pot_PRIs_raw = bin_centers(pot_idx);
        
        if isempty(pot_PRIs_raw)
            C = C + 1; 
            continue;
        end
        
        % --- HARMONIC SUPPRESSION FILTER ---
        is_harmonic = false(1, length(pot_PRIs_raw));
        for i = 1:length(pot_PRIs_raw)
            for j_harm = 1:(i-1)
                if ~is_harmonic(j_harm)
                    ratio = pot_PRIs_raw(i) / pot_PRIs_raw(j_harm);
                    if abs(ratio - round(ratio)) < (t_Bin / pot_PRIs_raw(j_harm))
                        is_harmonic(i) = true;
                        break;
                    end
                end
            end
        end
        PRI_p = pot_PRIs_raw(~is_harmonic); % T?p h?p {PRI_p} chu?n hóa
        pot_idx_filtered = pot_idx(~is_harmonic);
        
        % --- PRI TYPE CLASSIFICATION (JITTER VS. NON-JITTER) ---
        is_jittered = false(1, length(pot_idx_filtered));
        for i = 1:length(pot_idx_filtered)
            if i < length(pot_idx_filtered) && (pot_idx_filtered(i+1) - pot_idx_filtered(i) <= 2)
                is_jittered(i) = true;
                is_jittered(i+1) = true;
            end
        end
        
        jittered_PRIs = PRI_p(is_jittered);
        non_jittered_PRIs = PRI_p(~is_jittered);
        
        fprintf('\n[M?c sai l?ch C=%d] ?ã trích xu?t t?p h?p {PRI_p}:\n', C);
        if ~isempty(non_jittered_PRIs), fprintf('  + Non-Jittered {PRI_p} (ms): %s\n', num2str(non_jittered_PRIs * 1e3, '%.4f  ')); end
        if ~isempty(jittered_PRIs),     fprintf('  + Jittered {PRI_p} (ms):     %s\n', num2str(jittered_PRIs * 1e3, '%.4f  ')); end
        
        success_in_this_C = false;
        max_p_pri = max(PRI_p);
        
        % =================================================================
        % MÔ T? TÌM KI?M CHU?I ?A THAM S? (HÌNH 4 & HÌNH 5)
        % =================================================================
        p = 1; % p: Ch? s? c?a xung b?t ??u th? nghi?m chu?i (p-th input pulse)
        while p <= (length(Xin.TOA) - 3)
            
            % Kh?i t?o tr?ng thái thi?t l?p chu?i (Figure 4)
            k = 1; % k: S? l??ng xung hi?n t?i trong chu?i ?ang tìm ki?m
            temp_seq_idx_local = [p]; 
            matched_pri_indices = [];
            PRI_out = []; % M?ng l?u l?ch s? chu k? th?c t? ph?c v? ph??ng trình (13)
            
            % -------------------------------------------------------------
            % TR?NG THÁI 1: B??C 1 KH?I T?O CHU?I M?I (Figure 4)
            % -------------------------------------------------------------
            while k < 4
                found_next = false;
                last_toa = Xin.TOA(temp_seq_idx_local(end));
                
                best_Z = Inf;
                best_j = -1;
                best_i = -1;
                
                % Tìm ki?m song song trên các ?ng viên j và t?p PRI ti?m n?ng i
                for j = (temp_seq_idx_local(end) + 1):length(Xin.TOA)
                    if (Xin.TOA(j) - last_toa) > 1.5 * max_p_pri
                        break; % Gi?i h?n c?a s? tìm ki?m t?i ?a ?? t?i ?u tính toán
                    end
                    
                    % Ph??ng trình (7): Tính toán kho?ng cách fingerprint v?t lý M_k(j)
                    M_k = ((Xin.RF(j) - Xin.RF(temp_seq_idx_local(end)))^2) / (RF_0^2) + ...
                          ((Xin.PW(j) - Xin.PW(temp_seq_idx_local(end)))^2) / (PW_0^2);
                    
                    if M_k < M0
                        for i = 1:length(PRI_p)
                            current_PRI = PRI_p(i);
                            
                            % C?u hình dung sai th?i gian W d?a trên ??c tính ?i?u ch?
                            if ismember(current_PRI, jittered_PRIs)
                                W = current_PRI * 0.35;
                            else
                                W = 2e-5;
                            end
                            
                            % Ph??ng trình (6): Tính hàm kh?p th?i gian Z_k(i, j)
                            Z_k = abs(current_PRI - (Xin.TOA(j) - last_toa));
                            
                            if Z_k < W
                                if Z_k < best_Z
                                    best_Z = Z_k;
                                    best_j = j;
                                    best_i = i;
                                end
                            end
                        end
                    end
                end
                
                % ?ánh giá k?t qu?: "Tìm th?y xung ti?p theo?" (Find the next pulse?)
                if best_j ~= -1
                    % Tr??ng h?p: CÓ (Yes) -> C?p nh?t ch? s? theo Figure 4
                    k = k + 1;
                    temp_seq_idx_local(end+1) = best_j; 
                    PRI_out(end+1) = Xin.TOA(best_j) - last_toa; 
                    matched_pri_indices(end+1) = best_i; 
                    found_next = true;
                else
                    % Tr??ng h?p: KHÔNG (No) -> ??t chu?i m?i
                    break;
                end
            end
            
            % Ki?m tra ?i?u ki?n k?t thúc B??c 1: ?? 4 xung liên ti?p (success = 1)
            if k == 4
                success = 1;
            else
                success = 0;
            end
            
            % -------------------------------------------------------------
            % TR?NG THÁI 2: B??C 2 BÁM SÁT CHU?I VÀ BÙ M?T XUNG (Figure 5)
            % -------------------------------------------------------------
            if success == 1
                % Khóa ki?u ?i?u ch? c?a ?ài phát d?a trên xung m?i cu?i
                locked_jitter = ismember(PRI_p(matched_pri_indices(end)), jittered_PRIs);
                
                M = 1; % M: H? s? ??m xung m?t liên ti?p (M-multiplier)
                TOA_ref = Xin.TOA(temp_seq_idx_local(end)); % M?c th?i gian tham chi?u
                last_found_local_idx = temp_seq_idx_local(end);
                T_1 = Xin.TOA(end); % T_1: Th?i ?i?m k?t thúc khung d? li?u
                
                while TOA_ref <= T_1
                    best_Z = Inf;
                    best_j = -1;
                    best_i = -1;
                    
                    % Quét tìm xung j th?a mãn trong vùng n?i r?ng ??ng theo M
                    for j = (last_found_local_idx + 1):length(Xin.TOA)
                        if (Xin.TOA(j) - TOA_ref) > 1.5 * M * max_p_pri
                            break;
                        end
                        
                        % Tính kho?ng cách v?t lý M_k(j) t? xung th?c thu g?n nh?t
                        M_k = ((Xin.RF(j) - Xin.RF(last_found_local_idx))^2) / (RF_0^2) + ...
                              ((Xin.PW(j) - Xin.PW(last_found_local_idx))^2) / (PW_0^2);
                        
                        if M_k < M0
                            for i = 1:length(PRI_p)
                                current_PRI = PRI_p(i);
                                is_curr_jitter = ismember(current_PRI, jittered_PRIs);
                                
                                if is_curr_jitter
                                    W_limit = (current_PRI * 0.35) * M;
                                    % Ph??ng trình (12): ??ng b? hóa theo ?? giãn Jittered
                                    Z_k = abs(M * current_PRI - (Xin.TOA(j) - TOA_ref));
                                else
                                    W_limit = 2e-5;
                                    % Ph??ng trình (6): Áp d?ng m?c TOA_ref ?ã d?ch chuy?n ??ng
                                    Z_k = abs(current_PRI - (Xin.TOA(j) - TOA_ref));
                                end
                                
                                if Z_k < W_limit
                                    if Z_k < best_Z
                                        best_Z = Z_k;
                                        best_j = j;
                                        best_i = i;
                                    end
                                end
                            end
                        end
                    end
                    
                    % ?ánh giá k?t qu? bám sát: "Tìm ki?m thành công?" (Successful search?)
                    if best_j ~= -1
                        % Tr??ng h?p: CÓ (Yes) -> Tìm th?y xung th?c t?
                        k = k + 1;
                        temp_seq_idx_local(end+1) = best_j;
                        
                        if ismember(PRI_p(best_i), jittered_PRIs)
                            PRI_out(end+1) = (Xin.TOA(best_j) - TOA_ref) / M;
                        else
                            PRI_out(end+1) = Xin.TOA(best_j) - TOA_ref;
                        end
                        
                        matched_pri_indices(end+1) = best_i;
                        last_found_local_idx = best_j;
                        TOA_ref = Xin.TOA(best_j);
                        M = 1; % Reset h? s? ??m xung m?t v? m?c ??nh
                    else
                        % Tr??ng h?p: KHÔNG (No) -> X? lý m?t xung theo nhánh r? Figure 5
                        if locked_jitter == 1
                            % Nhánh JITTER = 1: Ch? t?ng b? nhân M
                            M = M + 1;
                        else
                            % Nhánh JITTER = 0 (Non-Jittered): T?ng M và d?ch chuy?n m?c TOA_ref
                            M = M + 1;
                            target_index = k - M + 1; % ??nh v? chính xác ph?n t? (k-M) trong MATLAB
                            
                            if target_index >= 1 && target_index <= length(PRI_out)
                                shift_val = PRI_out(target_index);
                                TOA_ref = TOA_ref + shift_val; % T?nh ti?n m?c toán h?c theo ph??ng trình (13)
                            else
                                % Ph??ng phòng v? d? phòng b?ng PRI t?nh g?n nh?t
                                TOA_ref = TOA_ref + PRI_p(matched_pri_indices(end));
                            end
                        end
                        
                        % Ng??ng d?ng bám (H?y bám n?u s? xung m?t liên ti?p v??t quá gi?i h?n)
                        if M > 15 
                            break; 
                        end 
                    end
                end
                
                % -------------------------------------------------------------
                % TR?NG THÁI 3: KI?M TRA ?? DÀI VÀ TRÍCH XU?T CHU?I
                % -------------------------------------------------------------
                if length(temp_seq_idx_local) >= 10
                    % Chuy?n ??i ch? s? n?i b? (local) sang ch? s? lu?ng d? li?u toàn c?c (global)
                    global_indices = active_idx(temp_seq_idx_local);
                    extracted_flags(global_indices) = true; 
                    separated_sequences{end+1} = global_indices; 
                    
                    avg_pri = mean(PRI_p(matched_pri_indices));
                    fprintf('>> ?ã bóc tách thành công 1 ?ài phát! (C=%d, PRI trung bình=%.4f ms, S? xung=%d)\n', ...
                            C, avg_pri * 1e3, length(global_indices));
                    
                    success_in_this_C = true;
                    
                    % Thu h?p lu?ng d? li?u Xin ngay l?p t?c ?? gi?i phóng b? ??m
                    active_idx = find(~extracted_flags);
                    Xin.TOA = TOA_sorted(active_idx);
                    Xin.RF  = RF_sorted(active_idx);
                    Xin.PW  = PW_sorted(active_idx);
                    
                    % ??t l?i p = 1 ?? quét l?i t? ??u trên lu?ng d? li?u s?ch m?i
                    p = 1; 
                    continue; 
                end
            end
            p = p + 1;
        end
        
        % ?i?u khi?n b?c nâng h? C vòng l?p ngoài
        C = iif(success_in_this_C, 1, C + 1);
    end
end

% --- Hàm b? tr? n?i tuy?n ---
function res = iif(cond, true_val, false_val)
    if cond, res = true_val; else, res = false_val; end
end