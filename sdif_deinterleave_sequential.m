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
    extracted_flags = false(1, N_pulses); % Boolean mask for processed pulses
    separated_sequences = {};             % Cell array to store extracted pulse indices
    C = 1;                                % Initialize to the first difference level

    while C <= C_max
        % Isolate remaining active (unextracted) pulses
        active_idx = find(~extracted_flags);
        if length(active_idx) < 10 % 
            break; % Terminate if remaining pulse count is insufficient
        end 
        
        % Extract active parameter vectors using index subsets
        active_TOA = active_idx_subset(TOA_sorted, active_idx);
        active_RF  = active_idx_subset(RF_sorted, active_idx);
        active_PW  = active_idx_subset(PW_sorted, active_idx);
        
        if length(active_TOA) <= C
            break; % Terminate if active vector length is smaller than difference level
        end
        
        % --- COMPUTE SDIF HISTOGRAM AT DIFFERENCE LEVEL C ---
        diff_TOA = abs(active_TOA(1 : end-C) - active_TOA(1+C : end));
        diff_TOA = round(diff_TOA, 7); % Round to 0.1 microseconds to suppress floating-point jitter
        max_diff = max(diff_TOA);
        edges = 0 : t_Bin : (max_diff + t_Bin);
        [N_counts, ~] = histcounts(diff_TOA, edges);
        bin_centers = edges(1:end-1) + t_Bin/2; 
        
        % --- CALCULATE DYNAMIC OPTIMAL THRESHOLD CURVE ---
        E_active = length(active_TOA);
        N_bins = length(bin_centers);
        tau = 1:N_bins; 
        Threshold = x_emp * (E_active - C) * exp(-tau ./ (k_emp * N_bins));
        
        % --- EXTRACT POTENTIAL RADAR PRI CANDIDATES ---
        pot_idx = find(N_counts > Threshold);
        pot_PRIs_raw = bin_centers(pot_idx);
        
        if isempty(pot_PRIs_raw)
            C = C + 1; % Increment difference rank if no peaks breach the threshold
            continue;
        end
        
        % --- HARMONIC SUPPRESSION FILTER ---
        is_harmonic = false(1, length(pot_PRIs_raw));
        for i = 1:length(pot_PRIs_raw)
            for j = 1:(i-1)
                if ~is_harmonic(j)
                    ratio = pot_PRIs_raw(i) / pot_PRIs_raw(j);
                    % Check if the candidate is an integer multiple within bin tolerance
                    if abs(ratio - round(ratio)) < (t_Bin / pot_PRIs_raw(j))
                        is_harmonic(i) = true;
                        break;
                    end
                end
            end
        end
        pot_PRIs = pot_PRIs_raw(~is_harmonic);
        pot_idx_filtered = pot_idx(~is_harmonic);
        
        % --- PRI TYPE CLASSIFICATION (JITTER VS. NON-JITTER) ---
        % Classify based on the adjacency/density of triggered bins
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
        
        % Print real-time diagnostics to the Command Window for debugging
        fprintf('\n[Difference Level C=%d] Potential PRIs (PRI_p) exceeding threshold:\n', C);
        if ~isempty(non_jittered_PRIs), fprintf('  + Non-Jittered PRIs (ms): %s\n', num2str(non_jittered_PRIs * 1e3, '%.4f  ')); end
        if ~isempty(jittered_PRIs),     fprintf('  + Jittered PRIs (ms):     %s\n', num2str(jittered_PRIs * 1e3, '%.4f  ')); end
        
        success_in_this_C = false;
        
        % --- MULTI-PARAMETER SEQUENCE SEARCH ENGINE ---
        for i = 1:length(all_pot_PRIs)
            current_PRI = all_pot_PRIs(i);
            
            % Dynamically scale the time tolerance window (W) based on modulation class
            W = iif(ismember(current_PRI, jittered_PRIs), current_PRI * 0.35, 2e-5);
            
            p = 1; 
            while p <= length(active_TOA)
                k = p;
                temp_seq_idx_local = k; 
                pulses_found = 1;
                
                % PHASE 1: Establish Anchor Sequence (First 4 pulses)
                while pulses_found < 4
                    found_next = false;
                    for j = (k+1):length(active_TOA)
                        TOA_diff = active_TOA(j) - active_TOA(k);
                        Z_k = abs(current_PRI - TOA_diff);
                        
                        % Compute multi-parameter normalized fingerprint distance metric (M_k)
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
                
                % PHASE 2: Track Pulse Stream with Missing Pulse Compensation
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
                                M_multiplier = 1; % Reset multiplier upon pulse recovery
                                found_next = true;
                                break;
                            end
                        end
                        
                        if ~found_next
                            M_multiplier = M_multiplier + 1;
                            if M_multiplier > 15 
                                break; % Break sequence if missing consecutive dropouts exceed 15
                            end 
                        end
                    end
                    
                    % PHASE 3: Minimum Sequence Length Acceptance Test
                    if length(temp_seq_idx_local) >= 10
                        % Map local active indices back to global stream indices
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
        
        % Hasani (2021) Rank control update logic loop
        C = iif(success_in_this_C, 1, C + 1); % Reset C to 1 on success to capture residual bins; increment on failure
    end
end

% --- Inline functional helpers for clean code structure ---
function val = active_idx_subset(arr, idx), val = arr(idx); end
function res = iif(cond, true_val, false_val), if cond, res = true_val; else, res = false_val; end; end