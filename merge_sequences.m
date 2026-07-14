function merged_sequences = merge_sequences(pdwData, separated_sequences, mergeParams)
    fprintf('\n--- STARTING ADVANCED HARMONIC SEQUENCE MERGING ---\n');
    
    % --- Extract packed structural parameter vectors ---
    TOA_sorted = pdwData.TOA_sorted;
    RF_sorted = pdwData.RF_sorted;
    PW_sorted = pdwData.PW_sorted;

    % --- Load post-processing boundary thresholds ---
    rf_merge_tol = mergeParams.rf_merge_tol;
    pw_merge_tol = mergeParams.pw_merge_tol;
    pri_ratio_tol = mergeParams.pri_ratio_tol;

    % --- Initialize sequence track arrays ---
    merged_sequences = separated_sequences;
    has_changed = true;

    % --- Iterative scan loop until no further sequences satisfy the merge criteria ---
    while has_changed
        has_changed = false;
        num_seq = length(merged_sequences);
        
        i = 1;
        while i <= num_seq
            j = i + 1;
            while j <= num_seq
                idx_i = merged_sequences{i};
                idx_j = merged_sequences{j};
                
                % --- Compute the mean physical fingerprints of both sequences for feature matching ---
                mean_rf_i = mean(RF_sorted(idx_i)); mean_pw_i = mean(PW_sorted(idx_i));
                mean_rf_j = mean(RF_sorted(idx_j)); mean_pw_j = mean(PW_sorted(idx_j));
                
                % Estimate base baseline PRI using the median of sequential differences
                if length(idx_i) >= 2, pri_i = median(diff(TOA_sorted(idx_i))); else, pri_i = 1e-3; end
                if length(idx_j) >= 2, pri_j = median(diff(TOA_sorted(idx_j))); else, pri_j = 1e-3; end
                
                % --- CONDITION 1: Match intrinsic physical characteristics (RF & PW consistency) ---
                if abs(mean_rf_i - mean_rf_j) < rf_merge_tol && ...
                   abs(mean_pw_i - mean_pw_j) < pw_merge_tol
                    
                    ratio = pri_i / pri_j;
                    is_harmonic_relation = false;
                    
                    % --- CONDITION 2: Evaluate mathematical harmonic ratios or grid layout failures due to dropouts ---
                    % Case A: PRIs match closely (within 10% tolerance due to jitter or minor pulse dropouts)
                    if abs(pri_i - pri_j) / min(pri_i, pri_j) < 0.10
                        is_harmonic_relation = true;
                    % Case B: Sequence i is an integer multiple of sequence j (True Harmonic relationship)
                    elseif ratio > 1.5 && abs(ratio - round(ratio)) < pri_ratio_tol
                        is_harmonic_relation = true;
                    % Case C: Sequence j is an integer multiple of sequence i (Inverse Harmonic relationship)
                    elseif ratio < 0.7 && abs((1/ratio) - round(1/ratio)) < pri_ratio_tol
                        is_harmonic_relation = true;
                    end
                    
                    if is_harmonic_relation
                        % --- GUARD LIMIT: Protect independent parallel co-channel emitters ---
                        % Tentatively blend the TOA sets of both tracks and sort chronologically
                        combined_toa = sort(TOA_sorted(union(idx_i, idx_j)));
                        min_pri_base = min(pri_i, pri_j);
                        
                        % If pulses become excessively dense (<40% base PRI), they belong to separate
                        % co-channel systems operating concurrently. Reject false merge.
                        if any(diff(combined_toa) < 0.4 * min_pri_base)
                            j = j + 1; 
                            continue; % Abort merge path for parallel co-channel emitters
                        end
                        
                        % --- SUCCESSFUL TRACK CONCATENATION ---
                        fprintf('>> [MERGE RECOVERED] Successfully merged harmonic fragments: Sequence %d and %d (RF ~ %.2f GHz, PRI: %.2f ms / %.2f ms)\n', ...
                            i, j, mean_rf_i/1e9, pri_i*1e3, pri_j*1e3);
                        
                        % Merge indices and sort array entries
                        merged_sequences{i} = sort(union(idx_i, idx_j));
                        merged_sequences(j) = []; % Purge sequence j now absorbed by i
                        
                        has_changed = true;
                        num_seq = num_seq - 1; % Readjust counter bounds
                        break; % Escape loop level to restart scan matrix with updated topology
                    else
                        j = j + 1;
                    end
                else
                    j = j + 1;
                end
            end
            if has_changed, break; end
            i = i + 1;
        end
    end
    fprintf('--- END OF ADVANCED POST-PROCESSING: Resolved into %d complete emitter sources ---\n\n', length(merged_sequences));
end