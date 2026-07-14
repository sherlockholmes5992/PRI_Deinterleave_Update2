function merged_sequences = merge_sequences(pdwData, separated_sequences, mergeParams)
    fprintf('\n--- STARTING ADVANCED HARMONIC SEQUENCE MERGING ---\n');
    
    TOA_sorted = pdwData.TOA_sorted;
    RF_sorted = pdwData.RF_sorted;
    PW_sorted = pdwData.PW_sorted;

    rf_merge_tol = mergeParams.rf_merge_tol;
    pw_merge_tol = mergeParams.pw_merge_tol;
    pri_ratio_tol = mergeParams.pri_ratio_tol;

    merged_sequences = separated_sequences;
    has_changed = true;

    while has_changed
        has_changed = false;
        num_seq = length(merged_sequences);
        
        i = 1;
        while i <= num_seq
            j = i + 1;
            while j <= num_seq
                idx_i = merged_sequences{i};
                idx_j = merged_sequences{j};
                
                % Tính toán ??c tr?ng trung bình c?a 2 chu?i ?? ??i chi?u vân tay
                mean_rf_i = mean(RF_sorted(idx_i)); mean_pw_i = mean(PW_sorted(idx_i));
                mean_rf_j = mean(RF_sorted(idx_j)); mean_pw_j = mean(PW_sorted(idx_j));
                
                if length(idx_i) >= 2, pri_i = median(diff(TOA_sorted(idx_i))); else, pri_i = 1e-3; end
                if length(idx_j) >= 2, pri_j = median(diff(TOA_sorted(idx_j))); else, pri_j = 1e-3; end
                
                % ?i?u ki?n 1: Kh?p hoàn toàn thông s? v?t lý (RF & PW)
                if abs(mean_rf_i - mean_rf_j) < rf_merge_tol && ...
                   abs(mean_pw_i - mean_pw_j) < pw_merge_tol
                    
                    ratio = pri_i / pri_j;
                    is_harmonic_relation = false;
                    
                    % Ki?m tra m?i quan h? toán h?c sóng hài ho?c ??t gãy l??i do m?t xung
                    if abs(pri_i - pri_j) / min(pri_i, pri_j) < 0.10
                        is_harmonic_relation = true;
                    elseif ratio > 1.5 && abs(ratio - round(ratio)) < pri_ratio_tol
                        is_harmonic_relation = true;
                    elseif ratio < 0.7 && abs((1/ratio) - round(1/ratio)) < pri_ratio_tol
                        is_harmonic_relation = true;
                    end
                    
                    if is_harmonic_relation
                        % Ch?t ch?n b?o v? ?ài ??ng kênh song song
                        combined_toa = sort(TOA_sorted(union(idx_i, idx_j)));
                        min_pri_base = min(pri_i, pri_j);
                        
                        if any(diff(combined_toa) < 0.4 * min_pri_base)
                            j = j + 1; 
                            continue; % T? ch?i g?p n?u xung ??t l??i th?i gian th?c t?
                        end
                        
                        fprintf('>> [MERGE RECOVERED] G?p thành công m?nh sóng hài: Chu?i %d và %d (RF ~ %.2f GHz, PRI: %.2f ms / %.2f ms)\n', ...
                            i, j, mean_rf_i/1e9, pri_i*1e3, pri_j*1e3);
                        
                        merged_sequences{i} = sort(union(idx_i, idx_j));
                        merged_sequences(j) = []; 
                        
                        has_changed = true;
                        num_seq = num_seq - 1;
                        break;
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
    fprintf('--- K?T THÚC H?U X? LÝ CAO C?P: Phân tách thành %d ngu?n phát hoàn ch?nh ---\n\n', length(merged_sequences));
end