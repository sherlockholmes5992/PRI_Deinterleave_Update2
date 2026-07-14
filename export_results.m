function export_results(pdwData, separated_sequences, extracted_flags)
    TOA_sorted = pdwData.TOA_sorted;
    RF_sorted = pdwData.RF_sorted;
    PW_sorted = pdwData.PW_sorted;
    True_ID = pdwData.True_ID;
    N_pulses = pdwData.N_pulses;

    % T?o b?ng t?ng quan d? li?u th?c t? nh?n ???c
    TOA_ms = TOA_sorted(:) * 1e3;    
    RF_GHz = RF_sorted(:) / 1e9;     
    PW_us  = PW_sorted(:) * 1e6;     
    Emitter_ID = True_ID(:);         

    PDW_Table = table(TOA_ms, RF_GHz, PW_us, Emitter_ID, ...
        'VariableNames', {'TOA_ms', 'RF_GHz', 'PW_us', 'Emitter_ID'});

    fprintf('\n>> ?Ã ??NG B? PDW_TABLE TH?C T? (Hi?n th? 20 xung ??u tiên):\n');
    disp(head(PDW_Table, 20));
    assignin('base', 'PDW_Table', PDW_Table);
    try openvar('PDW_Table'); catch; end

    % ?ánh giá ?? chính xác toán h?c d?a trên Ground Truth ban ??u
    fprintf('\n--- ACCURACY EVALUATION (Based on Ground Truth) ---\n');
    total_extracted_pulses = sum(extracted_flags);
    fprintf('Total input pulses: %d\n', N_pulses);
    fprintf('Successfully extracted pulses: %d (Rate: %.2f%%)\n', total_extracted_pulses, (total_extracted_pulses/N_pulses)*100);

    fprintf('\n=======================================================\n');
    fprintf('--- DETAILED TABLE FOR EACH EXTRACTED SEQUENCE ---\n');
    fprintf('=======================================================\n');

    num_extracted = length(separated_sequences);

    for i = 1:num_extracted
        seq_idx = separated_sequences{i};
        
        TOA_seq = TOA_sorted(seq_idx) * 1e3; 
        RF_seq  = RF_sorted(seq_idx) / 1e9;  
        PW_seq  = PW_sorted(seq_idx) * 1e6;  
        ID_seq  = True_ID(seq_idx);          
        
        Extracted_Table = table(TOA_seq(:), RF_seq(:), PW_seq(:), ID_seq(:), ...
            'VariableNames', {'TOA_ms', 'RF_GHz', 'PW_us', 'True_ID'});
        
        % ?ánh giá ?? thu?n khi?t (Purity) t? ??ng
        unique_IDs = unique(ID_seq);
        if length(unique_IDs) == 1
            status = sprintf('Perfect (100%% of pulses belong to original ID %d)', unique_IDs);
        else
            status = sprintf('PULSE MIX-UP! (Contains IDs: %s) -> NEEDS REVIEW', num2str(unique_IDs(:)'));
        end
        
        fprintf('\n>> EXTRACTED SEQUENCE %d (Contains %d pulses)\n', i, length(seq_idx));
        fprintf('   Status: %s\n', status);
        disp(head(Extracted_Table, 10)); 
        
        % ??y bi?n b?ng ra màn hình n?n Workspace ??c l?p ?? ti?n kích ?úp ki?m tra
        table_name = sprintf('Table_Extracted_Seq_%d', i);
        assignin('base', table_name, Extracted_Table);
    end
end