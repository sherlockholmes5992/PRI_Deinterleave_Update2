function export_results(pdwData, separated_sequences, extracted_flags)
    % --- Unpack Pulse Description Word (PDW) parameters from structured data ---
    TOA_sorted = pdwData.TOA_sorted;
    RF_sorted = pdwData.RF_sorted;
    PW_sorted = pdwData.PW_sorted;
    True_ID = pdwData.True_ID;
    N_pulses = pdwData.N_pulses;

    % --- Format and scale received pulse parameters for real-world metric conversion ---
    TOA_ms = TOA_sorted(:) * 1e3;    % Convert Time of Arrival to milliseconds (ms)
    RF_GHz = RF_sorted(:) / 1e9;     % Convert Radio Frequency to Gigahertz (GHz)
    PW_us  = PW_sorted(:) * 1e6;     % Convert Pulse Width to microseconds (us)
    Emitter_ID = True_ID(:);         % Retain original simulated baseline ground truth labels

    % Create global interleaved observation dataset summary table
    PDW_Table = table(TOA_ms, RF_GHz, PW_us, Emitter_ID, ...
        'VariableNames', {'TOA_ms', 'RF_GHz', 'PW_us', 'Emitter_ID'});

    fprintf('\n>> REAL-TIME OBSERVED SIGNAL MATRIX SYNCHRONIZED (Displaying first 20 pulses):\n');
    disp(head(PDW_Table, 20));
    assignin('base', 'PDW_Table', PDW_Table); % Export main dataset to MATLAB workspace environment
    try openvar('PDW_Table'); catch; end      % Launch spreadsheet variable editor UI layout

    % --- Evaluate global algorithm deinterleaving accuracy metrics against ground truth path ---
    fprintf('\n--- ACCURACY EVALUATION (Based on Ground Truth) ---\n');
    total_extracted_pulses = sum(extracted_flags);
    fprintf('Total input pulses: %d\n', N_pulses);
    fprintf('Successfully extracted pulses: %d (Rate: %.2f%%)\n', total_extracted_pulses, (total_extracted_pulses/N_pulses)*100);

    fprintf('\n=======================================================\n');
    fprintf('--- DETAILED DATA Matrix FOR EACH EXTRACTED TRACK ---\n');
    fprintf('=======================================================\n');

    num_extracted = length(separated_sequences);

    % Loop through each successfully isolated radar track cluster block
    for i = 1:num_extracted
        seq_idx = separated_sequences{i};
        
        % Extract parameter track subset data matrices
        TOA_seq = TOA_sorted(seq_idx) * 1e3; 
        RF_seq  = RF_sorted(seq_idx) / 1e9;  
        PW_seq  = PW_sorted(seq_idx) * 1e6;  
        ID_seq  = True_ID(seq_idx);          
        
        Extracted_Table = table(TOA_seq(:), RF_seq(:), PW_seq(:), ID_seq(:), ...
            'VariableNames', {'TOA_ms', 'RF_GHz', 'PW_us', 'True_ID'});
        
        % --- Perform automated cluster purity and leakage metric calculations ---
        unique_IDs = unique(ID_seq);
        if length(unique_IDs) == 1
            status = sprintf('Perfect (100%% of pulses belong to original ID %d)', unique_IDs);
        else
            % Flag a pulse mix-up error if multiple original source tags cross-contaminate the track
            status = sprintf('PULSE MIX-UP FILTER FAULT! (Contains cross-contamination from IDs: %s) -> REVIEW CONFIG', num2str(unique_IDs(:)'));
        end
        
        fprintf('\n>> DEINTERLEAVED SEQUENCE CHANNEL %d (Contains %d recovered pulses)\n', i, length(seq_idx));
        fprintf('   Purity Evaluation: %s\n', status);
        disp(head(Extracted_Table, 10)); % Display top 10 rows for manual terminal check
        
        % Push the sorted standalone variable array to the MATLAB base workspace for inspector view
        table_name = sprintf('Table_Extracted_Seq_%d', i);
        assignin('base', table_name, Extracted_Table);
    end
end