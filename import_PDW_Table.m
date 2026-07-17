function [pdwData, initialSDIF] = import_PDW_Table(PDW_Table)
    pdwData.TOA_sorted = PDW_Table.TOA_ms / 1e3;
    pdwData.RF_sorted = PDW_Table.RF_GHz * 1e9;
    pdwData.PW_sorted = PDW_Table.PW_us / 1e6
    pdwData.True_ID = PDW_Table.Emitter_ID;
    pdwData.N_pulses = length(pdwData.TOA_sorted);
    
    % --- Capture initial SDIF Histogram state at difference level C=1 for Plot 1 ---
    diff_TOA_C1 = abs(pdwData.TOA_sorted(1 : end-1) - pdwData.TOA_sorted(2 : end));
    t_Bin_C1 = 1e-5; 
    max_diff_C1 = max(diff_TOA_C1);
    edges_C1 = 0 : t_Bin_C1 : (max_diff_C1 + t_Bin_C1);
    [N_counts_C1, ~] = histcounts(diff_TOA_C1, edges_C1);
    bin_centers_C1 = edges_C1(1:end-1) + t_Bin_C1/2; 

    x_emp_C1 = 0.5; k_emp_C1 = 0.1 / length(bin_centers_C1); 
    tau_C1 = 1:length(bin_centers_C1); 
    Threshold_C1 = x_emp_C1 * (pdwData.N_pulses - 1) * exp(-tau_C1 * k_emp_C1 * length(bin_centers_C1));
    pot_idx_C1 = find(N_counts_C1 > Threshold_C1);

    initialSDIF.bin_centers_C1 = bin_centers_C1;
    initialSDIF.N_counts_C1 = N_counts_C1;
    initialSDIF.Threshold_C1 = Threshold_C1;
    initialSDIF.pot_idx_C1 = pot_idx_C1;
    initialSDIF.max_diff_C1 = max_diff_C1;

end

