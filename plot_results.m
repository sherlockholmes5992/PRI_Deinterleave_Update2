function plot_results(pdwData, separated_sequences, initialSDIF, num_emitters)
    fprintf('\n--- STARTING VISUAL VERIFICATION PLOTS ---\n');
    
    % --- Unpack Pulse Description Word (PDW) parameters ---
    TOA_sorted = pdwData.TOA_sorted;
    RF_sorted = pdwData.RF_sorted;
    PW_sorted = pdwData.PW_sorted;
    True_ID = pdwData.True_ID;
    N_pulses = pdwData.N_pulses;

    % --- Unpack initial histogram parameters (at first difference rank C=1) ---
    bin_centers_C1 = initialSDIF.bin_centers_C1;
    N_counts_C1 = initialSDIF.N_counts_C1;
    Threshold_C1 = initialSDIF.Threshold_C1;
    pot_idx_C1 = initialSDIF.pot_idx_C1;
    max_diff_C1 = initialSDIF.max_diff_C1;
    num_extracted = length(separated_sequences);

    % -------------------------------------------------------------------------
    % PLOT 1: SDIF Histogram (C=1) - Verify candidate PRI spectral lines exceeding threshold
    % -------------------------------------------------------------------------
    figure('Name', 'Verification 1: SDIF Histogram at C=1', 'Position', [100, 100, 800, 450]);
    stem(bin_centers_C1 * 1e3, N_counts_C1, 'Marker', 'none', 'LineWidth', 1.5, 'Color', [0.2 0.4 0.8]);
    hold on; grid on;
    plot(bin_centers_C1 * 1e3, Threshold_C1, 'r--', 'LineWidth', 2, 'DisplayName', 'Optimal Threshold');
    if ~isempty(pot_idx_C1)
        plot(bin_centers_C1(pot_idx_C1) * 1e3, N_counts_C1(pot_idx_C1), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 5, 'DisplayName', 'Potential PRI');
    end
    title('Sequential Difference Histogram (SDIF) at C = 1');
    xlabel('TOA Difference (ms)'); ylabel('Count (Bin)');
    legend('Location', 'northeast'); xlim([0, min(max_diff_C1*1e3, 3.5)]);

    % -------------------------------------------------------------------------
    % PLOT 2: PDW Scatter (RF vs PW) - Emitter fingerprint cluster map
    % -------------------------------------------------------------------------
    figure('Name', 'Verification 2: PDW Space (RF vs PW)', 'Position', [150, 150, 700, 450]);
    colors = lines(num_emitters); hold on; grid on;
    for emitter_idx = 1:num_emitters 
        idx = (True_ID == emitter_idx);
        scatter(RF_sorted(idx) / 1e9, PW_sorted(idx) * 1e6, 36, colors(emitter_idx, :), 'filled', 'DisplayName', sprintf('Original Radar %d', emitter_idx));
    end
    title('PDW Parameter Scatter Plot (Emitter Fingerprints)');
    xlabel('Carrier Frequency - RF (GHz)'); ylabel('Pulse Width - PW (\mu s)');
    legend('Location', 'best');

    % -------------------------------------------------------------------------
    % PLOT 3: Pulse Timeline - Extraction progress of overlapping pulse streams
    % -------------------------------------------------------------------------
    figure('Name', 'Verification 3: Pulse Timeline', 'Position', [200, 200, 1000, 400]);
    hold on; grid on;
    stem(TOA_sorted * 1e3, zeros(1, N_pulses), 'k', 'Marker', 'none', 'LineWidth', 1);
    c_map = lines(num_extracted);
    for i = 1:num_extracted
        seq_idx = separated_sequences{i};
        extracted_TOA = TOA_sorted(seq_idx);
        stem(extracted_TOA * 1e3, ones(1, length(seq_idx)) * i, 'Color', c_map(i,:), 'Marker', '^', 'MarkerFaceColor', c_map(i,:), 'MarkerSize', 5, 'LineWidth', 1);
    end
    ylim([-1, num_extracted + 1]); yticks(0:num_extracted);
    ytick_labels = {'Interleaved Stream'};
    for i = 1:num_extracted, ytick_labels{end+1} = sprintf('Seq %d', i); end
    yticklabels(ytick_labels);
    xlabel('Time of Arrival - TOA (ms)'); title('Visualization of Interleaved Pulses and Extraction Results');

    % -------------------------------------------------------------------------
    % PLOT 4: PRI Variation - Fixed/Staggered/Jittered modulation behavior verification profiles
    % -------------------------------------------------------------------------
    if num_extracted > 0
        figure('Name', 'Verification 4: PRI Variation', 'Position', [250, 250, 1200, 600]);
        plots_to_draw = min(num_extracted, 6); 
        for i = 1:plots_to_draw
            seq_idx = separated_sequences{i}; 
            extracted_PRI_variations = diff(TOA_sorted(seq_idx)); 
            subplot(2, 3, i);
            plot(1:length(extracted_PRI_variations), extracted_PRI_variations * 1e3, '-b', 'LineWidth', 1.2);
            grid on; ylim([0, 3]); 
            title(sprintf('Extracted Seq %d (Pulses: %d)', i, length(seq_idx)));
            xlabel('Pulse Index'); ylabel('PRI (ms)');
        end
        
        % Dynamic title rendering block handling pre-R2018b fallback paths
        try
            sgtitle('PRI Variation of Extracted Sequences (Simulation of Fig 6)', 'FontSize', 14, 'FontWeight', 'bold');
        catch
            axes('Position', [0 0 1 1], 'Visible', 'off');
            text(0.5, 0.96, 'PRI Variation of Extracted Sequences (Simulation of Fig 6)', 'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold');
        end
    end
end