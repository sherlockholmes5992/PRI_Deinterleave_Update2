function [separated_sequences, extracted_flags] = radar_core(pdwData, algoParams)
%#codegen

    % G?i hàm x? lý lõi SDIF
    [separated_sequences, extracted_flags] = sdif_deinterleave(pdwData, algoParams);

end