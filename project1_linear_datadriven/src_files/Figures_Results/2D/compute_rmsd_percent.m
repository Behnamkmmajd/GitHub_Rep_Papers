function rmsd_percent = compute_rmsd_percent(reference_series, reconstructed_series)
%COMPUTE_RMSD_PERCENT Compute normalized RMSD (%) relative to a reference series.

    reference_series = reference_series(:);
    reconstructed_series = reconstructed_series(:);

    valid = isfinite(reference_series) & isfinite(reconstructed_series);
    if ~any(valid)
        rmsd_percent = NaN;
        return;
    end

    reference_series = reference_series(valid);
    reconstructed_series = reconstructed_series(valid);

    reference_rms = sqrt(mean(reference_series .^ 2));
    if reference_rms == 0
        rmsd_percent = NaN;
        return;
    end

    rmsd_percent = 100 * sqrt(mean((reconstructed_series - reference_series) .^ 2)) / reference_rms;
end