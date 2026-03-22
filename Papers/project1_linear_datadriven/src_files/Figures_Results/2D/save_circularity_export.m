function save_circularity_export(Re, Bo, time_nd, circ, metadata)
%SAVE_CIRCULARITY_EXPORT Save analysis-ready circularity arrays for manuscript review.

    dns_circularity = circ(:, 3);
    rmsd_percent = struct( ...
        'e_pod', compute_rmsd_percent(dns_circularity, circ(:, 1)), ...
        'm_pod', compute_rmsd_percent(dns_circularity, circ(:, 2)), ...
        'e_dmd', compute_rmsd_percent(dns_circularity, circ(:, 4)), ...
        'm_dmd', compute_rmsd_percent(dns_circularity, circ(:, 5)), ...
        'l_pod', compute_rmsd_percent(dns_circularity, circ(:, 6)), ...
        'l_dmd', compute_rmsd_percent(dns_circularity, circ(:, 7)));

    rmsd_definition = 'RMSD_percent = 100 * sqrt(mean((q_rec - q_dns).^2)) / sqrt(mean(q_dns.^2)).';

    export_dir = fullfile(fileparts(mfilename('fullpath')), 'exports');
    if ~exist(export_dir, 'dir')
        mkdir(export_dir);
    end

    export_file = fullfile(export_dir, sprintf('circularity_analysis_Re%d_Bo%d.mat', Re, Bo));
    save(export_file, 'Re', 'Bo', 'time_nd', 'circ', 'metadata', 'rmsd_percent', 'rmsd_definition', '-v7.3');

    fprintf('Saved circularity export: %s\n', export_file);
end