function save_rise_velocity_export_3D(Re, Bo, time_nd, export_data, metadata)
%SAVE_RISE_VELOCITY_EXPORT_3D Save analysis-ready 3D rise-velocity arrays and RMSD values.

    export_dir = fullfile(fileparts(mfilename('fullpath')), 'exports');
    if ~exist(export_dir, 'dir')
        mkdir(export_dir);
    end

    variant_names = fieldnames(export_data.centroid);
    variant_names = setdiff(variant_names, {'dns'}, 'stable');

    rmsd_percent = struct('centroid', struct(), 'mean', struct());
    dns_centroid = export_data.centroid.dns;
    dns_mean = export_data.mean.dns;

    for idx = 1:numel(variant_names)
        variant_name = variant_names{idx};
        rmsd_percent.centroid.(variant_name) = compute_rmsd_percent(dns_centroid, export_data.centroid.(variant_name));
        rmsd_percent.mean.(variant_name) = compute_rmsd_percent(dns_mean, export_data.mean.(variant_name));
    end

    rmsd_definition = 'RMSD_percent = 100 * sqrt(mean((q_rec - q_dns).^2)) / sqrt(mean(q_dns.^2)).';
    export_file = fullfile(export_dir, sprintf('rise_velocity_analysis_3D_Re%d_Bo%d.mat', Re, Bo));
    save(export_file, 'Re', 'Bo', 'time_nd', 'export_data', 'metadata', 'rmsd_percent', 'rmsd_definition', '-v7.3');

    fprintf('Saved 3D rise-velocity export: %s\n', export_file);
end