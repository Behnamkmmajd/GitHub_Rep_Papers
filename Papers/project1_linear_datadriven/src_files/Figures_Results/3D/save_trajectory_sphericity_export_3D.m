function save_trajectory_sphericity_export_3D(Re, Bo, time_nd, export_data, metadata)
%SAVE_TRAJECTORY_SPHERICITY_EXPORT_3D Save 3D trajectory and sphericity arrays and RMSD values.

    export_dir = fullfile(fileparts(mfilename('fullpath')), 'exports');
    if ~exist(export_dir, 'dir')
        mkdir(export_dir);
    end

    rmsd_percent = struct();
    rmsd_percent.sphericity = struct( ...
        'm_pod', compute_rmsd_percent(export_data.sphericity.dns, export_data.sphericity.m_pod), ...
        'm_dmd', compute_rmsd_percent(export_data.sphericity.dns, export_data.sphericity.m_dmd));

    dns_traj = export_data.trajectory.dns;
    pod_traj = export_data.trajectory.m_pod;
    dmd_traj = export_data.trajectory.m_dmd;
    dns_traj_mag = sqrt(sum(dns_traj .^ 2, 2));
    pod_err_mag = sqrt(sum((pod_traj - dns_traj) .^ 2, 2));
    dmd_err_mag = sqrt(sum((dmd_traj - dns_traj) .^ 2, 2));
    dns_traj_rms = sqrt(mean(dns_traj_mag .^ 2));
    if dns_traj_rms == 0
        rmsd_percent.trajectory = struct('m_pod', NaN, 'm_dmd', NaN);
    else
        rmsd_percent.trajectory = struct( ...
            'm_pod', 100 * sqrt(mean(pod_err_mag .^ 2)) / dns_traj_rms, ...
            'm_dmd', 100 * sqrt(mean(dmd_err_mag .^ 2)) / dns_traj_rms);
    end

    rmsd_definition = 'Scalar RMSD_percent = 100 * sqrt(mean((q_rec - q_dns).^2)) / sqrt(mean(q_dns.^2)). Trajectory RMSD uses the Euclidean centroid-position error normalized by the RMS centroid-position magnitude of DNS.';
    export_file = fullfile(export_dir, sprintf('trajectory_sphericity_analysis_3D_Re%d_Bo%d.mat', Re, Bo));
    save(export_file, 'Re', 'Bo', 'time_nd', 'export_data', 'metadata', 'rmsd_percent', 'rmsd_definition', '-v7.3');

    fprintf('Saved 3D trajectory/sphericity export: %s\n', export_file);
end