function save_interface_export(Re, Bo, snapshot_data, metadata)
%SAVE_INTERFACE_EXPORT Save analysis-ready interface arrays for manuscript review.

    export_dir = fullfile(fileparts(mfilename('fullpath')), 'exports');
    if ~exist(export_dir, 'dir')
        mkdir(export_dir);
    end

    export_file = fullfile(export_dir, sprintf('interface_analysis_Re%d_Bo%d.mat', Re, Bo));
    save(export_file, 'Re', 'Bo', 'snapshot_data', 'metadata', '-v7.3');

    fprintf('Saved interface export: %s\n', export_file);
end