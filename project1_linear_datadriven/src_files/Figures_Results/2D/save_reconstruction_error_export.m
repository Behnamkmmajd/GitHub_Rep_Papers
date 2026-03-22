function save_reconstruction_error_export(export_data, metadata, export_name)
%SAVE_RECONSTRUCTION_ERROR_EXPORT Save analysis-ready reconstruction-error arrays.

    if nargin < 3 || isempty(export_name)
        export_name = 'reconstruction_error_analysis_2D.mat';
    end

    export_dir = fullfile(fileparts(mfilename('fullpath')), 'exports');
    if ~exist(export_dir, 'dir')
        mkdir(export_dir);
    end

    export_file = fullfile(export_dir, export_name);
    save(export_file, 'export_data', 'metadata', '-v7.3');

    fprintf('Saved reconstruction-error export: %s\n', export_file);
end