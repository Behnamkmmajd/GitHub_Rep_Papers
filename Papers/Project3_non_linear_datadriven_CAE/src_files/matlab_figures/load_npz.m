function data = load_npz(filepath)
% LOAD_NPZ  Load a Python .npz (NumPy compressed archive) file into MATLAB.
%
%   data = load_npz(filepath)
%
%   .npz files are ZIP archives where each entry is a .npy file.
%   This function extracts them to a temporary directory, reads each array,
%   and returns a struct with one field per array.
%
%   Example:
%       data = load_npz('../data/radial_advection_data.npz');
%       U = data.train_snapshots;   % (15, 16, 256, 256)
%       mu = data.mu_train;         % (15,)
%
%   Supported dtypes: float64, float32, int64, int32, int16, uint8, bool.
%
%   Author: CNN-ROM Tutorial
%   Date:   2026

    % Create temporary directory for extraction
    tmpdir = tempname;
    mkdir(tmpdir);
    cleanup = onCleanup(@() rmdir(tmpdir, 's'));
    
    % Unzip the .npz file (it's a standard ZIP archive)
    unzip(filepath, tmpdir);
    
    % Find all .npy files
    npy_files = dir(fullfile(tmpdir, '*.npy'));
    
    data = struct();
    for k = 1:length(npy_files)
        fname = npy_files(k).name;
        varname = fname(1:end-4);  % Remove .npy extension
        
        % Replace invalid MATLAB field name characters
        varname = strrep(varname, '-', '_');
        varname = strrep(varname, '.', '_');
        
        fullpath = fullfile(tmpdir, fname);
        data.(varname) = read_npy(fullpath);
    end
end


function arr = read_npy(filename)
% READ_NPY  Read a single .npy file.
%
%   The .npy format specification:
%   - 6 bytes magic: \x93NUMPY
%   - 1 byte major version
%   - 1 byte minor version
%   - 2 bytes (v1) or 4 bytes (v2) header length
%   - ASCII header (Python dict with 'descr', 'fortran_order', 'shape')
%   - Raw data

    fid = fopen(filename, 'r');
    if fid == -1
        error('Cannot open file: %s', filename);
    end
    cleanupFile = onCleanup(@() fclose(fid));
    
    % Read magic string (6 bytes)
    magic = fread(fid, 6, 'uint8')';
    assert(magic(1) == 147 && all(magic(2:6) == [78 85 77 80 89]), ...
        'Not a valid .npy file');
    
    % Read version
    major = fread(fid, 1, 'uint8');
    minor = fread(fid, 1, 'uint8');
    
    % Read header length
    if major == 1
        header_len = fread(fid, 1, 'uint16');
    else
        header_len = fread(fid, 1, 'uint32');
    end
    
    % Read header (ASCII Python dict)
    header = char(fread(fid, header_len, 'char')');
    
    % Parse dtype from header
    descr_match = regexp(header, '''descr''\s*:\s*''([^'']+)''', 'tokens');
    dtype_str = descr_match{1}{1};
    
    % Parse shape from header
    shape_match = regexp(header, '''shape''\s*:\s*\(([^)]*)\)', 'tokens');
    shape_str = strtrim(shape_match{1}{1});
    if isempty(shape_str)
        shape = [1, 1];  % Scalar
    else
        shape_str = strrep(shape_str, ',', ' ');
        shape = sscanf(shape_str, '%d')';
        if length(shape) == 1
            shape = [1, shape];  % 1D array -> row vector
        end
    end
    
    % Parse fortran_order (column-major vs row-major)
    fortran_order = contains(header, '''fortran_order'': True');
    
    % Map NumPy dtype to MATLAB type
    % Strip byte order character (<, >, =, |)
    dtype_clean = regexprep(dtype_str, '^[<>=|]', '');
    
    % Check for Unicode string dtype first (U followed by number)
    U_match = regexp(dtype_clean, '^U(\d+)$', 'tokens');
    if ~isempty(U_match)
        char_len = str2double(U_match{1}{1});
        n_elements = prod(shape);
        % Each Unicode char is 4 bytes (UTF-32)
        raw = fread(fid, n_elements * char_len, '*uint32');
        raw = reshape(raw, [char_len, n_elements]);
        arr = cell(1, n_elements);
        for si = 1:n_elements
            s = char(raw(:, si)');
            s(s == 0) = [];  % Remove null padding
            arr{si} = s;
        end
        if n_elements == 1
            arr = arr{1};
        end
        return;
    end

    switch dtype_clean
        case 'f8'
            matlab_type = 'double';
        case 'f4'
            matlab_type = 'single';
        case 'i8'
            matlab_type = 'int64';
        case 'i4'
            matlab_type = 'int32';
        case 'i2'
            matlab_type = 'int16';
        case 'u1'
            matlab_type = 'uint8';
        case 'b1'
            matlab_type = 'uint8';  % Boolean
        otherwise
            error('Unsupported dtype: %s', dtype_str);
    end
    
    % Read the data
    n_elements = prod(shape);
    arr = fread(fid, n_elements, ['*' matlab_type]);
    
    % Reshape
    if ~fortran_order
        % NumPy default is C-order (row-major), MATLAB is column-major
        % Need to reverse shape and permute
        arr = reshape(arr, fliplr(shape));
        arr = permute(arr, length(shape):-1:1);
    else
        arr = reshape(arr, shape);
    end
    
    % Handle scalar
    if isscalar(arr)
        return;
    end
    
    % Squeeze leading singleton for 1D arrays
    if shape(1) == 1 && length(shape) == 2
        arr = arr(:)';  % Return as row vector
    end
end
