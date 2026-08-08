function Target = PreprocessTargetData(Target, Include_min_COV, Resample, Number, Stat_3rd_Col)
% Imports and preprocesses the target dispersion curve: reads the file, converts
% slowness to velocity, builds the std/COV column, and optionally log-resamples.
% Author      : Mrinal Bhaumik
% Affiliation : Utah State University, 2026

%% Import data files ------------------------------------------------------------------------------
if ~isempty(Target.fv_file)
    fv_file        = Target.fv_file;
    Data_fv        = readmatrix(fv_file);
    Target.fv_file = Data_fv;
    Target.name    = fv_file;
end

%% Preprocess fv (dispersion) data ----------------------------------------------------------------
% Convert slowness to velocity if present.
isSlowness = any(Target.fv_file(:, 2) < 1);
if isSlowness
    Target.fv_file(:, 2) = 1 ./ Target.fv_file(:, 2);
end

[~, numCols] = size(Target.fv_file);
if numCols >= 3
    % A 3rd column exists: convert it to the standard deviation of PHASE VELOCITY.
    % Its meaning is set by Stat_3rd_Col ('COV', else treated as 'STD') and by whether
    % col-2 was slowness. For v = 1/s, error propagation gives  std_v = v^2 * std_s.
    v     = Target.fv_file(:, 2);       % phase velocity (col-2 already converted above)
    c3    = Target.fv_file(:, 3);       % raw 3rd column (COV or STD, of velocity or slowness)
    isCOV = strcmpi(Stat_3rd_Col, 'COV');

    if isSlowness
        if isCOV
            % COV of slowness -> std of velocity:  std_v = v * COV_s  (linear).
            % Large COVs (>1) get the nonlinear reciprocal correction, per point.
            COV_s        = c3;
            big          = COV_s > 1;
            COV_s(big)   = COV_s(big) - sqrt(COV_s(big).^2 - 2*COV_s(big) + 2);
            stdv         = v .* COV_s;
        else
            % STD of slowness -> std of velocity:  std_v = v^2 * std_s.
            stdv = (v.^2) .* c3;
        end
    else
        if isCOV
            stdv = v .* c3;             % COV of velocity -> std_v = v * COV
        else
            stdv = c3;                  % already the std of velocity -> use as-is
        end
    end
    Target.fv_file(:, 3) = stdv;

elseif ~isempty(Include_min_COV)
    % No 3rd column: synthesise a std column from a default COV of velocity.
    covDefault     = Include_min_COV(1);
    stdCol         = Target.fv_file(:, 2) * covDefault;
    Target.fv_file = [Target.fv_file, stdCol];
end

%% Optional log-scale resampling ------------------------------------------------------------------
% Resamples the data on a log-frequency scale.
if ~strcmpi(Resample, 'yes')
    return;                 % No resampling needed
end
if isempty(Number)
    disp('Resample "Number" is missing, considering default 40 samples')
    Number = 40;
end
Target.fv_file = Get_Resample(Target.fv_file);

%% Nested: log-frequency resampling ---------------------------------------------------------------
    function Data = Get_Resample(Data)
        hasStd   = size(Data, 2) >= 3;      % is a 3rd (std) column present?

        freq     = Data(:, 1);      % original frequency
        velocity = Data(:, 2);
        min_freq = min(freq);       % or user defined value
        max_freq = max(freq);       % or user defined value

        % Ensure data is sorted
        [freq, sortIdx] = sort(freq);
        velocity = velocity(sortIdx);

        logFreq = logspace(log10(min_freq), log10(max_freq), Number)';

        % Interpolate using shape-preserving piecewise cubic interpolation
        velocityInterp = interp1(freq, velocity, logFreq, 'pchip', 'extrap');

        % Return resampled data (carry the std column through only if it exists)
        if hasStd
            stdDev    = Data(sortIdx, 3);
            stdInterp = interp1(freq, stdDev, logFreq, 'pchip', 'extrap');
            Data      = [logFreq, velocityInterp, stdInterp];
        else
            Data      = [logFreq, velocityInterp];
        end
    end

end
