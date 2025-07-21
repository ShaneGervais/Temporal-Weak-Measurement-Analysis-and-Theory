%% MATLAB Script for Arrival Time Fitting and Signal Speed Calculation
% This script processes signal arrival times, estimates measurement uncertainty,
% and calculates the signal propagation speed in a BNC cable.

function run_analysis()
    clc; clear; close all;
    warning off;
    % Define base file path
    base_path = "vitesse_dans_les_cables_DET"; % <-- Change this to your data directory
    
    % Automatically detect CSV files and group by distance
    files = dir(fullfile(base_path, '*.csv'));
    file_names = {files.name};
    
    % Extract distances from filenames assuming format:
    % "fixed_delay_4ps_ET_single_DET_distance_cm_number.csv"
    distance_pattern = '_(\d+)_cm';
    distances = cellfun(@(x) regexp(x, distance_pattern, 'tokens', 'once'), file_names, 'UniformOutput', false);
    
    % Ensure all distances are extracted correctly
    valid_indices = ~cellfun('isempty', distances);
    extracted_distances = cellfun(@(x) str2double(x{1}), distances(valid_indices));
    file_names = file_names(valid_indices); % Ensure file_names has matching valid indices
    
    % Unique distances
    unique_distances = unique(extracted_distances);
    fprintf("Detected unique distances: %s\n", mat2str(unique_distances));
    
    % Group files by distance
    distance_groups = struct('distance', {}, 'files', {});
    for i = 1:length(unique_distances)
        dist_value = unique_distances(i);
        matching_files = file_names(extracted_distances == dist_value);
        distance_groups(end+1).distance = dist_value;
        distance_groups(end).files = matching_files;
        
        % Debugging output
        fprintf("Distance: %d cm, Files: %s\n", dist_value, strjoin(matching_files, ', '));
    end
    
    % Choose fits to apply (only polynomials to save computation time)
    selected_fits = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17]; % Polynomial 2, 3, and 4 only
    
    % Run processing function
    process_arrival_time_fitting(base_path, distance_groups, selected_fits);
end

function process_arrival_time_fitting(base_path, distance_groups, selected_fits)
    % Constants
    c = 299792458; % Speed of light in vacuum (m/s)
    speed_theory = 0.66 * c; % Expected speed in BNC cable
    
    % Define available fits
    fit_types = {"poly2", "poly3", "poly4", "poly5", "poly6", "poly7", "poly9", "fourier1", "fourier2", "fourier3", "fourier4", "fourier5", "fourier6", "gauss1", "gauss2", "gauss3", "gauss4"};
    fit_names = {"Polynomial2", "Polynomial3", "Polynomial4", "Polynomial5", "Polynomial6", "Polynomial7", "Polynomial8", "Polynomial9", "Fourier1", "Fourier2", "Fourier3", "Fourier4", "Fourier5", "Fourier6", "Gauss1", "Gauss2", "Gauss3", "Gauss4"};
    
    % Initialize result storage
    arrival_time_table = table();
    velocity_table = table();
    
    % Process each selected fit type
    for f = 1:length(selected_fits)
        fit_type = fit_types{selected_fits(f)};
        fit_name = fit_names{selected_fits(f)};
        
        distances_list = [];
        mean_arrival_times = [];
        std_arrival_times = [];
        fit_quality_scores = [];
        
        for d = 1:length(distance_groups)
            distance = distance_groups(d).distance;
            files = distance_groups(d).files;
            
            if isempty(files)
                warning("No files found for distance %d cm", distance);
                continue;
            end
            
            peak_times = [];
            fit_quality = [];
            
            for i = 1:length(files)
                try
                    file_path = fullfile(base_path, files{i});
                    fprintf("Processing file: %s\n", file_path);
                    [peak_time, quality] = find_arrival_time(file_path, fit_type, distance);
                    if ~isempty(peak_time)
                        peak_times = [peak_times; peak_time];
                        fit_quality = [fit_quality; quality];
                    end
                catch ME
                    warning("Error processing file %s: %s", files{i}, ME.message);
                end
            end
            
            if ~isempty(peak_times)
                distances_list = [distances_list; distance];
                mean_arrival_times = [mean_arrival_times; mean(peak_times)];
                std_arrival_times = [std_arrival_times; std(peak_times)];
                fit_quality_scores = [fit_quality_scores; mean(fit_quality)];
            end
        end
        
        num_rows = length(distances_list);
        fit_names_col = repmat({fit_name}, num_rows, 1);
        arrival_time_table = [arrival_time_table; table(fit_names_col, distances_list, mean_arrival_times, std_arrival_times, fit_quality_scores, 'VariableNames', {'Fit Type', 'Distance (cm)', 'Mean Arrival Time (ps)', 'Std (ps)', 'Goodness of Fit'})];
        
        % Compute velocity fit
        if num_rows > 1
            delays = mean_arrival_times - mean_arrival_times(1);
            p = polyfit(distances_list(2:end), delays(2:end), 1);
            measured_speed = 1 / p(1) * 1e10; % Convert to m/s
            percentage_error = abs((measured_speed - speed_theory) / speed_theory) * 100;
            
            velocity_table = [velocity_table; table({fit_name}, measured_speed, speed_theory, percentage_error, 'VariableNames', {'Fit Type', 'Measured Speed (m/s)', 'Theoretical Speed (m/s)', 'Percentage Error (%)'})];
        end
    end
    
    disp("Arrival Time Table:");
    disp(arrival_time_table);
    writetable(arrival_time_table, fullfile(base_path, 'arrival_time_results.csv'));
    
    disp("Velocity Fit Table:");
    disp(velocity_table);
    writetable(velocity_table, fullfile(base_path, 'velocity_fit_results.csv'));
end

function [peak_time, fit_quality] = find_arrival_time(file_path, fit_type, distance)
    % Load data
    if ~isfile(file_path)
        warning("File not found: %s", file_path);
        peak_time = [];
        fit_quality = [];
        return;
    end
    data = load(file_path);
    if isempty(data)
        warning("File is empty: %s", file_path);
        peak_time = [];
        fit_quality = [];
        return;
    end
    
    amplitude = data(:,1);
    total_duration_ps = 100000;  % total time span in ps
    time = linspace(0, total_duration_ps, length(amplitude));
    
    % 1) Compute and smooth the gradient
    smoothed_grad = smooth(gradient(amplitude, time));
    
    % 2) Find a coarse peak index from the gradient with a high threshold
    %    so we know roughly where the main feature is
    threshold = max(abs(smoothed_grad)) * 0.40; 
    [~, locs] = findpeaks(abs(smoothed_grad), 'MinPeakHeight', threshold);
    if isempty(locs)
        warning("No large gradient peaks found in %s", file_path);
        peak_time = [];
        fit_quality = [];
        return;
    end
    
    % Take the first big gradient peak as a reference
    peak_idx = locs(1);
    coarse_peak_time = time(peak_idx);
    
    % 3) Define ±5 ns region around that peak
    half_window_ps = 5e3;         % ±5 ns = ±5000 ps
    start_time = coarse_peak_time - half_window_ps;
    end_time   = coarse_peak_time + half_window_ps;
    % Clip to data boundaries
    start_idx = find(time >= start_time, 1, 'first');
    end_idx   = find(time <= end_time,  1, 'last');
    
    if isempty(start_idx) || isempty(end_idx)
        warning("Peak window out of bounds for %s", file_path);
        peak_time = [];
        fit_quality = [];
        return;
    end
    
    % Extract the windowed region
    x_win = time(start_idx:end_idx)';
    y_win = abs(smoothed_grad(start_idx:end_idx));
    
    % 4) Filter out the bottom 20% of the amplitude range in that window
    y_max = max(y_win);
    if y_max <= 0
        warning("Window amplitude is zero/negative in %s", file_path);
        peak_time = [];
        fit_quality = [];
        return;
    end
    amp_cutoff = 0.40 * y_max;  % only keep top 80–100%
    keep_mask  = (y_win >= amp_cutoff);
    x_filt = x_win(keep_mask);
    y_filt = y_win(keep_mask);
    
    if length(x_filt) < 5
        % Not enough points to fit
        warning("Too few points left after filtering in %s", file_path);
        peak_time = [];
        fit_quality = [];
        return;
    end
    
    % 5) Fit the polynomial (poly2, poly3, or poly4) to the filtered data
    fit_obj = fit(x_filt, y_filt, fit_type);
    
    % 6) Make a higher‐resolution time axis *just* over the filtered domain
    x_hr = linspace(min(x_filt), max(x_filt), 2000);
    y_hr = feval(fit_obj, x_hr);
    
    % 7) Find arrival by derivative = 0 inside that domain
    %    (Because we’re using polynomials, we can just differentiate fit_obj.)
    d_fit = differentiate(fit_obj, x_hr);
    [~, zero_idx] = min(abs(d_fit));   % find the derivative crossing zero
    peak_time = x_hr(zero_idx);
    
    % 8) Evaluate fit quality (NRMSE) on *filtered* region
    fit_quality = goodnessOfFit(feval(fit_obj, x_filt), y_filt, 'NRMSE');
    
    % % 9) Plot for debugging
    % figure('Name', sprintf('Distance %d cm using %s', distance, fit_type), ...
    %        'NumberTitle', 'off');
    % hold on;
    % % Original window in blue
    % plot(x_win, y_win, 'bo', 'MarkerSize', 4, 'DisplayName', 'Original Window');
    % % Filtered points in black
    % plot(x_filt, y_filt, 'ko', 'MarkerSize', 5, 'DisplayName', 'Filtered (80%-100%)');
    % % Fitted curve in red
    % plot(x_hr, y_hr, 'r-', 'LineWidth', 1.3, 'DisplayName', 'Fit');
    % % Mark detected “peak time” in magenta
    % plot(peak_time, feval(fit_obj, peak_time), 'mo', 'MarkerFaceColor','m', ...
    %      'DisplayName', 'Detected Peak');
    % 
    % title(sprintf('Fit for Distance %d cm using %s', distance, fit_type), 'Interpreter','none');
    % xlabel('Time (ps)');
    % ylabel('Amplitude');
    % legend('Location','best');
    % grid on;
    % hold off;
    
    fprintf("Detected peak at time: %.2f ps for distance: %d cm\n", peak_time, distance);
end
