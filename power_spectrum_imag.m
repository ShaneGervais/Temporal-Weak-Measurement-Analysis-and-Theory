clear; clc; close all;

% ——— Parameters ———
polarization_states = {'h','v','r','a','l','d'};
data_folder        = './';
threshold          = 0.01;            % Min PSD to consider as a peak
searchHalfWidth    = 20;           % Points around each peak for local fit
f_center           = 703.12e6;     % Center frequency (Hz)
f_span             = 1.4062e9;     % Total frequency span (Hz)
dt                 = 4e-12;        % 4 ps sampling interval
fsamp              = 1/dt;         % Sample rate for pspectrum
tolerance_Hz = 50e6;

% Prepare figure
figure; hold on;
colors = lines(numel(polarization_states));

% Storage for each polarization’s fitted peaks (in Hz)
peakMap = containers.Map(polarization_states, repmat({[]},1,numel(polarization_states)));

% ——— Loop over each polarization file ———
for idx = 1:numel(polarization_states)
    pol  = polarization_states{idx};
    file = fullfile(data_folder, [pol '_int_fft.csv']);
    if ~exist(file,'file')
        warning("Missing %s", file);
        continue;
    end

    raw_y = load(file);
    raw_y = raw_y(:);

    % 1) Compute PSD
    [y, ~] = pspectrum(raw_y, fsamp);

    % 2) Build a correct frequency axis in Hz
    freq_Hz = linspace(...
        f_center - f_span/2, ...
        f_center + f_span/2, ...
        numel(y) ...
    );

    % 3) Plot raw PSD
    scatter(freq_Hz, y, 8, colors(idx,:), 'filled', 'DisplayName', upper(pol));

    % 4) Find coarse peaks
    [pks, locs] = findpeaks(y, 'MinPeakHeight', threshold);

    thisPeaks = [];
    for k = 1:numel(locs)
        centerIdx = locs(k);
        i1 = max(1,centerIdx - searchHalfWidth);
        i2 = min(numel(y), centerIdx + searchHalfWidth);

        x_fit = freq_Hz(i1:i2).';      % column vector
        y_fit = y(i1:i2);

        if numel(x_fit) < 5
            continue;
        end

        % 5) Lightly smooth & fit a single Gaussian
        y_sm = smooth(y_fit, 5);
        y_sm = y_sm(:);  
        w = sqrt(y_sm);
        opts = fitoptions('poly4');
        opts.Normalize = 'on';
        opts.Robust    = 'LAR';   % ensure column
        opts.Weights = w;
            % … after you do:
        spOpts = fitoptions('smoothingspline','SmoothingParam',0.9);
        g      = fit(x_fit, y_sm, 'smoothingspline', spOpts);
        %g    = fit(x_fit, y_sm, 'poly4', opts);
        xg   = linspace(min(x_fit), max(x_fit), 20000);
        yg   = feval(g, xg);

        % 1) figure out where the raw data peak is in your window
        %    (centerIdx is the locs(k) you already have)
        x_data_peak = freq_Hz(centerIdx);
        y_data_peak = y(centerIdx);

        % 2) what does your poly4 predict at that exact x?
        y_fit_at_data_peak = feval(g, x_data_peak);

        % 3) compute the offset needed to lift the poly up onto the data
        offset = y_data_peak - y_fit_at_data_peak;

        % 4) apply that same offset to your entire fitted curve
        yg = yg + offset;

        % 5) re‑find the maximum on the *offset* curve
        [true_peak_amp, imax] = max(yg);
        true_peak_freq        = xg(imax);

        % 6) plot it
        %plot(xg, yg, 'Color',colors(idx,:), 'LineWidth',1.5);
        xline(true_peak_freq, '--', 'Color',colors(idx,:), 'LineWidth',1.5);

        fprintf("Pol: %s | Peak Freq: %.6f MHz | Amp: %.3f\n", ...
            upper(pol), true_peak_freq/1e6, true_peak_amp);
        thisPeaks(end+1) = true_peak_freq;

    end

    peakMap(pol) = thisPeaks;
end

% ——— Final formatting of the plot ———
xlabel('Fréquence (Hz)', 'FontSize', 14, 'Interpreter','latex');
ylabel('Amplitude (u.a.)', 'FontSize', 14, 'Interpreter','latex');
legend('show', 'Location', 'best');
grid on;

% ——— Cluster all peaks into “zones” and build result table ———
allPeaks = [];
for idx = 1:numel(polarization_states)
    pol   = polarization_states{idx};
    peaks = peakMap(pol);
    for pf = peaks
        allPeaks(end+1, :) = [pf, idx];   % [frequency_Hz, pol_index]
    end
end
allPeaks = sortrows(allPeaks, 1);

% grouping tolerance: 20 MHz

zoneList     = {};
currZone     = allPeaks(1, :);
currRef      = allPeaks(1, 1);

for r = 2:size(allPeaks,1)
    f = allPeaks(r,1);
    if abs(f - currRef) <= tolerance_Hz
        currZone = [currZone; allPeaks(r,:)];
    else
        zoneList{end+1} = currZone;
        currZone        = allPeaks(r,:);
        currRef         = f;
    end
end
zoneList{end+1} = currZone;

% Build a cell array of [Zone, Pol, Freq_MHz, Delay_kHz]
results = {};
for z = 1:numel(zoneList)
    zone = zoneList{z};

    % pick the ‘v’ peak as reference
    ref_freq = NaN;
    for row = 1:size(zone,1)
        if strcmp(polarization_states{zone(row,2)}, 'v')
            ref_freq = zone(row,1);
            break;
        end
    end

    for idx = 1:numel(polarization_states)
        pol      = polarization_states{idx};
        freq_pos = NaN;
        for row = 1:size(zone,1)
            if strcmp(polarization_states{zone(row,2)}, pol)
                freq_pos = zone(row,1);
                break;
            end
        end
        delay_Hz = NaN;
        if ~isnan(ref_freq) && ~isnan(freq_pos)
            delay_Hz = freq_pos - ref_freq;
        end

        results(end+1, :) = { ...
            z, ...
            upper(pol), ...
            freq_pos   / 1e6, ...  % MHz
            delay_Hz   / 1e3       % kHz
        };
    end
end

% convert to table and write CSV
resultTable = cell2table(results, ...
    'VariableNames', {'Zone','Polarisation','Frequence_MHz','Delai_kHz'});
disp(resultTable);
writetable(resultTable, 'delais_par_rapport_V.csv');
