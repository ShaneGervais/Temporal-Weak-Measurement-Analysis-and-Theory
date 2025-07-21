%% Imaginary part data analysis
clc; clear; close all;

%% --- DATA PARAMETERS ---
folder = '../imag_0127_ult';       
sampling_interval_ps = 4;             % Sampling interval in picoseconds
time_step = sampling_interval_ps * 1e-12;  % Convert to seconds

%% --- LOAD FILES ---
files     = dir(fullfile(folder, '*.csv'));
angle_map = containers.Map('KeyType','double','ValueType','any');



for k = 1:numel(files)
    fname = files(k).name;
    t     = regexp(fname, '_(\d+)_deg_', 'tokens');
    if isempty(t), continue; end
    angle = str2double(t{1}{1});
    
    if isKey(angle_map, angle)
        lst = angle_map(angle);
    else
        lst = {};
    end
    lst{end+1} = fullfile(folder, fname);
    angle_map(angle) = lst;
end

%% --- FREQUENCY EXTRACTION ---
angles         = sort(cell2mat(keys(angle_map)))';
chirp_freq_MHz = nan(size(angles));
%% --- FREQUENCY EXTRACTION ---
angles         = sort(cell2mat(keys(angle_map)))';

% Fix: define ref_idx early
ref_idx = find(angles==48,1);
if isempty(ref_idx)
    error('No 48° data point found for normalization.');
end

chirp_freq_MHz = nan(size(angles));

figure;
for i = 1:numel(angles)
    files_i    = angle_map(angles(i));
    all_traces = [];
    for j = 1:numel(files_i)
        opts           = detectImportOptions(files_i{j});
        opts.DataLines = [3, Inf];
        M              = readmatrix(files_i{j}, opts);
        all_traces(:,end+1) = M(:,2);
    end
    avg_trace = mean(all_traces, 2);
    %plot(avg_trace)
    %hold on;
    N         = numel(avg_trace);

    % FFT → power & freq axis
    S = fft(avg_trace);
    P = abs(S).^2;

    %plot(P);
    %hold on;

    f = fftshift(fftfreq(N, time_step)).';
    P = fftshift(P);
    idx = f >= 0;
    chirp_freq_MHz(i) = sum(f(idx).*P(idx)) / sum(P(idx)) * 1e-6;
end
hold off;

%% --- SCALE STATES & NORMALIZE to V‑state ---
scaled_states  = angles * 2;
ref_idx        = find(angles==48,1);
if isempty(ref_idx)
    error('No 48° data point found for normalization.');
end
f_ref         = chirp_freq_MHz(ref_idx);
norm_freq_MHz = (chirp_freq_MHz - f_ref);

%% --- COSINE FIT ---
ft = fittype('A*(cos(B*x + C)).^2 ', ...
    'independent','x','coefficients',{'A','B','C'});
start = [ ...
    (max(norm_freq_MHz)-min(norm_freq_MHz))/2, ...   % A
    2*pi/(max(scaled_states)-min(scaled_states)), ...% B
    0 ...                                            % C
    %mean(norm_freq_MHz)
    ];                             % D
[cos_fit, gof] = fit(scaled_states, norm_freq_MHz, ft, 'StartPoint', start);

%% --- PLOTTING ---
figure; hold on;
plot(scaled_states, norm_freq_MHz, 'ko', 'LineWidth',1.5, 'DisplayName','Data');
xFit = linspace(min(scaled_states), max(scaled_states), 200);
yFit = feval(cos_fit, xFit);
plot(xFit, yFit, 'r--', 'LineWidth',2, 'DisplayName','Cosine Fit');
hold off; grid on;

xlabel('Input state, $|\psi(\theta)>$', ...
    'Interpreter','latex','FontSize',14);
ylabel('Weak value, $\Im\langle\hat{A}_W\rangle$ (a.u.)', ...
    'Interpreter','latex','FontSize',14);

legend('Location','best','Interpreter','none');

%% --- DISPLAY TABLE ---

T = table(angles(:), scaled_states(:), chirp_freq_MHz(:), norm_freq_MHz(:), ...
    'VariableNames',{'Angle_deg','ScaledState','Freq_MHz','NormFreq_MHz'});
disp(T);

%% --- FONCTIONS ---
function f = fftfreq(n, d)
    val = 1/(n*d);
    if mod(n,2)==0
        k = -n/2 : n/2-1;
    else
        k = -(n-1)/2 : (n-1)/2;
    end
    f = k * val;
end
