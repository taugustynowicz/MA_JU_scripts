%% Initial cleaning
clc, clear, close all

%% Define paths and load data

% Initialize FieldTrip
addpath('example_path');
ft_defaults;

% Define paths and file lists
output_dir = fullfile(pwd, 'output');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end
control_path = fullfile(pwd, 'mats', 'chirp_hap12_stderp_std', 'controls');
patient_path = fullfile(pwd, 'mats', 'chirp_hap12_stderp_std', 'patients');
control_files = dir(fullfile(control_path, '*.mat'));
patient_files = dir(fullfile(patient_path, '*.mat'));

% Load control group data
control_data = cell(1, length(control_files));
for i = 1:length(control_files)
    try
        tmp = load(fullfile(control_path, control_files(i).name));
        control_data{i} = tmp.preprocessed;
        % Verify trial count
        if isfield(control_data{i}, 'trial') && iscell(control_data{i}.trial)
            fprintf('Control data %d (%s): %d trials\n', i, control_files(i).name, length(control_data{i}.trial));
        else
            warning('Control data %d (%s) is missing trial field or is not a cell array!', i, control_files(i).name);
        end
    catch e
        warning('Failed to load control file %s: %s', control_files(i).name, e.message);
        control_data{i} = [];
    end
end

% Load patient group data
patient_data = cell(1, length(patient_files));
for i = 1:length(patient_files)
    try
        tmp = load(fullfile(patient_path, patient_files(i).name));
        patient_data{i} = tmp.preprocessed;
        % Verify trial count
        if isfield(patient_data{i}, 'trial') && iscell(patient_data{i}.trial)
            fprintf('Patient data %d (%s): %d trials\n', i, patient_files(i).name, length(patient_data{i}.trial));
        else
            warning('Patient data %d (%s) is missing trial field or is not a cell array!', i, patient_files(i).name);
        end
    catch e
        warning('Failed to load patient file %s: %s', patient_files(i).name, e.message);
        patient_data{i} = [];
    end
end

% Remove empty datasets
control_data = control_data(~cellfun(@isempty, control_data));
patient_data = patient_data(~cellfun(@isempty, patient_data));

%% Apply Low-Pass Filter at 40 Hz

cfg = [];
cfg.lpfilter = 'yes'; % Enable low-pass filter
cfg.lpfreq = 40; % Low-pass frequency cutoff at 40 Hz
cfg.lpfiltord = 4; % Filter order (e.g., 4th order Butterworth)
cfg.lpfilttype = 'but'; % Butterworth filter type
cfg.lpfiltdir = 'twopass'; % Two-pass filtering for zero-phase shift

% Filter control group data
for i = 1:length(control_data)
    fprintf('Applying 40 Hz low-pass filter to control data %d (%s)...\n', i, control_files(i).name);
    if isfield(control_data{i}, 'trial') && iscell(control_data{i}.trial)
        control_data{i} = ft_preprocessing(cfg, control_data{i});
    else
        warning('Control data %d (%s) is missing trial field or is not a cell array. Skipping filtering.', i, control_files(i).name);
    end
end

% Filter patient group data
for i = 1:length(patient_data)
    fprintf('Applying 40 Hz low-pass filter to patient data %d (%s)...\n', i, patient_files(i).name);
    if isfield(patient_data{i}, 'trial') && iscell(patient_data{i}.trial)
        patient_data{i} = ft_preprocessing(cfg, patient_data{i});
    else
        warning('Patient data %d (%s) is missing trial field or is not a cell array. Skipping filtering.', i, patient_files(i).name);
    end
end

%% Compute ERPs with channel averaging

% Define channels to average
selected_channels = {'FC1', 'FCz', 'C1', 'Cz', 'C2', 'CP1', 'CPz', 'CP2', 'CP4', 'P3', 'P1', 'Pz', 'P2', 'P4'};

cfg = [];
cfg.keeptrials = 'no'; % Average across trials
cfg.channel = selected_channels; % Select specified channels

for i = 1:length(control_data)
    % Check for missing channels
    missing_chans = setdiff(selected_channels, control_data{i}.label);
    if ~isempty(missing_chans)
        warning('Control data %d missing channels: %s', i, strjoin(missing_chans, ', '));
    end
    if isfield(control_data{i}, 'avg')
        fprintf('Control data %d is already time-locked. Skipping ft_timelockanalysis.\n', i);
    else
        % Select channels and compute ERP
        control_data{i} = ft_timelockanalysis(cfg, control_data{i});
        % Average across selected channels
        chan_idx = ismember(control_data{i}.label, selected_channels);
        control_data{i}.avg = mean(control_data{i}.avg(chan_idx, :), 1);
        if isfield(control_data{i}, 'var')
            control_data{i}.var = mean(control_data{i}.var(chan_idx, :), 1); % Average variance
        end
        if isfield(control_data{i}, 'dof')
            control_data{i}.dof = mean(control_data{i}.dof(chan_idx, :), 1); % Average dof
        end
        control_data{i}.label = {'avg_channel'};
        control_data{i}.dimord = 'chan_time';
    end
end
for i = 1:length(patient_data)
    % Check for missing channels
    missing_chans = setdiff(selected_channels, patient_data{i}.label);
    if ~isempty(missing_chans)
        warning('Patient data %d missing channels: %s', i, strjoin(missing_chans, ', '));
    end
    if isfield(patient_data{i}, 'avg')
        fprintf('Patient data %d is already time-locked. Skipping ft_timelockanalysis.\n', i);
    else
        % Select channels and compute ERP
        patient_data{i} = ft_timelockanalysis(cfg, patient_data{i});
        % Average across selected channels
        chan_idx = ismember(patient_data{i}.label, selected_channels);
        patient_data{i}.avg = mean(patient_data{i}.avg(chan_idx, :), 1);
        if isfield(patient_data{i}, 'var')
            patient_data{i}.var = mean(patient_data{i}.var(chan_idx, :), 1); % Average variance
        end
        if isfield(patient_data{i}, 'dof')
            patient_data{i}.dof = mean(patient_data{i}.dof(chan_idx, :), 1); % Average dof
        end
        patient_data{i}.label = {'avg_channel'};
        patient_data{i}.dimord = 'chan_time';
    end
end

% Compute Grand Averages for Visualization
cfg = [];
cfg.keepindividual = 'no'; % Average across subjects
control_grandavg = ft_timelockgrandaverage(cfg, control_data{:});
patient_grandavg = ft_timelockgrandaverage(cfg, patient_data{:});

%% Load channel labels from the first control file
tmp = load(fullfile(control_path, control_files(1).name));
chan_labels = tmp.preprocessed.label; % 31×1 cell array
fprintf('Data channel labels (%d channels):\n', length(chan_labels));
disp(chan_labels);

% Load BioSemi 64 layout
cfg_layout = [];
cfg_layout.layout = 'biosemi64.lay';
lay = ft_prepare_layout(cfg_layout);

% Match and subset to 31 channels
keep_idx = ismember(strtrim(lower(lay.label)), strtrim(lower(chan_labels)));
if sum(keep_idx) ~= 31
    error('Expected 31 channels, but %d were matched', sum(keep_idx));
end

% Create custom layout structure
cfg = [];
cfg.layout = struct();
cfg.layout.label = lay.label(keep_idx); % 31 channel labels
cfg.layout.pos = lay.pos(keep_idx, :); % Channel positions
cfg.layout.width = zeros(size(lay.width(keep_idx))); % Set marker width to 0
cfg.layout.height = zeros(size(lay.height(keep_idx))); % Set marker height to 0
cfg.showlabels = 'yes'; % Ensure labels are displayed
cfg.fontsize = 12; % Set font size for labels
cfg.color = 'k'; % Black labels
cfg.box = 'no'; % No boxes around channels
cfg.outline = 'yes'; % Show head outline
cfg.mask = 'yes'; % Show head mask
cfg.style = '2d'; % Explicitly set to 2D to avoid any 3D marker rendering

%% Perform cluster-based statistics on a single 0–4000 ms window

cfg = [];
cfg.method = 'montecarlo';
cfg.statistic = 'indepsamplesT';
cfg.correctm = 'cluster';
cfg.alpha = 0.05;
cfg.numrandomization = 1000;
cfg.design = [ones(1, length(control_data)) 2*ones(1, length(patient_data))]';
cfg.ivar = 1;
cfg.latency = [0 4]; % Single time window from 0 to 4000 ms
cfg.channel = {'avg_channel'};

% Perform permutation test
stats_result = ft_timelockstatistics(cfg, control_data{:}, patient_data{:});

fprintf('Time window 0–4 s:\n');
if isfield(stats_result, 'prob') && ~isempty(stats_result.prob)
    fprintf('  Cluster-based p-value = %.3f\n', min(stats_result.prob));
    if any(stats_result.prob < 0.05)
        fprintf('  Significant difference detected (p < 0.05).\n');
    else
        fprintf('  No significant difference (p >= 0.05).\n');
    end
else
    fprintf('  No valid statistical results for this window.\n');
end

% Print significant cluster time ranges
if isfield(stats_result, 'posclusters') && ~isempty(stats_result.posclusters)
    for c = 1:length(stats_result.posclusters)
        if stats_result.posclusters(c).prob < 0.05
            sig_time = stats_result.time(stats_result.posclusterslabelmat == c);
            if ~isempty(sig_time)
                fprintf('  Significant positive cluster %d: %.3f–%.3f s (p = %.3f)\n', ...
                    c, min(sig_time), max(sig_time), stats_result.posclusters(c).prob);
            end
        end
    end
end
if isfield(stats_result, 'negclusters') && ~isempty(stats_result.negclusters)
    for c = 1:length(stats_result.negclusters)
        if stats_result.negclusters(c).prob < 0.05
            sig_time = stats_result.time(stats_result.negclusterslabelmat == -c);
            if ~isempty(sig_time)
                fprintf('  Significant negative cluster %d: %.3f–%.3f s (p = %.3f)\n', ...
                    c, min(sig_time), max(sig_time), stats_result.negclusters(c).prob);
            end
        end
    end
end

% Save results
%save(fullfile(output_dir, 'chirp_hap1_erp_stats_avg_channels.mat'), 'stats_result', 'control_grandavg', 'patient_grandavg', '-v7.3');

% Create a figure with custom size (e.g., 1200x600 pixels)
figure('Position', [100, 100, 1200, 600]);
cfg = [];
cfg.channel = 'avg_channel';
cfg.xlim = [-0.2 4];
cfg.ylim = [-1.5 1.5];
cfg.showgrid = 'yes';
ft_singleplotER(cfg, control_grandavg, patient_grandavg);
hold on;

% Plot SEM as shaded areas
time = control_grandavg.time;
control_sem = sqrt(control_grandavg.var(strcmp(control_grandavg.label, cfg.channel), :)) / sqrt(length(control_data));
patient_sem = sqrt(patient_grandavg.var(strcmp(patient_grandavg.label, cfg.channel), :)) / sqrt(length(patient_data));
control_mean = control_grandavg.avg(strcmp(control_grandavg.label, cfg.channel), :);
patient_mean = patient_grandavg.avg(strcmp(patient_grandavg.label, cfg.channel), :);

% Shaded area for control group ±2 SEM
fill([time, fliplr(time)], [control_mean + 2*control_sem, fliplr(control_mean - 2*control_sem)], ...
    'b', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
% Shaded area for patient group ±2 SEM
fill([time, fliplr(time)], [patient_mean + 2*patient_sem, fliplr(patient_mean - 2*patient_sem)], ...
    'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none');

% Highlight significant clusters
if isfield(stats_result, 'posclusters') && ~isempty(stats_result.posclusters)
    for c = 1:length(stats_result.posclusters)
        if stats_result.posclusters(c).prob < 0.05
            sig_time = stats_result.time(stats_result.posclusterslabelmat == c);
            if ~isempty(sig_time)
                x = [min(sig_time), max(sig_time), max(sig_time), min(sig_time)];
                y = [-1.5, -1.5, 1.5, 1.5]; % Adjusted to match ylim
                fill(x, y, [0.8 0.8 0.8], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
            end
        end
    end
end
if isfield(stats_result, 'negclusters') && ~isempty(stats_result.negclusters)
    for c = 1:length(stats_result.negclusters)
        if stats_result.negclusters(c).prob < 0.05
            sig_time = stats_result.time(stats_result.negclusterslabelmat == -c);
            if ~isempty(sig_time)
                x = [min(sig_time), max(sig_time), max(sig_time), min(sig_time)];
                y = [-1.5, -1.5, 1.5, 1.5]; % Adjusted to match ylim
                fill(x, y, [0.8 0.8 0.8], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
            end
        end
    end
end

% Add vertical line at t=0 for stimulus onset
plot([0 0], [-1.5 1.5], 'k-', 'LineWidth', 1);

hold off;

title('ERP Comparison between the Patient and Control Group');
xlabel('Time (s)');
ylabel('Amplitude (µV)');
legend({'Control', 'Patient', 'Control ±2 SEM', 'Patient ±2 SEM', 'Significant Clusters'});
%saveas(gcf, fullfile(output_dir, 'erp_comparison_avg_channels_significant_2sem.eps'), 'epsc');
