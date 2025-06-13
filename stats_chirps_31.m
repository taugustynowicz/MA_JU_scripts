%% Initial cleaning
clc, clear, close all

%% Define paths and load data

% Initialize FieldTrip
addpath('d:\marek\github-repositories\fieldtrip\');
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
    tmp = load(fullfile(control_path, control_files(i).name));
    control_data{i} = tmp.preprocessed;
    % Verify trial count
    if isfield(control_data{i}, 'trial') && iscell(control_data{i}.trial)
        fprintf('Control data %d (%s): %d trials\n', i, control_files(i).name, length(control_data{i}.trial));
    else
        warning('Control data %d (%s) is missing trial field or is not a cell array!', i, control_files(i).name);
    end
end

% Load patient group data
patient_data = cell(1, length(patient_files));
for i = 1:length(patient_files)
    tmp = load(fullfile(patient_path, patient_files(i).name));
    patient_data{i} = tmp.preprocessed;
    % Verify trial count
    if isfield(patient_data{i}, 'trial') && iscell(patient_data{i}.trial)
        fprintf('Patient data %d (%s): %d trials\n', i, patient_files(i).name, length(patient_data{i}.trial));
    else
        warning('Patient data %d (%s) is missing trial field or is not a cell array!', i, patient_files(i).name);
    end
end

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

%% Compute ERPs

cfg = [];
cfg.keeptrials = 'no'; % Average across trials
% Uncomment the next line if you need to select STD stimuli
% cfg.trials = find(control_data{i}.trialinfo == 1); % Example for STD stimuli

for i = 1:length(control_data)
    if isfield(control_data{i}, 'avg')
        fprintf('Control data %d is already time-locked. Skipping ft_timelockanalysis.\n', i);
    else
        control_data{i} = ft_timelockanalysis(cfg, control_data{i});
    end
end
for i = 1:length(patient_data)
    if isfield(patient_data{i}, 'avg')
        fprintf('Patient data %d is already time-locked. Skipping ft_timelockanalysis.\n', i);
    else
        patient_data{i} = ft_timelockanalysis(cfg, patient_data{i});
    end
end

% Compute Grand Averages for Visualization
cfg = [];
cfg.keepindividual = 'no'; % Average across subjects
control_grandavg = ft_timelockgrandaverage(cfg, control_data{:});
patient_grandavg = ft_timelockgrandaverage(cfg, patient_data{:});

%% Custom 31-channel layout

% Load channel labels from the first control file
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
% Verify that exactly 31 channels are matched
if sum(keep_idx) ~= 31
    error('Expected 31 channels, but %d were matched', sum(keep_idx));
end

% Create custom layout structure
cfg = [];
cfg.layout = struct();
cfg.layout.label = lay.label(keep_idx);
cfg.layout.pos = lay.pos(keep_idx, :);
cfg.layout.width = lay.width(keep_idx);
cfg.layout.height = lay.height(keep_idx);
cfg.showlabels = 'yes'; % Display channel labels on plot

%% Perform cluster-based statistics

% Prepare neighbours
cfg_neighb = [];
cfg_neighb.method = 'distance';
cfg_neighb.layout = cfg.layout;
neighbours = ft_prepare_neighbours(cfg_neighb);
fprintf('Neighbours structure created with %d channels.\n', length(neighbours));

% Cluster-based permutation test
cfg = [];
cfg.method = 'montecarlo'; % Permutation test
cfg.statistic = 'indepsamplesT'; % Independent samples t-test
cfg.correctm = 'cluster'; % Cluster-based correction
cfg.clusteralpha = 0.05; % Alpha for cluster entry
cfg.clusterstatistic = 'maxsum'; % Statistic for cluster selection
cfg.minnbchan = 2; % Minimum number of neighboring channels
cfg.neighbours = neighbours; % Neighbours structure
cfg.tail = 0; % Two-tailed test
cfg.clustertail = 0;
cfg.alpha = 0.05; % Alpha for final test
cfg.numrandomization = 1000; % Number of permutations
cfg.design = [ones(1, length(control_data)) 2*ones(1, length(patient_data))]'; % Group labels
cfg.ivar = 1; % Independent variable (group)
cfg.latency = [0 4]; % Explicitly set to 0 to 4000 ms
cfg.channel = {'all'}; % All channels

% Run statistics
stat = ft_timelockstatistics(cfg, control_data{:}, patient_data{:});

% Report clusters
if isfield(stat, 'posclusters') && ~isempty(stat.posclusters)
    fprintf('Found %d positive clusters:\n', length(stat.posclusters));
    for i = 1:length(stat.posclusters)
        fprintf('Positive cluster %d: p = %.3f\n', i, stat.posclusters(i).prob);
        if stat.posclusters(i).prob < 0.05 && isfield(stat, 'posclusterslabelmat') && ~isempty(stat.posclusterslabelmat)
            cluster_mask = stat.posclusterslabelmat == i;
            if any(cluster_mask(:))
                channels = stat.label(any(cluster_mask, 2));
                times = stat.time(any(cluster_mask, 1));
                fprintf('  Channels: %s\n', strjoin(channels, ', '));
                fprintf('  Time range: %.3f–%.3f s\n', min(times), max(times));
            else
                fprintf('  No valid channels or time points in cluster %d.\n', i);
            end
        end
    end
else
    fprintf('No positive clusters detected.\n');
end
if isfield(stat, 'negclusters') && ~isempty(stat.negclusters)
    fprintf('Found %d negative clusters:\n', length(stat.negclusters));
    for i = 1:length(stat.negclusters)
        fprintf('Negative cluster %d: p = %.3f\n', i, stat.negclusters(i).prob);
        if stat.negclusters(i).prob < 0.05 && isfield(stat, 'negclusterslabelmat') && ~isempty(stat.negclusterslabelmat)
            cluster_mask = stat.negclusterslabelmat == i;
            if any(cluster_mask(:))
                channels = stat.label(any(cluster_mask, 2));
                times = stat.time(any(cluster_mask, 1));
                fprintf('  Channels: %s\n', strjoin(channels, ', '));
                fprintf('  Time range: %.3f–%.3f s\n', min(times), max(times));
            else
                fprintf('  No valid channels or time points in cluster %d.\n', i);
            end
        end
    end
else
    fprintf('No negative clusters detected.\n');
end

% Print structure of stat for debugging
fprintf('\nStructure of stat:\n');
disp(fieldnames(stat));

% Print details of stat.posclusterslabelmat
if isfield(stat, 'posclusterslabelmat') && ~isempty(stat.posclusterslabelmat)
    fprintf('posclusterslabelmat dimensions: %s\n', mat2str(size(stat.posclusterslabelmat)));
    fprintf('posclusterslabelmat unique values: %s\n', mat2str(unique(stat.posclusterslabelmat(:))));
else
    fprintf('posclusterslabelmat is missing or empty.\n');
end

% Print details of stat.negclusterslabelmat
if isfield(stat, 'negclusterslabelmat') && ~isempty(stat.negclusterslabelmat)
    fprintf('negclusterslabelmat dimensions: %s\n', mat2str(size(stat.negclusterslabelmat)));
    fprintf('negclusterslabelmat unique values: %s\n', mat2str(unique(stat.negclusterslabelmat(:))));
else
    fprintf('negclusterslabelmat is missing or empty.\n');
end

% Print stat.time and stat.label for reference
if isfield(stat, 'time')
    fprintf('stat.time: %d time points, range: [%.3f, %.3f] s\n', length(stat.time), min(stat.time), max(stat.time));
else
    fprintf('stat.time is missing.\n');
end
if isfield(stat, 'label')
    fprintf('stat.label: %d channels\n', length(stat.label));
    disp(stat.label);
else
    fprintf('stat.label is missing.\n');
end

% Visualize significant positive cluster with SEM
if isfield(stat, 'posclusters') && ~isempty(stat.posclusters) && any([stat.posclusters.prob] < 0.05)
    % Find the first significant positive cluster
    sig_cluster_idx = find([stat.posclusters.prob] < 0.05, 1, 'first');
    if ~isempty(sig_cluster_idx) && isfield(stat, 'posclusterslabelmat') && ~isempty(stat.posclusterslabelmat)
        cluster_mask = stat.posclusterslabelmat == sig_cluster_idx;
        if any(cluster_mask(:))
            % Select channels
            selected_channels = stat.label(any(cluster_mask, 2));
            % Verify channels exist in grand averages
            valid_channels = intersect(selected_channels, control_grandavg.label);
            if ~isempty(valid_channels)
                cfg = [];
                cfg.channel = valid_channels;
                cfg.xlim = [-0.2 4]; % Explicitly set to 0 to 4000 ms
                cfg.ylim = [-2 2];
                if all(isfinite(cfg.xlim))
                    figure('Position', [100, 100, 1200, 600]);
                    ft_singleplotER(cfg, control_grandavg, patient_grandavg);
                    hold on; % Allow adding to the plot
                    
                    % Compute SEM for selected channels
                    control_mean = mean(control_grandavg.avg(ismember(control_grandavg.label, valid_channels), :), 1);
                    patient_mean = mean(patient_grandavg.avg(ismember(patient_grandavg.label, valid_channels), :), 1);
                    if isfield(control_grandavg, 'var') && isfield(patient_grandavg, 'var')
                        control_sem = sqrt(mean(control_grandavg.var(ismember(control_grandavg.label, valid_channels), :), 1)) / sqrt(length(control_data));
                        patient_sem = sqrt(mean(patient_grandavg.var(ismember(patient_grandavg.label, valid_channels), :), 1)) / sqrt(length(patient_data));
                        
                        % Plot SEM as shaded areas
                        time = control_grandavg.time;
                        fill([time, fliplr(time)], [control_mean + 2*control_sem, fliplr(control_mean - 2*control_sem)], ...
                            'b', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
                        fill([time, fliplr(time)], [patient_mean + 2*patient_sem, fliplr(patient_mean - 2*patient_sem)], ...
                            'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
                    else
                        fprintf('Variance field missing in grand averages for positive cluster %d SEM plotting.\n', sig_cluster_idx);
                    end
                    
                    % Find significant time bins
                    sig_times = any(cluster_mask, 1); % Logical array for significant time points
                    time_points = stat.time; % Time vector
                    ylim = cfg.ylim; % Y-axis limits for shading
                    
                    % Identify contiguous significant time segments
                    diff_sig = diff([0 sig_times 0]); % Pad with zeros to detect edges
                    starts = find(diff_sig == 1); % Start indices of significant segments
                    ends = find(diff_sig == -1) - 1; % End indices of significant segments

                    % After identifying significant segments
                    for seg = 1:length(starts)
                        t_start = time_points(starts(seg));
                        t_end = time_points(ends(seg));
                        fprintf('Significant time segment %d: %.3f s to %.3f s\n', seg, t_start, t_end);
                    end
                    
                    % Plot grey patches for each significant segment
                    for seg = 1:length(starts)
                        t_start = time_points(starts(seg));
                        t_end = time_points(ends(seg));
                        patch([t_start t_end t_end t_start], [ylim(1) ylim(1) ylim(2) ylim(2)], ...
                            [0.7 0.7 0.7], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
                    end
                    
                    % Add vertical line at t=0 for stimulus onset
                    plot([0 0], cfg.ylim, 'k-', 'LineWidth', 1);
                    
                    hold off;
                    title(sprintf('ERP for Significant Positive Cluster %d (p = %.3f)', sig_cluster_idx, stat.posclusters(sig_cluster_idx).prob));
                    xlabel('Time (s)');
                    ylabel('Amplitude (µV)');
                    legend({'Control', 'Patient', 'Control ±2 SEM', 'Patient ±2 SEM', 'Significant Windows'});
                    saveas(gcf, fullfile(output_dir, sprintf('significantfeat_cluster_erp_pos_%d_0to4000ms.eps', sig_cluster_idx)), 'epsc');
                else
                    fprintf('Invalid time range for positive cluster %d.\n', sig_cluster_idx);
                end
            else
                fprintf('No valid channels found for positive cluster %d in grand averages.\n', sig_cluster_idx);
            end
        else
            fprintf('No valid cluster mask for positive cluster %d.\n', sig_cluster_idx);
        end
    else
        fprintf('No significant positive clusters found with p < 0.05.\n');
    end
else
    fprintf('No positive clusters detected.\n');
end

% Visualize significant negative clusters with SEM
if isfield(stat, 'negclusters') && ~isempty(stat.negclusters) && any([stat.negclusters.prob] < 0.05)
    % Find the first significant negative cluster
    sig_cluster_idx = find([stat.negclusters.prob] < 0.05, 1, 'first');
    if ~isempty(sig_cluster_idx) && isfield(stat, 'negclusterslabelmat') && ~isempty(stat.negclusterslabelmat)
        cluster_mask = stat.negclusterslabelmat == sig_cluster_idx;
        if any(cluster_mask(:))
            % Select channels and time points
            selected_channels = stat.label(any(cluster_mask, 2));
            selected_times = stat.time(any(cluster_mask, 1));
            % Verify channels exist in grand averages
            valid_channels = intersect(selected_channels, control_grandavg.label);
            if ~isempty(valid_channels)
                cfg = [];
                cfg.channel = valid_channels;
                cfg.xlim = [min(selected_times), max(selected_times)]; % Use cluster-specific time range
                cfg.ylim = [-1.2 1.2]; % Set y-axis limits to -1.2 to +1.2 µV
                if all(isfinite(cfg.xlim))
                    figure('Position', [100, 100, 1200, 600]);
                    ft_singleplotER(cfg, control_grandavg, patient_grandavg);
                    hold on; % Allow adding to the plot
                    
                    % Compute SEM for selected channels
                    control_mean = mean(control_grandavg.avg(ismember(control_grandavg.label, valid_channels), :), 1);
                    patient_mean = mean(patient_grandavg.avg(ismember(patient_grandavg.label, valid_channels), :), 1);
                    if isfield(control_grandavg, 'var') && isfield(patient_grandavg, 'var')
                        control_sem = sqrt(mean(control_grandavg.var(ismember(control_grandavg.label, valid_channels), :), 1)) / sqrt(length(control_data));
                        patient_sem = sqrt(mean(patient_grandavg.var(ismember(patient_grandavg.label, valid_channels), :), 1)) / sqrt(length(patient_data));
                        
                        % Plot SEM as shaded areas
                        time = control_grandavg.time;
                        fill([time, fliplr(time)], [control_mean + 2*control_sem, fliplr(control_mean - 2*control_sem)], ...
                            'b', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
                        fill([time, fliplr(time)], [patient_mean + 2*patient_sem, fliplr(patient_mean - 2*patient_sem)], ...
                            'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
                    else
                        fprintf('Variance field missing in grand averages for negative cluster %d SEM plotting.\n', sig_cluster_idx);
                    end
                    
                    % Add vertical line at t=0 for stimulus onset
                    plot([0 0], cfg.ylim, 'k-', 'LineWidth', 1);
                    
                    hold off;
                    title(sprintf('ERP for Significant Negative Cluster %d (p = %.3f, %.3f–%.3f ms)', sig_cluster_idx, stat.negclusters(sig_cluster_idx).prob, min(selected_times), max(selected_times)));
                    xlabel('Time (s)');
                    ylabel('Amplitude (µV)');
                    legend({'Control', 'Patient', 'Control ±2 SEM', 'Patient ±2 SEM'});
                    saveas(gcf, fullfile(output_dir, sprintf('significant_cluster_erp_neg_%d_0to4000ms.eps', sig_cluster_idx)), 'epsc');
                else
                    fprintf('Invalid time range for negative cluster %d.\n', sig_cluster_idx);
                end
            else
                fprintf('No valid channels found for negative cluster %d in grand averages.\n', sig_cluster_idx);
            end
        else
            fprintf('No valid cluster mask for negative cluster %d.\n', sig_cluster_idx);
        end
    else
        fprintf('No significant negative clusters found with p < 0.05.\n');
    end
else
    fprintf('No negative clusters detected.\n');
end

%% Default ERP Plot with SEM (All Channels)
% Plot grand average ERP for all channels, regardless of cluster significance
cfg = [];
cfg.channel = {'all'};
cfg.xlim = [-0.2 4];
cfg.ylim = [-2 2];
figure('Position', [100, 100, 1200, 600]);
ft_singleplotER(cfg, control_grandavg, patient_grandavg);
hold on;

% Compute SEM for all channels
control_mean = mean(control_grandavg.avg, 1);
patient_mean = mean(patient_grandavg.avg, 1);
if isfield(control_grandavg, 'var') && isfield(patient_grandavg, 'var')
    control_sem = sqrt(mean(control_grandavg.var, 1)) / sqrt(length(control_data));
    patient_sem = sqrt(mean(patient_grandavg.var, 1)) / sqrt(length(patient_data));
    
    % Plot SEM as shaded areas
    time = control_grandavg.time;
    fill([time, fliplr(time)], [control_mean + 2*control_sem, fliplr(control_mean - 2*control_sem)], ...
        'b', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
    fill([time, fliplr(time)], [patient_mean + 2*patient_sem, fliplr(patient_mean - 2*patient_sem)], ...
        'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
else
    fprintf('Variance field missing in grand averages for default ERP SEM plotting.\n');
end

% Add vertical line at t=0 for stimulus onset
plot([0 0], cfg.ylim, 'k-', 'LineWidth', 1);

hold off;
title('Grand Average ERP Across All Channels');
xlabel('Time (s)');
ylabel('Amplitude (µV)');
legend({'UWS', 'MCSe', 'UWS ±2 SEM', 'MCSe ±2 SEM'});
%saveas(gcf, fullfile(output_dir, 'grand_average_erp_all_channels_0to4000ms.eps'), 'epsc');

%% Topographic Plot of Significant Clusters

% Compute difference ERP for visualization
cfg = [];
cfg.operation = 'subtract';
cfg.parameter = 'avg';
diff_grandavg = ft_math(cfg, control_grandavg, patient_grandavg);

% Plot topographic map for significant positive cluster
if isfield(stat, 'posclusters') && ~isempty(stat.posclusters) && any([stat.posclusters.prob] < 0.05)
    % Find the first significant positive cluster
    sig_cluster_idx = find([stat.posclusters.prob] < 0.05, 1, 'first');
    if ~isempty(sig_cluster_idx) && isfield(stat, 'posclusterslabelmat') && ~isempty(stat.posclusterslabelmat)
        cluster_mask = stat.posclusterslabelmat == sig_cluster_idx;
        if any(cluster_mask(:))
            % Select channels and time points
            selected_channels = stat.label(any(cluster_mask, 2));
            selected_times = stat.time(any(cluster_mask, 1));
            % Verify channels exist in grand averages
            valid_channels = intersect(selected_channels, diff_grandavg.label);
            if ~isempty(valid_channels)
                cfg = [];
                cfg.layout = struct();
                cfg.layout.label = lay.label(keep_idx);
                cfg.layout.pos = lay.pos(keep_idx, :);
                cfg.layout.width = lay.width(keep_idx);
                cfg.layout.height = lay.height(keep_idx);
                cfg.xlim = [min(selected_times), max(selected_times)]; % Time window of the cluster
                cfg.zlim = [-2 2]; % Amplitude range for difference
                cfg.parameter = 'avg'; % Plot the average ERP
                cfg.comment = sprintf(' ', sig_cluster_idx, stat.posclusters(sig_cluster_idx).prob);
                cfg.marker = 'labels'; % Show channel labels
                cfg.colorbar = 'yes'; % Show colorbar
                figure('Name', sprintf('Topoplot Positive Cluster %d', sig_cluster_idx), 'Color', 'w', 'Position', [100, 100, 600, 600]);
                ft_topoplotER(cfg, diff_grandavg);
                % Post-process to color significant channel labels red and increase font size
                ax = gca;
                text_objects = findobj(ax, 'Type', 'Text');
                for i = 1:length(text_objects)
                    if ismember(lower(strtrim(text_objects(i).String)), lower(valid_channels))
                        text_objects(i).Color = [1 0 0]; % Set label color to red for significant channels
                    end
                    text_objects(i).FontSize = 10; % Increase font size for all channel labels
                end
                % Add title to the colorbar
                cb = colorbar; % Get handle to the colorbar
                title(cb, 'Amplitude (µV)', 'FontSize', 10); % Add title to colorbar
                title(sprintf('Topoplot of Positive Cluster %d (p = %.3f)', sig_cluster_idx, stat.posclusters(sig_cluster_idx).prob), 'FontSize', 12);
                saveas(gcf, fullfile(output_dir, sprintf('topoplot_pos_cluster_%d_0to4000ms.eps', sig_cluster_idx)), 'epsc');
            else
                fprintf('No valid channels found for positive cluster %d in diff_grandavg.\n', sig_cluster_idx);
            end
        else
            fprintf('No valid cluster mask for positive cluster %d.\n', sig_cluster_idx);
        end
    else
        fprintf('No significant positive clusters found with p < 0.05 or missing posclusterslabelmat.\n');
    end
else
    fprintf('No positive clusters detected for topographic plotting.\n');
end

% Plot topographic map for significant negative cluster
if isfield(stat, 'negclusters') && ~isempty(stat.negclusters) && any([stat.negclusters.prob] < 0.05)
    % Find the first significant negative cluster
    sig_cluster_idx = find([stat.negclusters.prob] < 0.05, 1, 'first');
    if ~isempty(sig_cluster_idx) && isfield(stat, 'negclusterslabelmat') && ~isempty(stat.negclusterslabelmat)
        cluster_mask = stat.negclusterslabelmat == sig_cluster_idx;
        if any(cluster_mask(:))
            % Select channels and time points
            selected_channels = stat.label(any(cluster_mask, 2));
            selected_times = stat.time(any(cluster_mask, 1));
            % Verify channels exist in grand averages
            valid_channels = intersect(selected_channels, diff_grandavg.label);
            if ~isempty(valid_channels)
                cfg = [];
                cfg.layout = struct();
                cfg.layout.label = lay.label(keep_idx);
                cfg.layout.pos = lay.pos(keep_idx, :);
                cfg.layout.width = lay.width(keep_idx);
                cfg.layout.height = lay.height(keep_idx);
                cfg.xlim = [min(selected_times), max(selected_times)]; % Time window of the cluster
                cfg.zlim = [-2 2]; % Amplitude range for difference
                cfg.parameter = 'avg'; % Plot the average ERP
                cfg.comment = sprintf('Negative Cluster %d (p = %.3f)', sig_cluster_idx, stat.negclusters(sig_cluster_idx).prob);
                cfg.marker = 'labels'; % Show channel labels
                cfg.highlight = 'on'; % Highlight significant channels
                cfg.highlightchannel = valid_channels; % Channels in the cluster
                cfg.highlightsymbol = '*'; % Symbol for significant channels
                cfg.highlightcolor = [1 0 0]; % Red for significant channels
                cfg.colorbar = 'yes'; % Show colorbar
                figure('Name', sprintf('Topoplot Negative Cluster %d', sig_cluster_idx), 'Color', 'w');
                ft_topoplotER(cfg, diff_grandavg);
                title(sprintf('Topoplot of Negative Cluster %d (p = %.3f)', sig_cluster_idx, stat.negclusters(sig_cluster_idx).prob), 'FontSize', 12);
                saveas(gcf, fullfile(output_dir, sprintf('topoplot_neg_cluster_%d_0to4000ms.eps', sig_cluster_idx)), 'epsc');
            else
                fprintf('No valid channels found for negative cluster %d in diff_grandavg.\n', sig_cluster_idx);
            end
        else
            fprintf('No valid cluster mask for negative cluster %d.\n', sig_cluster_idx);
        end
    else
        fprintf('No significant negative clusters found with p < 0.05 or missing negclusterslabelmat.\n');
    end
else
    fprintf('No negative clusters detected for topographic plotting.\n');
end