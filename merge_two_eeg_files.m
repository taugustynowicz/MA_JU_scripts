clc, clear, close all

% Load the first .mat file
load('D:\marek\chirps_ftrip\workspace_hapt_2025_erp\mats\chirp_hap12_stderp_std\chirp_hap12_stderp_std_p172_01a_preprocessed.mat'); % Assumes the variable is named 'preprocessed'
preprocessed1 = preprocessed;

% Load the second .mat file
load('D:\marek\chirps_ftrip\workspace_hapt_2025_erp\mats\chirp_hap12_stderp_std\chirp_hap12_stderp_std_p172_01b_preprocessed.mat'); % Assumes the variable is named 'preprocessed'
preprocessed2 = preprocessed;

% Verify compatibility of headers and labels
if ~isequal(preprocessed1.hdr.label, preprocessed2.hdr.label) || preprocessed1.fsample ~= preprocessed2.fsample
    error('Headers or sampling rates are not compatible between the two files.');
end

% Get the number of samples in each trial of the first dataset
nsamples1 = cellfun(@length, preprocessed1.time); % Number of samples per trial
total_samples1 = sum(nsamples1); % Total samples across all trials in first dataset

% Concatenate trial data
merged_data.trial = [preprocessed1.trial, preprocessed2.trial];

% Concatenate time data
merged_data.time = [preprocessed1.time, preprocessed2.time];

% Concatenate trialinfo
merged_data.trialinfo = [preprocessed1.trialinfo; preprocessed2.trialinfo];

% Adjust sampleinfo for the second dataset
sampleinfo2_adjusted = preprocessed2.sampleinfo + total_samples1; % Shift sample indices
merged_data.sampleinfo = [preprocessed1.sampleinfo; sampleinfo2_adjusted];

% Copy other fields from the first dataset (assuming they are the same)
merged_data.hdr = preprocessed1.hdr;
merged_data.label = preprocessed1.label;
merged_data.fsample = preprocessed1.fsample;

% Merge cfg, keeping the most relevant configuration
merged_data.cfg = preprocessed1.cfg; % You may need to customize this if cfgs differ
merged_data.cfg.previous = {preprocessed1.cfg, preprocessed2.cfg}; % Store both cfgs for reference

% Save the merged data
save('chhap12_p172_01_zol26_merged_SG_linear_ERP.mat', 'merged_data');

% Optional: Clear temporary variables
clear preprocessed1 preprocessed2 sampleinfo2_adjusted;

