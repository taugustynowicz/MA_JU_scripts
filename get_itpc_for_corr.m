%% Initial cleaning
clc, clear, close all

%% Define paths and load data
% Define the directory path
folderPath = 'D:\marek\chirps_ftrip\workspace_hapt_2025_tfr\mats\patients_clust1';

% Add the folder to MATLAB's search path
addpath(folderPath);

% Specify the file name to load
fileName = 'chirp_hap12_patients_stdtfr_std_cluster1_nobascor_envelope_results.mat';

% Construct the full file path
fullFilePath = fullfile(folderPath, fileName);

% Load the .mat file
load(fullFilePath);

%% Process env_data
% Initialize arrays to store averages
numSubjects = length(env_data); % 23 subjects
avgGroup1 = zeros(numSubjects, 1); % For time points 5-9
avgGroup2 = zeros(numSubjects, 1); % For time points 10-17

% Loop through each subject
for subj = 1:numSubjects
    % Extract data for the current subject
    subjectData = env_data{subj}.data;
    
    % Calculate averages for the specified time point groups
    avgGroup1(subj) = mean(subjectData(5:9)); % Average of indices 5-9
    avgGroup2(subj) = mean(subjectData(10:17)); % Average of indices 10-17
end

%% Create table and save to CSV
% Verify subj_id exists and has the correct number of elements
if ~exist('subj_id', 'var')
    error('Variable subj_id does not exist in the loaded .mat file.');
end

% Assign subj_id to subjectIDs
subjectIDs = subj_id;

% Ensure subjectIDs is a column cell array
if isrow(subjectIDs)
    subjectIDs = subjectIDs'; % Transpose to column (23x1)
end

% Verify size match
if length(subjectIDs) ~= numSubjects
    error('subj_id has %d elements, but env_data has %d subjects. They must match.', ...
          length(subjectIDs), numSubjects);
end

% Create a table with subject IDs and averages
resultsTable = table(subjectIDs, avgGroup1, avgGroup2, ...
    'VariableNames', {'SubjectID', 'Avg_TimePoints_5to9', 'Avg_TimePoints_10to17'});

% Define output CSV file path
outputFilePath = fullfile(folderPath, 'subject_averages.csv');

% Write the table to a CSV file with explicit comma delimiter
writetable(resultsTable, outputFilePath, 'Delimiter', ',');

% Display confirmation
disp(['Results saved to: ' outputFilePath]);