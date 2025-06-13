clear; clc; close all;

% Parameters
duration = 4;           % Total time in seconds
fs = 1000;              % Sampling frequency (Hz)
f_start = 1;            % Starting frequency (Hz)
f_end = 35;             % Ending frequency (Hz)

% Time vector
t = 0:1/fs:duration-1/fs;

% Linear frequency sweep (chirp signal)
k = (f_end - f_start) / duration;  % Frequency sweep rate
phase = 2*pi * (f_start * t + k * t.^2 / 2);
sine_wave = sin(phase);

% Instantaneous frequency for verification
inst_freq = f_start + k * t;

% Fixed frequency sine wave at 25 Hz
fixed_freq = 25;  % Hz
fixed_sine = sin(2*pi * fixed_freq * t);
fixed_freq_array = fixed_freq * ones(size(t));  % Constant frequency array for plotting

% Create the plot
figure('Position', [100, 100, 1200, 1200]);

% Plot the frequency sweep sine wave
subplot(3,1,1);
plot(t, sine_wave, 'k-', 'LineWidth', 1);
xlabel('Time (s)');
title('Standard stimuli (Ascending frequency from 1 to 35 Hz)');
grid on;
xlim([0 duration]);
set(gca, 'YTickLabel', []);

% Plot the fixed frequency sine wave
subplot(3,1,2);
plot(t, fixed_sine, 'k-', 'LineWidth', 1);
xlabel('Time (s)');
title('Target stimuli (Constant frequency of 25 Hz)');
grid on;
xlim([0 duration]);
set(gca, 'YTickLabel', []);

% Plot the instantaneous frequency for both signals
subplot(3,1,3);
plot(t, inst_freq, 'b-', 'LineWidth', 1, 'DisplayName', 'Standard (1 to 35 Hz)');
hold on;
plot(t, fixed_freq_array, 'r-', 'LineWidth', 1, 'DisplayName', 'Target (25 Hz)');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Instantaneous Frequency');
grid on;
xlim([0 duration]);
ylim([0 40]);  % Set y-limit to show frequency range clearly
legend('show');

% Optional: Save the data
% save('frequency_sweep_data.mat', 't', 'sine_wave', 'inst_freq', 'fixed_sine', 'fixed_freq_array');