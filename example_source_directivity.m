%% Example: Source Directivity with Acoustic Reciprocity Test
% This example demonstrates the source directivity feature of the RIR generator
% and validates it using acoustic reciprocity.
%
% Acoustic reciprocity principle: The impulse response from source A to receiver B
% should be equal to the impulse response from source B to receiver A, when both
% have the same directivity patterns and orientations.

clear; clc;

%% Basic room acoustic parameters
c = 340;                    % Sound velocity (m/s)
fs = 16000;                 % Sample frequency (samples/s)
L = [5 4 6];                % Room dimensions [x y z] (m)
beta = 0.3;                 % Reverberation time (s)
n = 4096;                   % Number of samples
order = 3;                  % Reflection order (limited for faster computation)
dim = 3;                    % Room dimension
hp_filter = 1;              % Enable high-pass filter

%% Position setup
% Position A and B for reciprocity test
pos_A = [1.5 1.5 1.5];      % Position A [x y z] (m)
pos_B = [3.5 2.5 1.5];      % Position B [x y z] (m)

%% Directivity patterns to test
directivity_types = {'omnidirectional', 'cardioid', 'hypercardioid', 'bidirectional'};

% Test orientations (azimuth, elevation) in radians
orientations = [0 0; pi/4 0; pi/2 0; 0 pi/6];  % 0°, 45°, 90°, and 30° elevation

fprintf('Testing Source Directivity with Acoustic Reciprocity\n');
fprintf('====================================================\n\n');

%% Test each directivity pattern
for d_idx = 1:length(directivity_types)
    mtype = directivity_types{d_idx};
    stype = directivity_types{d_idx};  % Same directivity for source and receiver
    
    fprintf('Testing %s directivity pattern:\n', mtype);
    fprintf('--------------------------------------\n');
    
    for o_idx = 1:size(orientations, 1)
        orientation = orientations(o_idx, :);
        
        % Convert to degrees for display
        azimuth_deg = orientation(1) * 180/pi;
        elevation_deg = orientation(2) * 180/pi;
        
        fprintf('Orientation: Azimuth=%.1f°, Elevation=%.1f°\n', azimuth_deg, elevation_deg);
        
        %% Forward direction: Source at A, Receiver at B
        r_forward = pos_B;          % Receiver position
        s_forward = pos_A;          % Source position
        mic_orientation_forward = orientation;    % Microphone orientation
        source_orientation_forward = orientation; % Source orientation
        
        h_forward = rir_generator(c, fs, r_forward, s_forward, L, beta, n, ...
                                 mtype, order, dim, mic_orientation_forward, hp_filter, ...
                                 stype, source_orientation_forward);
        
        %% Reverse direction: Source at B, Receiver at A (reciprocity test)
        r_reverse = pos_A;          % Receiver position (swapped)
        s_reverse = pos_B;          % Source position (swapped)
        mic_orientation_reverse = orientation;    % Microphone orientation
        source_orientation_reverse = orientation; % Source orientation
        
        h_reverse = rir_generator(c, fs, r_reverse, s_reverse, L, beta, n, ...
                                 mtype, order, dim, mic_orientation_reverse, hp_filter, ...
                                 stype, source_orientation_reverse);
        
        %% Compare impulse responses (reciprocity check)
        % Normalize for comparison
        h_forward_norm = h_forward / max(abs(h_forward));
        h_reverse_norm = h_reverse / max(abs(h_reverse));
        
        % Calculate correlation and error metrics
        correlation = max(xcorr(h_forward_norm, h_reverse_norm, 'normalized'));
        mse = mean((h_forward_norm - h_reverse_norm).^2);
        max_diff = max(abs(h_forward_norm - h_reverse_norm));
        
        fprintf('  Correlation: %.4f, MSE: %.6f, Max Diff: %.4f\n', ...
                correlation, mse, max_diff);
        
        %% Store results for plotting
        if d_idx == 2 && o_idx == 2  % Store cardioid with 45° orientation for detailed analysis
            h_forward_plot = h_forward;
            h_reverse_plot = h_reverse;
            time_axis = (0:length(h_forward)-1) / fs * 1000; % time in ms
        end
    end
    fprintf('\n');
end

%% Visualization
figure('Position', [100 100 1200 400]);

% Plot 1: Reciprocity comparison
subplot(1,3,1);
plot(time_axis, h_forward_plot, 'b-', 'LineWidth', 1.5);
hold on;
plot(time_axis, h_reverse_plot, 'r--', 'LineWidth', 1.5);
xlabel('Time (ms)');
ylabel('Amplitude');
title('Reciprocity Test: Cardioid 45°');
legend('Forward A→B', 'Reverse B→A', 'Location', 'northeast');
grid on;
xlim([0 50]);  % Show first 50ms

% Plot 2: Difference between forward and reverse
subplot(1,3,2);
plot(time_axis, h_forward_plot - h_reverse_plot, 'g-', 'LineWidth', 1);
xlabel('Time (ms)');
ylabel('Difference');
title('Reciprocity Error');
grid on;
xlim([0 50]);

% Plot 3: Frequency response comparison
subplot(1,3,3);
[H_forward, f] = freqz(h_forward_plot, 1, 1024, fs);
[H_reverse, ~] = freqz(h_reverse_plot, 1, 1024, fs);
semilogx(f, 20*log10(abs(H_forward)), 'b-', 'LineWidth', 1.5);
hold on;
semilogx(f, 20*log10(abs(H_reverse)), 'r--', 'LineWidth', 1.5);
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('Frequency Response Comparison');
legend('Forward', 'Reverse', 'Location', 'southwest');
grid on;
xlim([100 8000]);
