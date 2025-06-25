%% Example: Dipole Directivity with Acoustic Reciprocity Test
% This example focuses on bidirectional (dipole) directivity patterns
% and validates it using acoustic reciprocity.

clear; clc;

%% Room and simulation parameters
c = 340;                    % Sound velocity (m/s)
fs = 16000;                 % Sample frequency (samples/s)
L = [4 3 2.5];              % Room dimensions [x y z] (m)
beta = 0.25;                % Reverberation time (s)
n = 2048;                   % Number of samples
order = 2;                  % Reflection order
dim = 3;                    % Room dimension
hp_filter = 1;              % Enable high-pass filter

fprintf('Dipole (Bidirectional) Directivity\n');
fprintf('==================================\n\n');

%% Test 1: Acoustic Reciprocity with Dipoles
pos_A = [1 1.25 1.2];        % Position A
pos_B = [3 1.75 1.2];        % Position B (same z, different x and y)

fprintf('Test 1: Acoustic Reciprocity\n');
fprintf('----------------------------\n');

% Test different orientations
test_angles = [0, pi/4, pi/2, 3*pi/4];  % 0°, 45°, 90°, 135°
angle_names = {'0°', '45°', '90°', '135°'};

fprintf('Bidirectional (dipole) pattern reciprocity test:\n');
correlations = zeros(size(test_angles));

for i = 1:length(test_angles)
    angle = test_angles(i);
    orientation_A = [0, 0];  % Fixed
    orientation_B = [angle, 0];  % Only azimuth rotation
    
    % Forward: dipole source at A, dipole receiver at B
    h_AB = rir_generator(c, fs, pos_B, pos_A, L, beta, n, ...
                        'bidirectional', order, dim, orientation_B, hp_filter, ...
                        'bidirectional', orientation_A);
    
    % Reverse: dipole source at B, dipole receiver at A
    h_BA = rir_generator(c, fs, pos_A, pos_B, L, beta, n, ...
                        'bidirectional', order, dim, orientation_A, hp_filter, ...
                        'bidirectional', orientation_B);
    
    % Calculate correlation
    correlation = max(xcorr(h_AB/max(abs(h_AB)), h_BA/max(abs(h_BA)), 'normalized'));
    correlations(i) = correlation;
    
    fprintf('  Orientation %s: Correlation = %.4f\n', angle_names{i}, correlation);
end

%% Test 2: Bidirectional Source Orientation Effect
fprintf('\nTest 2: Bidirectional Source Orientation Effect\n');
fprintf('----------------------------------------------\n');

% Fixed positions for directivity comparison
r_test = [2.5 1.5 1.2];     % Receiver position
s_test = [1.5 1.5 1.2];     % Source position

order = 0; % Direct-path only

% Test bidirectional source at different orientations
dipole_angles = linspace(0, 2*pi, 90);
dipole_angles_deg = dipole_angles * 180/pi;
dipole_energies = zeros(size(dipole_angles));

fprintf('Bidirectional source at different orientations:\n');
for i = 1:length(dipole_angles)
    orientation = [dipole_angles(i), 0];
    
    h = rir_generator(c, fs, r_test, s_test, L, beta, n, ...
                     'omnidirectional', order, dim, [0 0], hp_filter, ...
                     'bidirectional', orientation);
    
    energy = sqrt(mean(h.^2));
    dipole_energies(i) = energy;
    fprintf('  %.0f°: Energy = %.4f\n', dipole_angles_deg(i), energy);
end

%% Generate detailed comparison for visualization
% Fixed positions for directivity comparison
r_test = [2.5 1.5 1.2];     % Receiver position
s_test = [1.5 1.5 1.2];     % Source position

order = -1; % Maximum reflection order

% Get impulse responses for different source types
h_omni = rir_generator(c, fs, r_test, s_test, L, beta, n, ...
                      'omnidirectional', order, dim, [0 0], hp_filter, ...
                      'omnidirectional', [pi/2 0]);

h_dipole_0 = rir_generator(c, fs, r_test, s_test, L, beta, n, ...
                          'omnidirectional', order, dim, [0 0], hp_filter, ...
                          'bidirectional', [0 0]);

h_dipole_90 = rir_generator(c, fs, r_test, s_test, L, beta, n, ...
                           'omnidirectional', order, dim, [0 0], hp_filter, ...
                           'bidirectional', [pi/2 0]);

% Time axis
t = (0:length(h_omni)-1) / fs * 1000;  % milliseconds

%% Visualization
figure('Position', [50 50 1200 900]);

% Plot 1: Bidirectional source orientation effect
subplot(2,2,1);
polar_angles = linspace(0, 2*pi, length(dipole_angles));
polarplot(polar_angles, dipole_energies, 'ro-', 'LineWidth', 2, 'MarkerSize', 6);
title('Dipole Energy vs Orientation');

% Plot 2: Theoretical bidirectional pattern
subplot(2,2,2);
theta_theory = linspace(0, 2*pi, 90);
dipole_pattern = abs(cos(theta_theory));  % |cos(θ)| pattern for bidirectional
polarplot(theta_theory, dipole_pattern, 'b--', 'LineWidth', 2);
title('Theoretical Dipole Pattern');

% Plot 3: Impulse responses comparison
subplot(2,2,3);
plot(t, h_omni, 'k-', 'LineWidth', 1.5);
hold on;
plot(t, h_dipole_0, 'r--', 'LineWidth', 1.5);
plot(t, h_dipole_90, 'b--', 'LineWidth', 1.5);
xlabel('Time (ms)');
ylabel('Amplitude');
title('Impulse Response Comparison');
legend('Omnidirectional', 'Dipole 0°', 'Dipole 90°', 'Location', 'northeast');
grid on;
xlim([0 100]);

% Plot 4: Frequency response comparison
subplot(2,2,4);
[H_omni, f] = freqz(h_omni, 1, 512, fs);
[H_dipole_0, ~] = freqz(h_dipole_0, 1, 512, fs);
[H_dipole_90, ~] = freqz(h_dipole_90, 1, 512, fs);
semilogx(f, 20*log10(abs(H_omni)), 'k-', 'LineWidth', 1.5);
hold on;
semilogx(f, 20*log10(abs(H_dipole_0)), 'r--', 'LineWidth', 1.5);
semilogx(f, 20*log10(abs(H_dipole_90)), 'b--', 'LineWidth', 1.5);
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('Frequency Response');
legend('Omnidirectional', 'Dipole 0°', 'Dipole 90°', 'Location', 'southwest');
grid on;
xlim([100 fs/2]);
