c = 340;                    % Sound velocity (m/s)
fs = 16000;                 % Sample frequency (samples/s)
r = [2 1.5 2];              % Receiver position [x y z] (m)
s = [2 3.5 2];              % Source position [x y z] (m)
L = [5 4 6];                % Room dimensions [x y z] (m)
beta = 0.4;                 % Reflections Coefficients
n = 2048;                   % Number of samples
mtype = 'omnidirectional';  % Type of microphone
order = 2;                  % Reflection order
dim = 3;                    % Room dimension
orientation = 0;            % Microphone orientation (rad)
hp_filter = 1;              % Enable high-pass filter

h = rir_generator(c, fs, r, s, L, beta, n, mtype, order, dim, orientation, hp_filter);