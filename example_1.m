c = 340;                    % Sound velocity (m/s)
fs = 16000;                 % Sample frequency (samples/s)
r = [-1 1 2];              % Receiver position [x y z] (m)
s = [2 3.5 2];              % Source position [x y z] (m)
L = [5 4 6];                % Room dimensions [x y z] (m)
beta = 0.4;                 % Reverberation time (s)
n = 4096;                   % Number of samples

h = rir_generator(c, fs, r, s, L, beta, n);


plot(h)
sound(h)
xlabel('Time')
ylabel('Audio Signal')