%% Problem 1 
clc; clear; close all;
R = 10e3;                 
C = 133e-9;               
f_drop = 1/(4*pi*R*C);    
fprintf('Drop (notch) frequency: %.2f Hz\n', f_drop);
% Part (a): 
f = 0:0.1:200;            
w = 2*pi*f;              
m_list = [0.8, 0.9];
Hmag = zeros(numel(m_list), numel(f));
Hphs = zeros(numel(m_list), numel(f));
for k = 1:numel(m_list)
   m = m_list(k);
   sRC2 = (2*1j*w*R*C).^2;                 
   sRC  =  2*1j*w*R*C;                     
   H = (1+m).*(sRC2 + 1) ./ (sRC2 + 4*(1-m).*sRC + 1);
   Hmag(k,:) = abs(H);
   Hphs(k,:) = unwrap(angle(H));
end
figure('Name','Part (a): H(w) magnitude and phase');
subplot(2,1,1);
plot(f, Hmag(1,:), 'LineWidth',1.4); hold on;
plot(f, Hmag(2,:), 'LineWidth',1.4);
grid on; xlabel('f [Hz]'); ylabel('|H(j\omega)|');
title('Magnitude of H(j\omega)');
legend('m = 0.8', 'm = 0.9', 'Location','best');
subplot(2,1,2);
plot(f, Hphs(1,:), 'LineWidth',1.4); hold on;
plot(f, Hphs(2,:), 'LineWidth',1.4);
grid on; xlabel('f [Hz]'); ylabel('angle H(j\omega) [rad]');
title('Phase of H(j\omega)');
legend('m = 0.8', 'm = 0.9', 'Location','best');

% Part (b): 
m = 0.9;
S = load('ecg_signal.mat');        
if isfield(S,'ecg')
   x = double(S.ecg(:));           
else
   error('ecg_signal.mat must contain variable "ecg".');
end
if isfield(S,'Fs')
   Fs = S.Fs;
elseif isfield(S,'t')
   t_ecg = S.t(:);
   Fs = 1/mean(diff(t_ecg));
else
   Fs = 500;                      
end
Tmax = 2.5;
N_des = round(Tmax*Fs);
if length(x) >= N_des
   x = x(1:N_des);
else
   x = [x; zeros(N_des - length(x),1)];
end
N = length(x);
t = (0:N-1)/Fs;                     
X = fftshift(fft(x));
fshift = (-N/2:N/2-1)*(Fs/N);       
wshift = 2*pi*fshift;               
sRC2 = (2*1j*wshift*R*C).^2;
sRC  =  2*1j*wshift*R*C;
H_w  = (1+m).*(sRC2 + 1) ./ (sRC2 + 4*(1-m).*sRC + 1);
Z = H_w .* X;
z = ifft(ifftshift(Z), 'symmetric');   
% Plots:
idx = (fshift >= -250) & (fshift <= 250);
figure('Name','Part (b): ECG filtering with twin-T notch (m=0.9)');
subplot(4,1,1);
plot(t, x, 'LineWidth',1.1); grid on;
xlim([0 2.5]); xlabel('t [s]'); ylabel('x(t)');
title('Input ECG x(t)');
subplot(4,1,2);
plot(fshift(idx), abs(X(idx)), 'LineWidth',1.1); grid on;
xlim([-250 250]); xlabel('f [Hz]'); ylabel('|X(f)|');
title('Spectrum of input ECG');
subplot(4,1,3);
plot(fshift(idx), abs(Z(idx)), 'LineWidth',1.1); grid on;
xlim([-250 250]); xlabel('f [Hz]'); ylabel('|Z(f)|');
title('Spectrum after twin-T notch (m=0.9)');
subplot(4,1,4);
plot(t, z, 'LineWidth',1.1); grid on;
xlim([0 2.5]); xlabel('t [s]'); ylabel('z(t)');
title('Filtered output z(t)');


%% Problem 2:
if ~exist('Fs','var') || ~exist('x','var') || ~exist('z','var') || ~exist('t','var')
   error('Run Part (b) first so x, z, Fs, t exist in the workspace.');
end
dt = 1/Fs;                
N  = numel(x);
E_x_time = trapz(t, abs(x).^2);
E_z_time = trapz(t, abs(z).^2);

X = fft(x);
Z = fft(z);
E_x_freq = (dt/N) * sum(abs(X).^2);
E_z_freq = (dt/N) * sum(abs(Z).^2);
fprintf('Energy x(t): time  = %.6g,   freq = %.6g,  |diff| = %.3e\n', ...
       E_x_time, E_x_freq, abs(E_x_time - E_x_freq));
fprintf('Energy z(t): time  = %.6g,   freq = %.6g,  |diff| = %.3e\n', ...
       E_z_time, E_z_freq, abs(E_z_time - E_z_freq));

Xc = fftshift(X);  Zc = fftshift(Z);
f  = (-N/2:N/2-1)*(Fs/N);          
Sx = (dt/N) * abs(Xc).^2;
Sz = (dt/N) * abs(Zc).^2;
band = (f >= -250) & (f <= 250);
figure('Name','Problem 2: Energy spectra');
subplot(2,1,1);
plot(f(band), Sx(band), 'LineWidth', 1.2); grid on;
xlabel('f [Hz]'); ylabel('Energy density');
title('Energy spectrum of X(f)');
subplot(2,1,2);
plot(f(band), Sz(band), 'LineWidth', 1.2); grid on;
xlabel('f [Hz]'); ylabel('Energy density');
title('Energy spectrum of Z(f)');

%% Problem 3:
clc; clear; close all;
% Parameters:
m_syms = [6 0 4 -6 2];     
fc     = 500e3;            
phi    = pi/3;            
fsym   = 50e3;             
Fs     = 10e6;             
reps   = 200;              
Ns   = round(Fs/fsym);                   
seq  = repmat(m_syms, 1, reps);           
m_t  = repelem(seq, Ns).';                
N    = numel(m_t);
t    = (0:N-1).' / Fs;                    
% Modulation:
s_t  = m_t .* cos(2*pi*fc*t);             
% Spectra:
S    = fftshift(fft(s_t));
f    = (-N/2:N/2-1).' * (Fs/N);           
% Demodulation:
lo_t = cos(2*pi*fc*t + phi);              
v_t  = s_t .* lo_t;                       
V    = fftshift(fft(v_t));
Hf   = 2.0 * (abs(f) < 500e3);
Vo   = Hf .* V;                           
vo_t = ifft(ifftshift(Vo), 'symmetric');  
theory_scale = cos(phi);   
nz = abs(m_t) > 0;
meas_scale = median(vo_t(nz) ./ m_t(nz)); 
fprintf('Expected scale = %.3f, measured ~ %.3f\n', theory_scale, meas_scale);
% Plots (4x1):
t_show = t <= 5e-4;   
figure('Name','DSB-SC modulation/demodulation');
subplot(4,1,1);
plot(t(t_show)*1e3, s_t(t_show), 'LineWidth', 1.1); grid on;
xlabel('t [ms]'); ylabel('s(t)'); title('Transmitted DSB-SC s(t)');
subplot(4,1,2);
plot(f/1e3, abs(S), 'LineWidth', 1.1); grid on; xlim([-1.5 1.5]*1e3);
xlabel('f [kHz]'); ylabel('|S(f)|'); title('Magnitude spectrum |S(f)|');
subplot(4,1,3);
plot(t(t_show)*1e3, vo_t(t_show), 'LineWidth', 1.1); grid on;
xlabel('t [ms]'); ylabel('v_o(t)'); title('Demodulated and LPF output v_o(t)');
subplot(4,1,4);
Vo_mag = abs(fftshift(fft(vo_t)));
plot(f/1e3, Vo_mag, 'LineWidth', 1.1); grid on; xlim([-200 200]); 
xlabel('f [kHz]'); ylabel('|V_o(f)|'); title('Magnitude spectrum |V_o(f)|');
% Test case
theory = cos(pi/3);                         
nz = abs(m_t)>0;                            
meas = median(vo_t(nz)./m_t(nz));
fprintf('Expected scale: %.3f, measured: %.3f\n', theory, meas);
