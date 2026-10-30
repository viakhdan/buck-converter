% Buck converter - PID controller design
 
clear; clc; close all;
 
% Power stage parameters
Vin  = 12;          % nominal input voltage [V]
Vout = 5;           % output voltage [V]
L    = 33e-6;       % inductance [H]
C    = 47e-6;       % output capacitance [F]
ESR  = 10e-3;       % capacitor ESR [Ohm]
R    = 1.25;        % nominal (full) load [Ohm]
fsw  = 100e3;       % switching frequency [Hz]
D    = Vout/Vin;    % nominal duty cycle
 
% Loop components
s   = tf('s');
Fm  = 1;                           % PWM modulator gain (sawtooth 0..1)
f_filt = 20e3;                     % measurement filter corner [Hz]
H   = 1/(s/(2*pi*f_filt) + 1);     % measurement low-pass filter
 
% plant transfer function (averaged CCM model)
Gvd = @(Vin_x, R_x) Vin_x*(1 + s*C*ESR) / ...
      (1 + s*(L/R_x + C*ESR) + s^2*L*C*(R_x + ESR)/R_x);
 
G = Gvd(Vin, R);       % nominal plant
P = Fm*G*H;            % plant seen by the controller
 
% plant parameters
f0   = 1/(2*pi*sqrt(L*C));          % LC resonance [Hz]
zeta = 1/(2*R)*sqrt(L/C);           % damping at nominal load
fz_esr = 1/(2*pi*ESR*C);            % ESR zero [Hz]
fprintf('Plant: DC gain = %.1f (%.1f dB), f0 = %.2f kHz, zeta = %.3f, ESR zero = %.0f kHz\n', Vin, 20*log10(Vin), f0/1e3, zeta, fz_esr/1e3);
 
% Controller design
fc  = 5e3;                 % target crossover frequency [Hz]
wc  = 2*pi*fc;
wz  = 2*pi*0.5*f0;         % real double zero at half the LC resonance
wp  = 2*pi*50e3;           % high-frequency pole (derivative filter)
 
Cbase = (1 + s/wz)^2 / (s*(1 + s/wp));
K     = 1/abs(squeeze(freqresp(Cbase*P, wc)));   % gain for crossover at fc
Cman  = K*Cbase;
 
Tf = 1/wp;
Ki = K;
Kp = 2*K/wz - K*Tf;
Kd = K/wz^2 - Kp*Tf;
N  = 1/Tf;                 % filter coefficient
Cpid = pid(Kp, Ki, Kd, Tf);
 
fprintf('\nPID gains (parallel form, for Simulink PID block):\n');
fprintf('  Kp = %.5g\n  Ki = %.5g\n  Kd = %.5g\n  N  = %.5g\n', Kp, Ki, Kd, N);
 
% Nominal analysis
figure('Name','Loop gain (nominal)');
margin(Cpid*P); grid on;
 
figure('Name','Closed-loop step (nominal)');
step(feedback(Cpid*G, H)); grid on;
title('Reference-to-output step response (nominal, linear model)');
 
% Corner check
Vin_list = [10 12 14];
R_list   = [1.25 5 10];
colors   = {'b', 'k', 'r'};        % one color per input voltage
styles   = {'-', '--', ':'};       % one line style per load
 
f_plot = logspace(1, log10(200e3), 2000);   % 10 Hz ... 200 kHz
w_plot = 2*pi*f_plot;
 
figure('Name','Loop gain - all corners');
axMag = subplot(2,1,1); hold(axMag, 'on'); grid(axMag, 'on');
axPh  = subplot(2,1,2); hold(axPh,  'on'); grid(axPh,  'on');
 
labels = {};
fprintf('\nCorner check:\n');
fprintf('  Vin [V]   R [Ohm]   PM [deg]   fc [kHz]   GM [dB]\n');
for i = 1:numel(Vin_list)
    for j = 1:numel(R_list)
        Lt = Cpid*Fm*Gvd(Vin_list(i), R_list(j))*H;
        [Gm, Pm, ~, Wcp] = margin(Lt);
        fprintf('  %5d     %6.2f    %7.1f    %7.2f    %6.1f\n', ...
                Vin_list(i), R_list(j), Pm, Wcp/(2*pi)/1e3, 20*log10(Gm));
 
        % Frequency response data, plotted manually so the legend is reliable
        [mag, ph] = bode(Lt, w_plot);
        mag = squeeze(mag);  ph = squeeze(ph);
        spec = [colors{i} styles{j}];
        semilogx(axMag, f_plot, 20*log10(mag), spec, 'LineWidth', 1.2);
        semilogx(axPh,  f_plot, ph,            spec, 'LineWidth', 1.2);
 
        labels{end+1} = sprintf('%d V, %.2f \\Omega', Vin_list(i), R_list(j));
    end
end
 
set([axMag axPh], 'XScale', 'log', 'XLim', [f_plot(1) f_plot(end)]);
yline(axMag, 0, 'k-', 'HandleVisibility', 'off');      % 0 dB line
yline(axPh, -180, 'k-', 'HandleVisibility', 'off');    % -180 deg line
ylabel(axMag, 'Magnitude (dB)');
ylabel(axPh,  'Phase (deg)');
xlabel(axPh,  'Frequency (Hz)');
title(axMag, 'Loop gain at all operating corners');
legend(axMag, labels, 'Location', 'southwest', 'NumColumns', 3);