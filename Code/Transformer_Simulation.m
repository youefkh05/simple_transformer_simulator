% Ideal Transformer Simulation
clear all;
close all;
clc;

% Parameters
f = 50;                 % Frequency (Hz)
Vp = 230;               % Primary voltage (V)
Np = 500;               % Primary turns
Ns = 100;               % Secondary turns
a = Np/Ns;              % Turns ratio
RL = 10;                % Load resistance (Ohm)

% Time vector
t = linspace(0, 0.1, 1000); % 10 cycles

% Primary voltage
Vp_t = Vp * sqrt(2) * sin(2*pi*f*t);

% Secondary voltage (ideal)
Vs_t = Vp_t / a;

% Secondary current
Is_t = Vs_t / RL;

% Primary current (ideal)
Ip_t = -Is_t / a;

% Plot results
figure;
subplot(2,1,1);
plot(t, Vp_t, 'b', t, Vs_t, 'r');
title('Primary and Secondary Voltages');
xlabel('Time (s)');
ylabel('Voltage (V)');
legend('Primary V_p', 'Secondary V_s');
grid on;

subplot(2,1,2);
plot(t, Ip_t, 'b', t, Is_t, 'r');
title('Primary and Secondary Currents');
xlabel('Time (s)');
ylabel('Current (A)');
legend('Primary I_p', 'Secondary I_s');
grid on;