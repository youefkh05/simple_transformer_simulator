% =============================================
% 75-kVA Distribution Transformer Analysis
% With Efficiency and Regulation Plots
% =============================================

clear all;
close all;
clc;

%% Transformer Parameters
S_rated = 75e3;         % Rated power (VA)
f = 60;                 % Frequency (Hz)
V1_rated = 4600;        % HV side rated voltage (V)
V2_rated = 240;         % LV side rated voltage (V)
a = V1_rated/V2_rated;  % Turns ratio

% Given parameters (in ohms)
R1 = 0.846;             % HV winding resistance
R2 = 0.00261;           % LV winding resistance
X1 = 26.8;              % HV leakage reactance
X2 = 0.0745;            % LV leakage reactance
Rc = 220000;            % Core loss resistance
Xm = 112000;            % Magnetizing reactance

% Base values
I2_rated = S_rated / V2_rated;  % Rated LV current
Z2_base = V2_rated / I2_rated;  % LV base impedance

%% Part a: HV Terminal Voltage vs Power Factor
pf_angles = acos([0.6 0.7 0.8 0.9 1.0 0.9 0.8 0.7 0.6]); % pf angles
pf_labels = {'0.6 lead', '0.7 lead', '0.8 lead', '0.9 lead', '1.0', ...
             '0.9 lag', '0.8 lag', '0.7 lag', '0.6 lag'};
theta_range = -pf_angles(1:4) pf_angles(5) pf_angles(6:9); % Angle range

V1_terminal = zeros(size(pf_angles));
regulation = zeros(size(pf_angles));

for i = 1:length(pf_angles)
    theta = theta_range(i);
    I2 = I2_rated * (cos(theta) + 1i*sin(theta)); % Full load current
    
    % Refer LV side quantities to HV side
    I2_prime = I2 / a;
    Zeq = (R1 + a^2*R2) + 1i*(X1 + a^2*X2); % Equivalent impedance
    
    % Calculate terminal voltage
    E1 = V2_rated * a + I2_prime * Zeq;
    V1_terminal(i) = abs(E1);
    
    % Calculate regulation
    regulation(i) = (V1_terminal(i) - V2_rated*a) / (V2_rated*a) * 100;
end

% Plot HV terminal voltage vs power factor
figure('Name', 'HV Terminal Voltage', 'NumberTitle', 'off');
plot(rad2deg(theta_range), V1_terminal, 'b-o', 'LineWidth', 2);
xlabel('Power Factor Angle (degrees)');
ylabel('HV Terminal Voltage (V)');
title('HV Terminal Voltage vs Power Factor at Full Load');
grid on;
xticks(rad2deg(theta_range));
xticklabels(pf_labels);
xlim([rad2deg(theta_range(1)) rad2deg(theta_range(end))]);

%% Part b: Efficiency and Regulation Plots
load_factors = [1.0, 0.5, 0.25]; % Full, half, and quarter load
colors = ['r', 'g', 'b'];
markers = ['o', 's', 'd'];

figure('Name', 'Efficiency and Regulation', 'NumberTitle', 'off', 'Position', [100 100 1200 500]);

% Efficiency subplot
subplot(1,2,1);
hold on;
for k = 1:length(load_factors)
    K = load_factors(k);
    efficiency = zeros(size(pf_angles));
    
    for i = 1:length(pf_angles)
        theta = theta_range(i);
        I2 = K * I2_rated * (cos(theta) + 1i*sin(theta));
        
        % Refer LV side quantities to HV side
        I2_prime = I2 / a;
        Zeq = (R1 + a^2*R2) + 1i*(X1 + a^2*X2);
        
        % Calculate voltages and currents
        E1 = V2_rated * a + I2_prime * Zeq;
        I1 = I2_prime + E1/Rc + E1/(1i*Xm);
        
        % Power calculations
        P_out = real(V2_rated * conj(I2));
        P_cu = abs(I2)^2 * (a^2*R2 + R1); % Copper losses
        P_core = abs(E1)^2 / Rc;           % Core losses
        P_in = P_out + P_cu + P_core;
        
        efficiency(i) = P_out / P_in * 100;
    end
    
    plot(rad2deg(theta_range), efficiency, [colors(k) '-' markers(k)], 'LineWidth', 1.5, ...
        'DisplayName', sprintf('%d%% Load', K*100));
end

xlabel('Power Factor Angle (degrees)');
ylabel('Efficiency (%)');
title('Efficiency vs Power Factor');
grid on;
xticks(rad2deg(theta_range));
xticklabels(pf_labels);
xlim([rad2deg(theta_range(1)) rad2deg(theta_range(end))]);
legend('Location', 'best');
ylim([95 100]); % Typical efficiency range for distribution transformers

% Regulation subplot
subplot(1,2,2);
hold on;
for k = 1:length(load_factors)
    K = load_factors(k);
    regulation = zeros(size(pf_angles));
    
    for i = 1:length(pf_angles)
        theta = theta_range(i);
        I2 = K * I2_rated * (cos(theta) + 1i*sin(theta));
        
        % Refer LV side quantities to HV side
        I2_prime = I2 / a;
        Zeq = (R1 + a^2*R2) + 1i*(X1 + a^2*X2);
        
        % Calculate terminal voltage
        E1 = V2_rated * a + I2_prime * Zeq;
        regulation(i) = (abs(E1) - V2_rated*a) / (V2_rated*a) * 100;
    end
    
    plot(rad2deg(theta_range), regulation, [colors(k) '-' markers(k)], 'LineWidth', 1.5, ...
        'DisplayName', sprintf('%d%% Load', K*100));
end

xlabel('Power Factor Angle (degrees)');
ylabel('Voltage Regulation (%)');
title('Voltage Regulation vs Power Factor');
grid on;
xticks(rad2deg(theta_range));
xticklabels(pf_labels);
xlim([rad2deg(theta_range(1)) rad2deg(theta_range(end))]);
legend('Location', 'best');

%% Tapping Range Analysis
max_V1 = max(V1_terminal);
min_V1 = min(V1_terminal);
nominal_V1 = V2_rated * a;

tapping_range = 5; % ±5%
min_allowed = nominal_V1 * (1 - tapping_range/100);
max_allowed = nominal_V1 * (1 + tapping_range/100);

fprintf('\n==== Tapping Range Analysis ====\n');
fprintf('Nominal HV voltage: %.1f V\n', nominal_V1);
fprintf('Maximum required HV voltage: %.1f V\n', max_V1);
fprintf('Minimum required HV voltage: %.1f V\n', min_V1);
fprintf('Allowed tapping range (±5%%): %.1f V to %.1f V\n', min_allowed, max_allowed);

if max_V1 > max_allowed || min_V1 < min_allowed
    fprintf('\nWARNING: The ±5%% tapping range is INSUFFICIENT to maintain\n');
    fprintf('the required secondary voltage under all load conditions.\n');
    
    % Calculate required tapping range
    required_min_tap = (min_V1 - nominal_V1)/nominal_V1 * 100;
    required_max_tap = (max_V1 - nominal_V1)/nominal_V1 * 100;
    
    fprintf('\nRequired tapping range: %.1f%% to %.1f%%\n', ...
        required_min_tap, required_max_tap);
else
    fprintf('\nThe ±5%% tapping range is SUFFICIENT to maintain\n');
    fprintf('the required secondary voltage under all load conditions.\n');
end