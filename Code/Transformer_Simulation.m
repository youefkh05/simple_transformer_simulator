clc;
clear;
close all;
% Transformer parameters
S_rated = 75e3;
V1_nom = 4600;
V2 = 240;
R1 = 0.846;
R2 = 0.00261;
X1 = 26.8;
X2 = 0.0745;
Rc = 220e3;
Xmu = 112e3;
% Turns ratio
a = V1_nom / V2;
% Referred secondary impedances to primary
R2_prime = R2 * a^2;
X2_prime = X2 * a^2;
R_total = R1 + R2_prime;
X_total = X1 + X2_prime;
% Full load secondary current
I2_full = S_rated / V2;
% power factor varies from 0.6 pf leading through unity pf to 0.6 pf lagging.
theta_deg = linspace(-54, 54, 100);
theta_values = deg2rad(theta_deg); % radians
pf = cos(theta_values);
% Part a: HV terminal voltage vs power factor at full-load
V1_values = zeros(size(pf));
for k = 1:length(theta_values)
theta = theta_values(k);
I_load_prime = (I2_full / a) * exp(-1j * theta);
V_drop = I_load_prime * (R_total + 1j * X_total);
V1 = V1_nom + V_drop;
V1_values(k) = abs(V1);
end
figure;
plot(theta_deg, V1_values);
xlabel('Power Factor Angle (degrees)');
ylabel('HV Terminal Voltage (V)');
title('a) HV Voltage vs Power Factor Angle at Full Load');
grid on;
% Part b: Efficiency and regulation for different load factors
load_factors = [1.0, 0.5, 0.25];
colors = ['b', 'g', 'r'];
labels = {'Full Load', 'Half Load', 'Quarter Load'};
% Efficiency plot
figure;
hold on;
for i = 1:length(load_factors)
lf = load_factors(i);
I2 = I2_full * lf;
I_load_prime_mag = I2 / a;
efficiencies = zeros(size(pf));
for k = 1:length(theta_values)
theta = theta_values(k);
I_load_prime = I_load_prime_mag * exp(-1j * theta);
V_drop = I_load_prime * (R_total + 1j * X_total);
V1 = V1_nom + V_drop;
V1_mag = abs(V1);
P_out = V2 * I2 * cos(theta);
copper_loss = (I_load_prime_mag^2) * R_total;
core_loss = (V1_mag^2) / Rc;
efficiency = P_out / (P_out + copper_loss + core_loss) * 100;
efficiencies(k) = efficiency;
end
plot(rad2deg(theta_values), efficiencies, colors(i), 'DisplayName', labels{i});
end
xlabel('Power Factor angle');
ylabel('Efficiency (%)');
title('b) Efficiency vs Power Factor');
legend;
grid on;
hold off;
% Regulation plot
figure;
hold on;
for i = 1:length(load_factors)
lf = load_factors(i);
I2 = I2_full * lf;
I_load_prime_mag = I2 / a;
regulations = zeros(size(pf));
for k = 1:length(theta_values)
theta = theta_values(k);
I_load_prime = I_load_prime_mag * exp(-1j * theta);
V_drop = I_load_prime * (R_total + 1j * X_total);
V1 = V1_nom + V_drop;
V1_mag = abs(V1);
regulation = (V1_mag - V1_nom) / V1_nom * 100;
regulations(k) = regulation;
end
plot(rad2deg(theta_values), regulations, colors(i), 'DisplayName', labels{i});
end
xlabel('Power Factor angle');
ylabel('Voltage Regulation (%)');
title('b) Voltage Regulation vs Power Factor');
legend;
grid on;
hold off;
%tapping check
figure;
hold on;
tapping_upper = V1_nom * 1.05;
tapping_lower = V1_nom * 0.95;
for i = 1:length(load_factors)
lf = load_factors(i);
I2 = I2_full * lf;
V1_values = zeros(size(theta_values));
for k = 1:length(theta_values)
theta = theta_values(k);
I_load_prime = (I2 / a) * exp(-1j * theta);
V_drop = I_load_prime * (R_total + 1j * X_total);
V1 = V1_nom + V_drop;
V1_values(k) = abs(V1);
end
% Plotting
plot(theta_deg, V1_values, colors(i), 'DisplayName', labels{i});
% Indices where voltage is within tapping range
idx_valid = find(V1_values <= tapping_upper & V1_values >= tapping_lower);
fprintf('Load factor: %.2f\n', lf);
if ~isempty(idx_valid)
min_valid_angle = theta_deg(min(idx_valid));
max_valid_angle = theta_deg(max(idx_valid));
fprintf(' ✅ HV within ±5%% tapping from %.1f° to %.1f°\n\n', min_valid_angle, max_valid_angle);
else
fprintf(' ❌ HV voltage exceeds ±5%% tapping for all power factor angles.\n\n');
end
end

yline(tapping_upper, '--r', 'Upper Tapping Limit', 'LabelHorizontalAlignment', 'left',...
'LabelVerticalAlignment', 'bottom');
yline(tapping_lower, '--b', 'Lower Tapping Limit', 'LabelHorizontalAlignment', 'left',...
'LabelVerticalAlignment', 'top');
xlabel('Power Factor Angle (degrees)');
ylabel('HV Terminal Voltage (V)');
title('d) HV Voltage vs Power Factor Angle at Full, Half, and Quarter Load');
legend show;
grid on;
hold off