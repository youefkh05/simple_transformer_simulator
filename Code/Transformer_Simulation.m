clc;
clear;
close all;

% Transformer Ratings
S_rated = 75000;        % VA
V_HV = 4600;            % V
V_LV = 240;             % V
a = V_HV / V_LV;        % Turns ratio
% Full load secondary current
I2_full = S_rated / V_LV;

% Transformer impedances (in ohms)
R1 = 0.846;
R2 = 0.00261;
X1 = 26.8;
X2 = 0.0745;
Rc = 220000;
Xm = 112000;


% Total equivalent impedance on HV side
[Req, Xeq, Zeq] = refer_secondary_to_HV(R2, X2, a, R1, X1);

%Part a) HV Voltage vs Power Factor
% power factor varies from 0.6 pf leading through unity pf to 0.6 pf lagging
% add some saftey factor
theta_deg =linspace(-60, 60, 100);
theta_rad = deg2rad(theta_deg);
pf = cos(theta_deg);

% plot High Voltage vs Power Factor
plot_HV_voltage_vs_power_factor(I2_full, a, Req, Xeq, V_HV, theta_deg);


% Part b: Efficiency and regulation for different load factors
load_factors = [1.0, 0.5, 0.25];

% plot efficiency vs Power Factor
plot_efficiency_vs_pf(I2_full, a, Req, Xeq, Rc, V_HV, V_LV, load_factors, theta_deg);

% plot regulation vs Power Factor
plot_regulation_vs_pf(I2_full, a, Req, Xeq, V_HV, load_factors, theta_deg);


%tapping check
figure;
hold on;
tapping_upper = V_HV * 1.05;
tapping_lower = V_HV * 0.95;
for i = 1:length(load_factors)
lf = load_factors(i);
I2 = I2_full * lf;
V1_values = zeros(size(theta_rad));
for k = 1:length(theta_rad)
theta = theta_rad(k);
I_load_prime = (I2 / a) * exp(-1j * theta);
V_drop = I_load_prime * (Req + 1j * Xeq);
V1 = V_HV + V_drop;
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

%************* Functions ********%
function [Req, Xeq, Zeq] = refer_secondary_to_HV(R2, X2, a, R1, X1)
% refer_secondary_to_HV - Refer secondary side impedance to the HV side
% 
% Inputs:
%   R2 - Secondary side resistance (Ohms)
%   X2 - Secondary side reactance (Ohms)
%   a  - Turns ratio (V1/V2)
%   R1 - Primary side resistance (Ohms)
%   X1 - Primary side reactance (Ohms)
%
% Outputs:
%   Req - Total equivalent resistance on HV side (Ohms)
%   Xeq - Total equivalent reactance on HV side (Ohms)
%   Zeq - Total equivalent impedance on HV side (complex Ohms)

    % Refer secondary parameters to HV side
    R2_HV = R2 * a^2;
    X2_HV = X2 * a^2;

    % Total equivalent impedance on HV side
    Req = R1 + R2_HV;
    Xeq = X1 + X2_HV;
    Zeq = Req + 1i * Xeq;
end

function plot_HV_voltage_vs_power_factor(I2_full, a, Req, Xeq, V_HV, theta_deg)
% plot_HV_voltage_vs_power_factor - Plots the HV terminal voltage vs. power factor angle at full load
%
% Inputs:
%   I2_full - Full-load secondary current (A)
%   a       - Turns ratio (V1/V2)
%   Req     - Equivalent resistance on HV side (Ohms)
%   Xeq     - Equivalent reactance on HV side (Ohms)
%   V_HV    - Nominal high-voltage terminal (V)
%   theta_deg - power factor angle range
%
% This function calculates and plots the HV terminal voltage for a range of power factor angles
% from 0.6 leading (-54°) to 0.6 lagging (+54°).

    % Define power factor angle range (degrees and radians)
    theta_rad = deg2rad(theta_deg);

    % Preallocate array for V1 values
    V1_values = zeros(size(theta_rad));

    % Compute HV terminal voltage for each angle
    for k = 1:length(theta_rad)
        theta = theta_rad(k);
        I_load_prime = (I2_full / a) * exp(-1j * theta);
        V_drop = I_load_prime * (Req + 1j * Xeq);
        V1 = V_HV + V_drop;
        V1_values(k) = abs(V1);
    end

    % Plotting
    figure;
    plot(theta_deg, V1_values, 'b', 'LineWidth', 1.5);
    xlabel('Power Factor Angle (degrees)');
    ylabel('HV Terminal Voltage (V)');
    title('a) HV Voltage vs Power Factor Angle at Full Load');
    grid on;

end

function plot_efficiency_vs_pf(I2_full, a, Req, Xeq, Rc, V_HV, V_LV, load_factors, theta_deg)
% plot_efficiency_vs_pf - Plots transformer efficiency vs power factor angle
%
% Inputs:
%   I2_full      - Full-load secondary current (A)
%   a            - Turns ratio (V1/V2)
%   Req          - Equivalent resistance on HV side (Ohms)
%   Xeq          - Equivalent reactance on HV side (Ohms)
%   Rc           - Core loss resistance referred to HV side (Ohms)
%   V_HV         - HV terminal voltage (V)
%   V_LV         - LV terminal voltage (V)
%   load_factors - Array of load factors (e.g., [1.0, 0.5, 0.25])
%   theta_deg    - Array of power factor angles in degrees (e.g., linspace(-54, 54, 100))

    theta_rad = deg2rad(theta_deg);
    pf = cos(theta_rad);  % Not used directly, but kept for interpretation

    % Plot setup
    figure;
    hold on;
    colors = lines(length(load_factors));  % dynamic color generation

    for i = 1:length(load_factors)
        lf = load_factors(i);
        I2 = I2_full * lf;
        I_load_prime_mag = I2 / a;
        efficiencies = zeros(size(theta_rad));

        for k = 1:length(theta_rad)
            theta = theta_rad(k);
            I_load_prime = I_load_prime_mag * exp(-1j * theta);
            V_drop = I_load_prime * (Req + 1j * Xeq);
            V1 = V_HV + V_drop;
            V1_mag = abs(V1);

            P_out = V_LV * I2 * cos(theta);
            copper_loss = (I_load_prime_mag^2) * Req;
            core_loss = (V1_mag^2) / Rc;
            P_in = P_out + copper_loss + core_loss;

            efficiency = P_out / P_in * 100;
            efficiencies(k) = efficiency;
        end

        plot(theta_deg, efficiencies, 'Color', colors(i,:), ...
             'DisplayName', sprintf('Load Factor = %.2f', lf), 'LineWidth', 1.5);
    end

    xlabel('Power Factor Angle (degrees)');
    ylabel('Efficiency (%)');
    title('b) Efficiency vs Power Factor Angle');
    legend('Location', 'best');
    grid on;
    hold off;

end

function plot_regulation_vs_pf(I2_full, a, Req, Xeq, V_HV, load_factors, theta_deg)
% plot_regulation_vs_pf - Plots transformer voltage regulation vs power factor angle
%
% Inputs:
%   I2_full      - Full-load secondary current (A)
%   a            - Turns ratio (V1/V2)
%   Req          - Equivalent resistance on HV side (Ohms)
%   Xeq          - Equivalent reactance on HV side (Ohms)
%   V_HV         - HV terminal voltage (V)
%   load_factors - Array of load factors (e.g., [1.0, 0.5, 0.25])
%   theta_deg    - Array of power factor angles in degrees (e.g., linspace(-54, 54, 100))

    theta_rad = deg2rad(theta_deg);

    % Plot setup
    figure;
    hold on;
    colors = lines(length(load_factors));  % dynamic color generation

    for i = 1:length(load_factors)
        lf = load_factors(i);
        I2 = I2_full * lf;
        I_load_prime_mag = I2 / a;
        regulations = zeros(size(theta_rad));

        for k = 1:length(theta_rad)
            theta = theta_rad(k);
            I_load_prime = I_load_prime_mag * exp(-1j * theta);
            V_drop = I_load_prime * (Req + 1j * Xeq);
            V1 = V_HV + V_drop;
            V1_mag = abs(V1);

            regulation = (V1_mag - V_HV) / V_HV * 100;
            regulations(k) = regulation;
        end

        plot(theta_deg, regulations, 'Color', colors(i,:), ...
             'DisplayName', sprintf('Load Factor = %.2f', lf), 'LineWidth', 1.5);
    end

    xlabel('Power Factor Angle (degrees)');
    ylabel('Voltage Regulation (%)');
    title('c) Voltage Regulation vs Power Factor Angle');
    legend('Location', 'best');
    grid on;
    hold off;

end

