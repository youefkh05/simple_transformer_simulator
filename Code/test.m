clc; clear;

% Transformer Ratings
S_rated = 75000;        % VA
V_HV = 4600;            % V
V_LV = 240;             % V
a = V_HV / V_LV;        % Turns ratio

% Transformer impedances (in ohms)
R1 = 0.846;
R2 = 0.00261;
X1 = 26.8;
X2 = 0.0745;
Rc = 220000;
Xm = 112000;

% Refer R1 and X1 to LV side
R1_LV = R1 / a^2;
X1_LV = X1 / a^2;

% Total equivalent impedance on LV side
Req = R1_LV + R2;
Xeq = X1_LV + X2;
Zeq = Req + 1i*Xeq;

% Core losses (approx.)
P_core = (V_LV^2) / Rc;

% Load current at full-load (LV side)
I_full = S_rated / V_LV;

% Power factor angle sweep: from 0.6 leading to 0.6 lagging
pf_vals = 0.6:0.01:1;
theta_lead = -acos(pf_vals); % leading
theta_lag = acos(pf_vals);  % lagging
theta = [theta_lead(end:-1:1), theta_lag(2:end)];
pf_total = [pf_vals(end:-1:1), pf_vals(2:end)];

% Load conditions
load_levels = [1, 0.5, 0.25]; % full, half, quarter
colors = {'r', 'g', 'b'};

% Initialize plots
figure; hold on; grid on;
title('HV Terminal Voltage vs Power Factor Angle');
xlabel('Power Factor'); ylabel('HV Terminal Voltage (V)');

figure_eff = figure; hold on; grid on;
title('Efficiency vs Power Factor'); xlabel('Power Factor'); ylabel('Efficiency (%)');

figure_reg = figure; hold on; grid on;
title('Voltage Regulation vs Power Factor');
xlabel('Power Factor'); ylabel('Regulation (%)');

for i = 1:length(load_levels)
    load_frac = load_levels(i);
    I_load = I_full * load_frac;
    
    V_HV_vals = zeros(1, length(theta));
    eff_vals = zeros(1, length(theta));
    reg_vals = zeros(1, length(theta));
    
    for k = 1:length(theta)
        angle = theta(k);
        pf = pf_total(k);
        
        % Load current phasor
        I = I_load * exp(-1i*angle);
        
        % Voltage drop across Zeq
        V_drop = I * Zeq;
        
        % Source voltage required to maintain 240 V at load
        V_source = 240 + V_drop;
        
        % Referred to HV side
        V_HV_required = abs(V_source) * a;
        V_HV_vals(k) = V_HV_required;
        
        % Output power (W)
        P_out = S_rated * load_frac * pf;
        
        % Copper losses
        P_cu = abs(I)^2 * real(Zeq);
        
        % Input power
        P_in = P_out + P_cu + P_core;
        
        % Efficiency
        eff_vals(k) = (P_out / P_in) * 100;
        
        % Voltage regulation
        V_no_load = abs(240 + 0); % no current, so no drop
        V_full_load = abs(V_source);
        reg_vals(k) = ((V_no_load - V_full_load) / V_full_load) * 100;
    end
    
    % Plot HV voltage
    figure(1);
    plot(pf_total, V_HV_vals, colors{i}, 'DisplayName', sprintf('Load: %.0f%%', load_frac*100));
    
    % Plot efficiency
    figure(figure_eff);
    plot(pf_total, eff_vals, colors{i}, 'DisplayName', sprintf('Load: %.0f%%', load_frac*100));
    
    % Plot regulation
    figure(figure_reg);
    plot(pf_total, reg_vals, colors{i}, 'DisplayName', sprintf('Load: %.0f%%', load_frac*100));
end

figure(1); legend;
figure(figure_eff); legend;
figure(figure_reg); legend;
