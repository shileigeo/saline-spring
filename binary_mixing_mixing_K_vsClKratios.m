%% data input
clear all;
close all;
clc;

load('river_Cl_K.mat');
load('sal_Cl_K.mat');

%% 1. fit Cl_K_mix_riv and 1/K_mix_riv
p = polyfit(Cl_K_mix_riv, 1 ./ K_mix_riv, 1);

slope_fit = p(1);
intercept_fit = p(2);

y_obs = 1 ./ K_mix_riv;
y_pred = polyval(p, Cl_K_mix_riv);

SS_res = sum((y_obs - y_pred).^2);
SS_tot = sum((y_obs - mean(y_obs)).^2);

R2_fit = 1 - SS_res / SS_tot;

fprintf('equation: 1/K = %.6f * Cl_K + %.6f\n', ...
    slope_fit, intercept_fit);
fprintf(' R2 = %.4f\n\n', R2_fit);

fprintf('equation: 1/K = %.6f * Cl_K + %.6f\n\n', ...
    slope_fit, intercept_fit);

%% 2. Cl_K_west_riv average composition A
Cl_west_mean = mean(Cl_K_west_riv);
Cl_west_std  = std(Cl_K_west_riv);

y_west_mean = slope_fit * Cl_west_mean + intercept_fit;
y_west_std  = abs(slope_fit) * Cl_west_std;

fprintf('Point A: (%.4f, %.6f), y error: %.6f\n\n', ...
    Cl_west_mean, y_west_mean, y_west_std);

%% 3. set paramaters
inv_K_EK1 = 1 / K_EK1;
inv_K_EK2 = 1 / K_EK2;

inv_K_sal_mean = mean(1 ./ K_sal);
inv_K_sal_std  = std(1 ./ K_sal);

inv_K_west_min = min(1 ./ K_west_riv);
inv_K_west_max = max(1 ./ K_west_riv);

Cl_sal_min = min(Cl_K_sal);
Cl_sal_max = max(Cl_K_sal);

n_target = 10000;


y_target = inv_K_sal_mean;


x_fixed = 2.3247;

fprintf('sal 1/K mean = %.6f\n', inv_K_sal_mean);
fprintf('sal Cl/K range = [%.4f, %.4f]\n', Cl_sal_min, Cl_sal_max);
fprintf('west river 1/K range = [%.6f, %.6f]\n\n', ...
    inv_K_west_min, inv_K_west_max);

%% 4. MC - EK1


X_EK1 = [];
slope_EK1 = [];
intercept_EK1 = [];
m_EK1 = [];

while length(X_EK1) < n_target

    n_try = 50000;


    X_try = Cl_sal_min + (Cl_sal_max - Cl_sal_min) * rand(n_try, 1);

    valid_denominator = abs(X_try - Cl_K_EK1) > eps;

    X_try = X_try(valid_denominator);


    slope_try = (y_target - inv_K_EK1) ./ (X_try - Cl_K_EK1);
    intercept_try = inv_K_EK1 - slope_try .* Cl_K_EK1;


    m_try = slope_try .* x_fixed + intercept_try;


    valid = m_try >= inv_K_west_min & m_try <= inv_K_west_max;

    X_EK1 = [X_EK1; X_try(valid)];
    slope_EK1 = [slope_EK1; slope_try(valid)];
    intercept_EK1 = [intercept_EK1; intercept_try(valid)];
    m_EK1 = [m_EK1; m_try(valid)];

end


X_EK1 = X_EK1(1:n_target);
slope_EK1 = slope_EK1(1:n_target);
intercept_EK1 = intercept_EK1(1:n_target);
m_EK1 = m_EK1(1:n_target);

X_EK1_mean = mean(X_EK1);
X_EK1_std  = std(X_EK1);

fprintf('EK1 m range: [%.6f, %.6f]\n', min(m_EK1), max(m_EK1));
fprintf('EK1: Cl_K = %.4f ± %.4f\n\n', X_EK1_mean, X_EK1_std);

%% 5. MC - EK2
X_EK2 = [];
slope_EK2 = [];
intercept_EK2 = [];
m_EK2 = [];

while length(X_EK2) < n_target

    n_try = 5000;


    X_try = Cl_sal_min + (Cl_sal_max - Cl_sal_min) * rand(n_try, 1);


    valid_denominator = abs(X_try - Cl_K_EK2) > eps;

    X_try = X_try(valid_denominator);


    slope_try = (y_target - inv_K_EK2) ./ (X_try - Cl_K_EK2);
    intercept_try = inv_K_EK2 - slope_try .* Cl_K_EK2;


    m_try = slope_try .* x_fixed + intercept_try;


    valid = m_try >= inv_K_west_min & m_try <= inv_K_west_max;

    X_EK2 = [X_EK2; X_try(valid)];
    slope_EK2 = [slope_EK2; slope_try(valid)];
    intercept_EK2 = [intercept_EK2; intercept_try(valid)];
    m_EK2 = [m_EK2; m_try(valid)];

end


X_EK2 = X_EK2(1:n_target);
slope_EK2 = slope_EK2(1:n_target);
intercept_EK2 = intercept_EK2(1:n_target);
m_EK2 = m_EK2(1:n_target);

X_EK2_mean = mean(X_EK2);
X_EK2_std  = std(X_EK2);

fprintf('EK2 m range: [%.6f, %.6f]\n', min(m_EK2), max(m_EK2));
fprintf('EK2: Cl_K = %.4f ± %.4f\n\n', X_EK2_mean, X_EK2_std);

%% 6. Predicted X at Y = inv_K_sal_mean 
X_fit_sal = (inv_K_sal_mean - intercept_fit) / slope_fit;
X_fit_sal_std = inv_K_sal_std / abs(slope_fit);

fprintf('equation: Cl_K = %.4f ± %.4f\n\n', ...
    X_fit_sal, X_fit_sal_std);

%% 7. draw
figure('Position', [100, 100, 760, 620]);

hold on;

x_min = min([Cl_K_west_riv(:); Cl_K_mix_riv(:); Cl_K_sal(:); Cl_K_EK1; Cl_K_EK2]);
x_max = max([Cl_K_west_riv(:); Cl_K_mix_riv(:); Cl_K_sal(:); Cl_K_EK1; Cl_K_EK2]);

x_range = linspace(x_min, x_max, 200);

step = max(1, floor(n_target / 80));


for i = 1:step:n_target
    y_mc = slope_EK1(i) .* x_range + intercept_EK1(i);
    plot(x_range, y_mc, '-', ...
        'Color', [0.25 0.75 0.75], ...
        'LineWidth', 0.4, ...
        'HandleVisibility', 'off');
end


for i = 1:step:n_target
    y_mc = slope_EK2(i) .* x_range + intercept_EK2(i);
    plot(x_range, y_mc, '-', ...
        'Color', [0.75 0.25 0.75], ...
        'LineWidth', 0.4, ...
        'HandleVisibility', 'off');
end


scatter(Cl_K_mix_riv, 1 ./ K_mix_riv, 55, ...
    'b', 'filled', ...
    'DisplayName', 'Major mixing river data');


y_fit = slope_fit .* x_range + intercept_fit;
plot(x_range, y_fit, 'b-', ...
    'LineWidth', 2, ...
    'DisplayName', 'equation');


yline(inv_K_sal_mean, 'k--', ...
    'LineWidth', 1.2, ...
    'DisplayName', 'sal mean 1/K');


plot(X_fit_sal, inv_K_sal_mean, 'bo', ...
    'MarkerSize', 10, ...
    'LineWidth', 2, ...
    'DisplayName', 'intersection');


errorbar(Cl_west_mean, y_west_mean, ...
    y_west_std, y_west_std, ...
    Cl_west_std, Cl_west_std, ...
    'ro', ...
    'MarkerSize', 8, ...
    'MarkerFaceColor', 'r', ...
    'LineWidth', 1.5, ...
    'DisplayName', 'Point A');

% EK1 和 EK2
plot(Cl_K_EK1, inv_K_EK1, 'g^', ...
    'MarkerSize', 10, ...
    'MarkerFaceColor', 'g', ...
    'DisplayName', 'EK1');

plot(Cl_K_EK2, inv_K_EK2, 'ms', ...
    'MarkerSize', 10, ...
    'MarkerFaceColor', 'm', ...
    'DisplayName', 'EK2');

% west river 
plot(Cl_K_west_riv, 1 ./ K_west_riv, 'rs', ...
    'MarkerSize', 5, ...
    'MarkerFaceColor', 'r', ...
    'DisplayName', 'west river');

% sal 
plot(Cl_K_sal, 1 ./ K_sal, 'kp', ...
    'MarkerSize', 8, ...
    'MarkerFaceColor', 'k', ...
    'DisplayName', 'sal');

xlabel('Cl/K', 'FontSize', 12);
ylabel('1/K', 'FontSize', 12);

title(sprintf('Monte Carlo simulation: uniform sal Cl/K intersections, n = %d', ...
    n_target), ...
    'FontSize', 13);

legend('Location', 'best');
grid on;
box on;
xlim([0 400]);ylim([0 0.05]);

%% 8. check
figure('Position', [150, 150, 700, 450]);

histogram(X_EK1, 30);
xlabel('EK1 simulated sal Cl/K intersection');
ylabel('Count');
title('Distribution of EK1 sal intersections');
grid on;
box on;

figure('Position', [200, 200, 700, 450]);

histogram(X_EK2, 30);
xlabel('EK2 simulated sal Cl/K intersection');
ylabel('Count');
title('Distribution of EK2 sal intersections');
grid on;
box on;


%% 9. summary
fprintf('\n========== summary ==========\n');
fprintf('equation:  Cl_K = %.4f ± %.4f\n', X_fit_sal, X_fit_sal_std);
fprintf('EK1:     Cl_K = %.4f ± %.4f\n', X_EK1_mean, X_EK1_std);
fprintf('EK2:     Cl_K = %.4f ± %.4f\n', X_EK2_mean, X_EK2_std);

fprintf('\n========== 约束信息 ==========\n');
fprintf(' y_target = %.6f\n', y_target);
fprintf('sal Cl/K range: [%.4f, %.4f]\n', Cl_sal_min, Cl_sal_max);
fprintf('west river 1/K range: [%.6f, %.6f]\n', ...
    inv_K_west_min, inv_K_west_max);

fprintf('\nEK1 m range: [%.6f, %.6f]\n', min(m_EK1), max(m_EK1));
fprintf('EK2 m range: [%.6f, %.6f]\n', min(m_EK2), max(m_EK2));
