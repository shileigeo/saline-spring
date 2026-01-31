clear; clc;

%% === 用户输入河流样品参数 ===
filepath = "riverdata.xlsx";
opts = detectImportOptions(filepath);
opts.VariableNamingRule = "preserve";   % 保留原始列名（更稳）
T = readtable(filepath, opts);

Sample_Cl_K  = T.("Cl_K");   % N×1
Sample_d41K  = T.("d41K");   % N×1

%% === 初始化参数与变量 ===
max_iter = 10000;  % 模拟次数
j = 1;
f1_list = [];      % 存储有效 f 值
rng('shuffle');    % 保证每次运行随机性

for i=1:height(T)
i
%% === 蒙特卡洛模拟混合端元 ===
while j <= max_iter
    % 模拟水体端元 (Cl/K 和 δ41K)
    Cl_K_w = 1.3 + rand * (7.36 - 1.3);
    d41K_w = -0.49 + rand * (0.33);  % 即 -0.49 到 -0.16
    
    % 模拟盐泉端元 (Cl/K 和 δ41K)
    Cl_K_s = 153.65 + rand * (366.71 - 153.65);
    d41K_s = -0.23 + rand * (0.6);  % 即 -0.23 到 0.37

    % 计算混合比例 f（来源于水体端元）
    f_num = Sample_Cl_K(i) - Cl_K_s;
    f_den = Cl_K_w - Cl_K_s;

    if f_den == 0  % 避免除零错误
        continue;
    end

    f = f_num / f_den;
    f_comp = 1 - f;

    % 计算理论 δ41K
    d41K_mix = f * d41K_w + (1 - f) * d41K_s;

    % 计算差值平方，用于筛选模拟点
    d_lN = (Sample_d41K(i) - d41K_mix)^2;

    if f > 0 && f < 1 && d_lN < 1e-4
        f1_list(j) = f;
        f2_list(j)=f_comp;
        Cl_K_w_list(j) = Cl_K_w;
        d41K_w_list(j) = d41K_w;
        Cl_K_s_list(j) = Cl_K_s;
        d41K_s_list(j) = d41K_s;
        j = j + 1;
    end
    
end

j=1;

%% === 输出结果 ===
if ~isempty(f1_list)
    f1_mean(i) = mean(f1_list);
    f1_std(i)  = std(f1_list);
     f2_mean(i) = mean(f2_list);
    f2_std(i)  = std(f2_list);
%     fprintf('\n模拟结果:\n');
%     fprintf('f₁（水体端元比例）平均值为 %.4f\n', f1_mean);
%     fprintf('标准差为 %.4f\n', f1_std);
else
    disp('未找到满足条件的 f 值，可能是输入参数与端元范围不匹配。');
end
end
%% === 可视化：δ41K vs Cl/K ===
% figure;
% hold on; grid on;
% scatter(Cl_K_w_list, d41K_w_list, 15, 'b', 'filled', 'DisplayName', '水体端元');
% scatter(Cl_K_s_list, d41K_s_list, 15, 'r', 'filled', 'DisplayName', '盐泉端元');
% scatter(Sample_Cl_K, Sample_d41K, 100, 'k', 'filled', 'DisplayName', '样品');
% 
% xlabel('Cl/K 比');
% ylabel('\delta^{41}K (‰)');
% title('河水样品与端元混合图');
% legend('Location','best');
