clear; clc;

%% === datainput ===
filepath = "riverdata.xlsx";
opts = detectImportOptions(filepath);
opts.VariableNamingRule = "preserve";   % 
T = readtable(filepath, opts);

Sample_Cl_K  = T.("Cl_K");   % N×1
Sample_d41K  = T.("d41K");   % N×1

%% === set paramaters ===
max_iter = 10000; 
j = 1;
f1_list = [];      
rng('shuffle');    

for i=1:height(T)
i
%% === endmembers ===
while j <= max_iter

    Cl_K_w = 1.3 + rand * (7.36 - 1.3);
    d41K_w = -0.49 + rand * (0.33);  %  -0.49 to -0.16
    

    Cl_K_s = 153.65 + rand * (366.71 - 153.65);
    d41K_s = -0.23 + rand * (0.6); 


    f_num = Sample_Cl_K(i) - Cl_K_s;
    f_den = Cl_K_w - Cl_K_s;

    if f_den == 0  
        continue;
    end

    f = f_num / f_den;
    f_comp = 1 - f;


    d41K_mix = f * d41K_w + (1 - f) * d41K_s;


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

%% === output ===
if ~isempty(f1_list)
    f1_mean(i) = mean(f1_list);
    f1_std(i)  = std(f1_list);
     f2_mean(i) = mean(f2_list);
    f2_std(i)  = std(f2_list);

else
    disp('cannot solve for f');
end
end

