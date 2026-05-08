function [x_s,P_s, R] = KS_R(u, y, x_0, P_0, A, B, C, Q, Sigma_0, nu_0)
steps = width(y);

Sigma = cell([1 steps]);
y = num2cell(y, 1);

Sigma{1} = Sigma_0;
iters = 10; %pocet iteraci IVB algoritmu

nu = nu_0 + steps;
for i = 1:iters
R = (nu*Sigma{i}^-1)^-1;   
[x_s, P_s] = RTSS(u, cell2mat(y), x_0, P_0, A, B, C, Q, R);
Sigma{i+1} = Sigma_0;
for k = 1:steps
    Sigma{i+1} = Sigma{i+1} + C*P_s{k}*C' + (y{k}-C*x_s{k})*(y{k}-C*x_s{k})';
end
end
end