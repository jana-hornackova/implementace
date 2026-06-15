function [x_f,P_f, R] = KF_R(u, y, x_0, P_0, A, B, C, Q, Sigma_0, s_0)
steps = width(y);
x_pred = cell([1 steps]);
x_f = cell([1 steps]);
P_pred = cell([1 steps]);
P_f = cell([1 steps]);
s = cell([1 steps]);
Sigma = cell([1 steps]);
y = num2cell(y, 1);


x_pred{1} = x_0;
P_pred{1} = P_0;
s{1} = s_0+1;
Sigma{1} = Sigma_0;

iters = 100; %pocet iteraci IVB algoritmu
I = eye(length(x_0));

for k = 1:steps
   s{k+1} = s{k}+1;
   %datovy krok
   P = cell(1, iters);
   x = cell(1, iters);
   Sigma_i = cell(1, iters);
   Sigma_i{1} = Sigma{k};
       
   for i = 2:iters
       R = (s{k}*Sigma_i{i-1}^-1)^-1; % \hat{R^-1}^-1
       K = P_pred{k}*C'/(R+C*P_pred{k}*C');
       P{i} = (I-K*C)*P_pred{k}*(I-K*C)'+K*R*K';
       x{i} = x_pred{k} + K*(y{k}-C*x_pred{k});

       Sigma_i{i} = Sigma{k} + C*P{i}*C' + (y{k}-C*x{i})*(y{k}-C*x{i})';
   end
   x_f{k} = x{end};
   P_f{k} = P{end};
   Sigma{k} = Sigma_i{end};

   %casovy krok
   P_pred{k+1} = Q + A*P_f{k}*A';
   x_pred{k+1} = A*x_f{k} + B*u{k};
   Sigma{k+1} = Sigma{k};
end
end