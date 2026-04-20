function [x_f,P_f] = STF_R(u, y, x_0, P_0, A, B, C, Q, nu, Sigma_0, s_0)

steps = width(y);
x_pred = cell([1 steps]);
x_f = cell([1 steps]);
P_pred = cell([1 steps]);
P_f = cell([1 steps]);
Sigma = cell([1 steps]);
s = cell([1 steps]);
y = num2cell(y, 1);

x_pred{1} = x_0;
P_pred{1} = P_0;
Sigma{1} = Sigma_0;
s{1} = s_0;

iters = 10; %pocet iteraci IVB algoritmu
mathring_y = height(y);
m = nu + mathring_y;

for k = 1:steps
   %datovy krok
   P = cell(1, iters);
   x = cell(1, iters);
   Sigma_i = cell(1, iters);
   Sigma_i{1} = Sigma{k};
   l = cell(1, iters);

   P{1} = P_pred{k};
   x{1} = x_pred{k};
   
   s{k+1} = s{k}+1;
   for i = 2:iters
       l{i} = nu+trace(s{k}*Sigma_i{i-1}^-1*((y{k}-C*x{i-1})*(y{k}-C*x{i-1})'+C*P{i-1}*C'));
       Sigma_i{i} = Sigma{k}+m/l{i}*((y{k}-C*x{i-1})*(y{k}-C*x{i-1})'+C*P{i-1}*C');
       
       P{i} = (P_pred{k}^-1+C'*s{k}/Sigma_i{i}*m/l{i}*C)^-1;
       x{i} = P{i}*(P_pred{k}^-1*x_pred{k} + C'*s{k}/Sigma_i{i}*m/l{i}*y{k});
   end
   x_f{k} = x{end};
   P_f{k} = P{end};
   Sigma{k} = Sigma_i{end};

   %casovy krok
   P_pred{k+1} = Q + A*P_f{k}*A';
   x_pred{k+1} = A*x_f{k} + B*u(k);
   Sigma{k+1} = Sigma{k};
end

x_f = cell2mat(x_f);

end