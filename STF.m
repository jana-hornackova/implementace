function [x_f,P_f] = STF(u, y, x_0, P_0, A, B, C, Q, nu, R)
steps = width(y);
x_pred = cell([1 steps]);
x_f = cell([1 steps]);
P_pred = cell([1 steps]);
P_f = cell([1 steps]);
y = num2cell(y, 1);

x_pred{1} = x_0;
P_pred{1} = P_0;

iters = 10; %pocet iteraci IVB algoritmu
mathring_y = height(y);
s = nu + mathring_y;

for k = 1:steps
   %datovy krok
   P = cell(1, iters);
   x = cell(1, iters);
   Sigma = cell(1, iters);
   
   P{1} = P_pred{k};
   x{1} = x_pred{k};
    
   for i = 2:iters
       Sigma{i} = nu + trace(((y{k}-C*x{i-1})*(y{k}-C*x{i-1})'+C*P{i-1}*C')/R);
       
       P{i} = (C'/R*s/Sigma{i}*C + P_pred{k}^-1)^-1;
       x{i} = P{i}*(P_pred{k}^-1*x_pred{k} + C'/R*s/Sigma{i}*y{k});
   end
   x_f{k} = x{end};
   P_f{k} = P{end};

   %casovy krok
   P_pred{k+1} = Q + A*P_f{k}*A';
   x_pred{k+1} = A*x_f{k} + B*u(k);
end

x_f = cell2mat(x_f);

end