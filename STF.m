function [x_f,P_f] = STF(u, y, x_0, P_0, A, B, C, Q, nu, R)
mathring_y = height(y);

steps = width(y);
x_pred = cell([1 steps]);
x_f = cell([1 steps]);
P_pred = cell([1 steps]);
P_f = cell([1 steps]);
y = num2cell(y, 1);

x_pred{1} = x_0;
P_pred{1} = P_0;

iters = 100; %pocet iteraci IVB algoritmu
m = nu + mathring_y;
I = eye(length(x_0));

for k = 1:steps
   %datovy krok
   P = cell(1, iters);
   x = cell(1, iters);
   l = cell(1, iters);
   
   P{1} = P_pred{k};
   x{1} = x_pred{k};
    
   for i = 2:iters
       l{i} = nu + trace(((y{k}-C*x{i-1})*(y{k}-C*x{i-1})'+C*P{i-1}*C')/R);
       
       exp_cov = (R^-1*m*l{i}^-1)^-1;
       K = P_pred{k}*C'/(exp_cov+C*P_pred{k}*C');
       P{i} = (I-K*C)*P_pred{k}*(I-K*C)'+K*exp_cov*K';
       x{i} = x_pred{k} + K*(y{k}-C*x_pred{k});
   end
   x_f{k} = x{end};
   P_f{k} = P{end};

   %casovy krok
   P_pred{k+1} = Q + A*P_f{k}*A';
   x_pred{k+1} = A*x_f{k} + B*u{k};
end

end