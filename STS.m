function [x_s,P_s] = STS(u, y, x_0, P_0, A, B, C, Q, nu, R)
mathring_y = height(y);
steps = width(y);
x_pred = cell([1 steps]);
x_f = cell([1 steps]);
P_pred = cell([1 steps]);
P_f = cell([1 steps]);
y = num2cell(y, 1);
Sigma = cell([1 steps]);
Sigma(:) = {nu};

x_pred{1} = x_0;
P_pred{1} = P_0;

iters = 10; %pocet iteraci IVB algoritmu

s = nu + mathring_y;
I = eye(length(x_0));

for i = 1:iters

%dopredna rekurze
for k = 1:steps
   %datovy krok
   exp_cov = (R^-1*s*Sigma{k}^-1)^-1;
   K = P_pred{k}*C'/(exp_cov+C*P_pred{k}*C');
   P_f{k} = (I-K*C)*P_pred{k}*(I-K*C)'+K*exp_cov*K';
   x_f{k} = x_pred{k} + K*(y{k}-C*x_pred{k});
   %P_f{k} = (C'/R*s/Sigma{k}*C + P_pred{k}^-1)^-1;
   %x_f{k} = P_f{k}*(P_pred{k}^-1*x_pred{k} + C'/R*s/Sigma{k}*y{k});

   %casovy krok
   P_pred{k+1} = Q + A*P_f{k}*A';
   x_pred{k+1} = A*x_f{k} + B*u(k);
end

x_pred{end} = [];
P_pred{end} = [];

%zpetna rekurze
x_s = cell(size(x_f));
x_s{end} = x_f{end};

P_s = cell(size(P_f));
P_s{end} = P_f{end};

for k = (length(x_f)-1) : -1 :1
    G_k = P_f{k} * A' /P_pred{k+1}; 
    x_s{k} = x_f{k} + G_k*(x_s{k+1}-x_pred{k+1});
    P_s{k} = P_f{k} + G_k*(P_s{k+1}-P_pred{k+1})*G_k';
end 

for k = 1:steps
    Sigma{k} = nu + trace(((y{k}-C*x_s{k})*(y{k}-C*x_s{k})'+C*P_s{k}*C')/R);
end

end
end