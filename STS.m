function [x_s,P_s] = STS(u, y, x_0, P_0, A, B, C, Q, nu, R)
steps = width(y);
x_pred = cell([1 steps]);
x_f = cell([1 steps]);
P_pred = cell([1 steps]);
P_f = cell([1 steps]);
y = num2cell(y, 1);

x_pred{1} = x_0;
P_pred{1} = P_0;

%dopredna rekurze

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

x_s = cell2mat(x_s);

end