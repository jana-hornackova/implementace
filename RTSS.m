function [x_s, P_s] = RTSS(u, y, x_0, P_0, A, B, C, Q, R)
steps = width(y);
x_pred = cell([1 steps]);
x_f = cell([1 steps]);
P_pred = cell([1 steps]);
P_f = cell([1 steps]);
y = num2cell(y, 1);

x_pred{1} = x_0;
P_pred{1} = P_0;

I = eye(length(x_0));
%dopredna rekurze
for k = 1:steps
   %datovy krok
   K = P_pred{k}*C'/(R+C*P_pred{k}*C');
   P_f{k} = (I-K*C)*P_pred{k}*(I-K*C)'+K*R*K';
   x_f{k} = x_pred{k} + K*(y{k}-C*x_pred{k}); 

   %casovy krok
   P_pred{k+1} = Q + A*P_f{k}*A';
   x_pred{k+1} = A*x_f{k} + B*u{k};
end

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
end
