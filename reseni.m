clc
clear
close all

%parametry simulace
Ts = 0.1;
TMAX = 10;
t = 0:Ts:TMAX;
steps = length(t);

%inicializacni hodnoty
x_0 = [0; 0; 1; 1];
P_0 = eye(4);
Sigma_0 = eye(2);
nu_0 = 5;

%system
A=kron([1 Ts; 0 1], eye(2));
B=zeros(4, 1);
C=[eye(2) zeros(2)];

%parametry sumu
%procesni sum
q = 1;
Q = q*kron([Ts^3/3 Ts^2/2;Ts^2/2 Ts], eye(2)); 

%sum mereni
%normalni
r = 10;
R_n = r * eye(2);
%studentuv
R_st = r * eye(2);
nu = 3;

%% simulace
%generovani vstupu
u = zeros([1 steps]);

%generovani procesniho sumu
w = mvnrnd(zeros(1, length(Q)), Q, steps)';

%generovani stavu
x_real = zeros([length(x_0) steps]);
x_real(:, 1) = x_0;
for i = 2:1:steps
    x_real(:, i) = A * x_real(:, i-1) + B*u(i-1) + w(:, i-1);
end

%generovani sumu mereni
e_n = mvnrnd(zeros(1, length(R_n)), R_n, steps)'; %normalni
e_st = mvtrnd(R_st, nu)'; %studentuv

%generovani vystupu
y = zeros([height(C) steps]);
for i = 1:1:steps
   y(:, i) = C*x_real(:, i); 
end

%vystupy zatizene sumem
y_n = y + e_n;
y_st = y + e_st;

%% estimace stavu


[x_f_KF, P_f_KF] = KF(u, y_n, x_0, P_0, A, B, C, Q, R_n);
[x_v_RTSS, P_v_RTSS] = RTSS(u, y_n, x_0, P_0, A, B, C, Q, R_n);
[x_f_stf, P_f_stf] = stf_VB(u, y_st, x_0, P_0, A, B, C, Q, nu, R_st);
[x_f_KF_R, P_f_KF_R] = KF_R(u, y_n, x_0, P_0, A, B, C, Q, Sigma_0, nu_0);

figure
stairs(t, y_n(1,:)')
hold on
stairs(t, x_real(1,:)')
hold on
stairs(t, x_f_stf(1,:)')
grid on
legend("Výstup","Reálná hodnota", "Filtrovaná hodnota")
%%
figure
stairs(t, y_n(1,:)')
hold on
stairs(t, x_real(1,:)')
hold on
stairs(t, x_f_KF(1,:)')
hold on
stairs(t, x_v_RTSS(1,:)')
grid on
legend("Výstup","Reálná hodnota", "Filtrovaná hodnota", "Vyhlazená hodnota")

figure
stairs(t, y_n(1,:)')
hold on
stairs(t, x_real(1,:)')
hold on
stairs(t, x_f_KF_R(1,:)')
hold on
legend("Výstup","Reálná hodnota", "Filtrovaná hodnota")

