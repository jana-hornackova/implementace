clc
clear
close all

%parametry simulace
Ts = 0.2;
TMAX = 50;
t = 0:Ts:TMAX;
steps = length(t);

%system
A=kron([1 Ts; 0 1], eye(2));
B=zeros(4, 1);
C=[eye(2) zeros(2)];

%parametry sumu
%procesni sum
q = 1;
Q = q*kron([Ts^3/3 Ts^2/2;Ts^2/2 Ts], eye(2)); 

%sum mereni
%studentuv
r = 10;
R_st = r * eye(2);
nu = 5;
%normalni
R_n = R_st * nu /(nu-2);


%% simulace - urceni MSE algoritmu
%inicializacni hodnoty
x_0 = [0; 0; 1; 1];
P_0 = eye(4);
Sigma_0 =100* eye(2);
nu_0 = 5;
s_0 = 3;

algoritmy  = ["KF"; "STF"; "KF R"; "STF R"; "RTSS"; "STS"; "KS R"; "STS R"];
odchylky_n = zeros([length(algoritmy) 2]);
odchylky_st = zeros([length(algoritmy) 2]);

iters = 10;
bar = waitbar(0, "0 %");

for j = 1:iters
%generovani vstupu
u = zeros([1 steps]);

%generovani procesniho sumu
w = mvnrnd(zeros(1, length(Q)), Q, steps)';

%generovani stavu
x_real = zeros([height(A) steps]);
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

% estimace stavu pro system zatizeny normalnim sumem
[x_f_KF, P_f_KF] = KF(u, y_n, x_0, P_0, A, B, C, Q, R_n);
[x_s_RTSS, P_s_RTSS] = RTSS(u, y_n, x_0, P_0, A, B, C, Q, R_n);
[x_f_STF, P_f_STF] = STF(u, y_n, x_0, P_0, A, B, C, Q, nu, R_st);
[x_s_STS, P_s_STS] = STS(u, y_n, x_0, P_0, A, B, C, Q, nu, R_st);
[x_f_KF_R, P_f_KF_R, Sigma_f_KF_R, nu_f_KF_R] = KF_R(u, y_n, x_0, P_0, A, B, C, Q, Sigma_0, nu_0);
[x_s_KS_R, P_s_KS_R] = KS_R(u, y_n, x_0, P_0, A, B, C, Q, Sigma_0, nu_0);
[x_f_STF_R, P_f_STF_R] = STF_R(u, y_n, x_0, P_0, A, B, C, Q, nu, Sigma_0, s_0);
[x_s_STS_R, P_s_STS_R] = STS_R(u, y_n, x_0, P_0, A, B, C, Q, nu, Sigma_0, s_0);

vysledky_n = {x_f_KF, x_f_STF, x_f_KF_R, x_f_STF_R, x_s_RTSS, x_s_STS, x_s_KS_R, x_s_STS_R};

for i = 1:numel(vysledky_n)
     [MSE_r, MSE_v] = MSE(vysledky_n{i}, x_real);
     odchylky_n(i, 1)= ((j-1)*odchylky_n(i,1)+MSE_r)/j; 
     odchylky_n(i, 2)= ((j-1)*odchylky_n(i,2)+MSE_v)/j;
end

% estimace stavu pro system zatizeny studentovym sumem
[x_f_KF, P_f_KF] = KF(u, y_st, x_0, P_0, A, B, C, Q, R_n);
[x_s_RTSS, P_s_RTSS] = RTSS(u, y_st, x_0, P_0, A, B, C, Q, R_n);
[x_f_STF, P_f_STF] = STF(u, y_st, x_0, P_0, A, B, C, Q, nu, R_st);
[x_s_STS, P_s_STS] = STS(u, y_st, x_0, P_0, A, B, C, Q, nu, R_st);
[x_f_KF_R, P_f_KF_R, Sigma_f_KF_R, nu_f_KF_R] = KF_R(u, y_st, x_0, P_0, A, B, C, Q, Sigma_0, nu_0);
[x_s_KS_R, P_s_KS_R] = KS_R(u, y_st, x_0, P_0, A, B, C, Q, Sigma_0, nu_0);
[x_f_STF_R, P_f_STF_R] = STF_R(u, y_st, x_0, P_0, A, B, C, Q, nu, Sigma_0, s_0);
[x_s_STS_R, P_s_STS_R] = STS_R(u, y_st, x_0, P_0, A, B, C, Q, nu, Sigma_0, s_0);

vysledky_st = {x_f_KF, x_f_STF, x_f_KF_R, x_f_STF_R, x_s_RTSS, x_s_STS, x_s_KS_R, x_s_STS_R};

for i = 1:numel(vysledky_st)
     [MSE_r, MSE_v] = MSE(vysledky_st{i}, x_real);
     odchylky_st(i, 1)= ((j-1)*odchylky_st(i,1)+MSE_r)/j; 
     odchylky_st(i, 2)= ((j-1)*odchylky_st(i,2)+MSE_v)/j;
end

waitbar(j/iters, bar, num2str(100*j/iters)+ " %");

end

close(bar)

%%
tabulka = table('RowNames',algoritmy);
tabulka.("Normalní šum") = odchylky_n;
tabulka.("Studentův šum") = odchylky_st;
tabulka

%% simulace - testování proti výpadkům měření

%generovani vstupu
u = zeros([1 steps]);

%generovani procesniho sumu
w = mvnrnd(zeros(1, length(Q)), Q, steps)';

%generovani stavu
x_real = zeros([height(A) steps]);
x_real(:, 1) = x_0;
for i = 2:1:steps
    x_real(:, i) = A * x_real(:, i-1) + B*u(i-1) + w(:, i-1);
end

%generovani sumu mereni
e = mvnrnd(zeros(1, length(R_n)), R_n, steps)'; 
t_var = mod(steps, 20):mod(steps, 20):steps; %časy ve kterých naroste variance
e_100 = mvnrnd(zeros(1, length(R_n)), 100*R_n, numel(t_var))';
e(:, t_var) = e_100;

%generovani vystupu
y = zeros([height(C) steps]);
for i = 1:1:steps
   y(:, i) = C*x_real(:, i); 
end

y_test = y + e;

%nalezeni optimalnich parametru
r = 100:1:400;
nu = 3:2:100;

bar = waitbar(0, "0 %");
iters = numel(r);
for i = 1:iters
R_test = r(i)*eye(2);

[x_f_KF, P_f_KF] = KF(u, y_test, x_0, P_0, A, B, C, Q, R_test);
MSE_KF(i) = MSE(x_f_KF, x_real);

for j = 1:numel(nu)
[x_f_STF, P_f_STF] = STF(u, y_test, x_0, P_0, A, B, C, Q, nu(j), R_test);
MSE_STF(i,j) = MSE(x_f_STF, x_real);
end

waitbar(i/iters, bar, num2str(100*i/iters)+ " %");
end
close(bar)
figure
plot(r, MSE_KF)
xlabel("r")
ylabel("MSE")
figure
imshow(MSE_STF, [])
xlabel("\nu")
ylabel("r")