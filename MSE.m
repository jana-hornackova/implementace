function [MSE_poloha, MSE_rychlost] = MSE(x_est,x_real)

chyba = x_real - x_est;
chyba = chyba.^2;
MSE_poloha = mean(chyba(1, :) + chyba(2,:));
MSE_rychlost = mean(chyba(3, :)+ chyba(4,:));
end