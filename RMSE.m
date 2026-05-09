function [RMSE_poloha, RMSE_rychlost] = RMSE(estimate,real)

chyba = real - estimate;
chyba = chyba.^2;
RMSE_poloha = mean(sqrt(chyba(1, :) + chyba(2,:)));
RMSE_rychlost = mean(sqrt(chyba(3, :)+ chyba(4,:)));
end