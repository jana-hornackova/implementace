function [MSE_poloha] = MSE(estimate,real)

chyba = real - estimate;
chyba = chyba.^2;
MSE_poloha = mean(chyba(1, :) + chyba(2,:));
%MSE_rychlost = mean(chyba(3, :)+ chyba(4,:));
end