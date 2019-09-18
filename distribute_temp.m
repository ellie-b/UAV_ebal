%% Calculate temperature from AWS temperature, absorbed radiation, and distributed absorbed radiation

function [dist_temp] = distribute_temp(SWabs,T)

Tc = T-273.15;
tf = 9.7721;
Tres_aws = Tc-(SWabs(1885,910)*tf); % cell containing AWS

dist_temp = (SWabs*tf) + Tres_aws + 273.15;