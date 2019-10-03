%% Plot temoerature matrix
Temp=dTK-273.15;
imagesc(Temp)
contour(Temp)
contourf(Temp)
pcolor(Temp)
surf(Temp, 'edgecolor', 'none'); view(2);

geotiffwrite('G:\Fountain_Fieldwork\yr2016\MeltModel\UAV_ebal_data\Temp_diagnosis\dT_Aug1_1885911.tif',Temp,R)

%% compare temperatures
Tjul = T(9014:9758);
aTp = aT';

plot(Tjul,'r')
hold on
plot(aT,'b')

%% dTK stats
mean(Temp(isinf(Temp)==0))
max(Temp(isinf(Temp)==0))
min(Temp(isinf(Temp)==0))

%% check calculations
calc_diagnosis = [SWabs_mod Taws dT zeros(size(Taws))];
for i = 1:744
    calc_diagnosis(i,4) = (SWabs_mod(i)*9.7721)+(Taws(i)-(SWabs_mod(i)*9.7721));
end

%% create row and column maps
rows = zeros(size(ice_albedo));
[r,c] = size(ice_albedo);
for i = 1:r
    rows(i,:)=i;
end

cols = zeros(size(ice_albedo));
for i = 1:c
    cols(:,i)=i;
end

% geotiffwrite('G:\Fountain_Fieldwork\yr2016\MeltModel\UAV_ebal_data\AWS_loc\columns.tif',cols,R)
% geotiffwrite('G:\Fountain_Fieldwork\yr2016\MeltModel\UAV_ebal_data\AWS_loc\rows.tif',rows,R)