%% Add hourly melt for July 21-24 (noon to noon)

location = 'G:\Fountain_Fieldwork\yr2016\MeltModel\UAV_ebal_data\DistributedMelt\';
list = dir(location);
filenames = {list.name}; % first two cells blank
sday = 203;
eday = 206;

total_melt = zeros(3294,7792);

nr = length(filenames);
for i=3:nr
    strings = split(filenames{i},{'_','.'});
    day = str2num(strings{2});
    hour = str2num(strings{3});
    if day==sday && hour>=12
        disp([day hour])
        [melt,R] = geotiffread([location filenames{i}]);
    elseif day==eday && hour<=12
        disp([day hour])
        [melt,R] = geotiffread([location filenames{i}]);
    elseif day>sday && day<eday
        disp([day hour])
        [melt,R] = geotiffread([location filenames{i}]);
    else
        melt = zeros(3294,7792);
    end
    
    total_melt = total_melt + melt;
    
    n_amelt(i-2)=melt(1885,911);
end

total_melt = total_melt/-1000;

imagesc(total_melt)
contour(total_melt)
contourf(total_melt)
pcolor(total_melt)
surf(total_melt, 'edgecolor', 'none'); view(2);

total_melt(total_melt<-1000) = -9999;
geotiffwrite('G:\Fountain_Fieldwork\yr2016\MeltModel\UAV_ebal_data\Final Outputs\jul21_24_ebal.tif',total_melt,R);
    