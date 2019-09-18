%% Function to define start and end of each month, 2004-2010

function [start_month,end_month,ndays_year]=calendar(firstyear,nmonths,nyears);

for n=1:nyears
    year=firstyear+n-1;
    % length of each month
    ndays(n,1)=31;
    if (mod(year,4)==0)    % leap year
        ndays(n,2)=29;
    else
        ndays(n,2)=28;
    end
    ndays(n,3)=31;
    ndays(n,4)=30;
    ndays(n,5)=31;
    ndays(n,6)=30;
    ndays(n,7)=31;
    ndays(n,8)=31;
    ndays(n,9)=30;
    ndays(n,10)=31;
    ndays(n,11)=30;
    ndays(n,12)=31;
    ndays_year(n)=sum(ndays(n,:));
end
    
% define the start and end dates for each month, decimal days
for n=1:nyears   % 1   % first year
    start_month(n,1)=1;       % Jan 1
    end_month(n,1)=31;        % Jan 31
    for m=2:nmonths
       start_month(n,m)=start_month(n,m-1)+ndays(n,m-1);
       end_month(n,m)=start_month(n,m)+ndays(n,m)-1;
    end
end

return

for n=2:nyears
    start_month(n,1)=end_month(n-1,12)+1;           % Jan 1
    end_month(n,1)=start_month(n,1)+ndays(n,1)-1;   % Jan 31
    for m=2:nmonths
       start_month(n,m)=start_month(n,m-1)+ndays(n,m-1);
       end_month(n,m)=start_month(n,m)+ndays(n,m)-1;
    end
end