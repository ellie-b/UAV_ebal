%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate glacier energy balance at the Bylot AWS site

% Read in the AWS data
% Headers are ignored by the matlab function xlsread
clear all       % clears matlab memory
data=xlsread('Bylot_AWSdata.xlsx');      % AWS data from 2015-2016
ndata=size(data);
ndat=ndata(1,1);

albedofile1 = 'j21albLLadj.tif';
albedofile2 = 'j23albLLadj.tif';

ndays=365;
nmonths=12;
nyears=2;
firstyear=2015;
[start_month,end_month]=calendar(firstyear,nmonths,nyears);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start procssing the AWS data
year=data(:,1);
dday=data(:,2);
Tn=data(:,5);
T=data(:,6);
Tx=data(:,7);
RHn=data(:,8);
RHx=data(:,9);
SWin=data(:,10);
SWout=data(:,11);
wind=data(:,17);

TK=T+273.15;
RH=(RHn+RHx)/2;

% convert SW data from MJ/m2 to W/m2
mj2w=1e6/3600;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Process the AWS data for the 1-hour energy balance

% constants for the energy balance model - these need to be set properly
zAWS=368;     % altitude of the AWS site m
PAWS=969;     % average air presure of the AWS site, mbar
ice_albedo=geotiffread(albedofile);

% for potential temperature and specific humidity calculations
grav=9.81;    % m/s2
Pref=1000;    % mbar
R=288.5;      % gas law constant
eps=0.622;    % constant for vapour pressure to specific humidity
evsurf=6.11;  % mbar, assumes saturated, melting
cp=1010;      % heat capacity of air, J/kg K
gamma=R/cp;
Lv=2514;      % latent heat of vaporization, J/g
Lf=335000;    % latent heat of vaporization, J/kg
rhow=1000;    % density of water, kg/m3
sigma=5.67e-8;   % Stefan-Boltzmann constant
epss=1.0;       % emissivity for melting snow,ice
% constants for turbulent fluxes
z=1.5;           % assume constant, 1.5 m, for now
z0=0.001;        % 1 mm, a customary number
z0H=z0/100;
z0E=z0H;         % roughness length scales (m)
kvk=0.4;         % von Karman
Tmelt=273.15;    % K

% calculate air density and specific humidity, AWS site
for n=1:ndat
   % calculate air density, ev, and qv at the AWS site
   curday=floor(dday(n));
   rho(n)=PAWS*100/(R*TK(n));    % kg/m3
   % humidity
   % saturation vapour pressure in mbar
   es(n)=6.112*exp(17.62*T(n)/(243.12+T(n)));  % mbar: WMO (2008), over water
   ev(n)=es(n)*RH(n)/100;    % actual vapour pressure, mbar;
   qv(n)=eps*ev(n)*1000/(PAWS-ev(n));   % specific humidity (g/kg) 
end

% miscellaneous QC or diagnostic thresholds
pdd_threshold=0.;
frost_threshold=0.;
SW_threshold=10;      % W/m2
max_albedo=0.9;

% set up solar radiation file location
location = 'July2016_all\';
list = dir(location);
filenames = {list.name}; % first two cells blank

%% calculate the 1-hour and daily energy balance at all points
Nt=24;      % integrate 1-hour data; Nt = time steps per day
cumPDD = 0;
totalmelt = 0;
% Just do calculations for July 2016 (dday 183-215)
jul1=9014;
aug1=9758;
for n=jul1:aug1
    k=n-jul1+1;    % index for hourly data
    day(k)=floor(dday(n));   % current day
   
    % hourly radiation balance
    dSWin=geotiffread([location filenames{k+2}]); % hourly incoming radiation
    SWnet=dSWin*mj2w.*(1-ice_albedo);
    
    % calculate temperature distribution
    dTK = distribute_temp(SWnet,TK(n));

    % model the LW radiation from T, ev and RH (Ebrahimi, JGR, 2015)
    epsa=min(0.445+0.0055*RH(n)+0.0052*ev(n),1);
    LWin=epsa*sigma*dTK.^4;
    LWout=epss*sigma*Tmelt^4;
    LWnet=LWin-LWout;
     
    Qstar=SWnet+LWnet;
       
    % Calculate the turbulent heat fluxes on the glacier
    % need the surface temperature
    TKs=Tmelt;    % assume a melting surface, unless data is available
    TCs = TKs - 273.15;
    % potential temperature 
    Tpot=dTK*(Pref/PAWS)^gamma;
    Tpots=TKs*(Pref/PAWS)^gamma;
    % Calculate turbulent fluxes with no stability corrections
    zm=1.5;     % measurement height, assume 1.5 m
    z0=0.001;   % 1 mm, a customary number
    % saturation vapour pressure in mbar
    evsurf=6.112*exp(22.46*TCs/(276.62+TCs));  % mbar: WMO (2008), over water
    % humidity
    qvs=eps*evsurf*1000/(PAWS-evsurf);     % specific humidity (g/kg) 
    z0H=z0/100;
    z0E=z0H;      % roughness length scales for QH and QE (m)
    denom=log(zm/z0)*log(zm/z0H);

    QH=rho(n)*cp*kvk^2*wind(n)*(Tpot-Tpots)/denom;
    QE=rho(n)*Lv*kvk^2*wind(n)*(qv(n)-qvs)/denom;
    Qnet=Qstar+QH+QE;
    Enet=Qnet*3600;      % J/m2
    
    Emelt=Enet; % melt energy available
    Emelt(Emelt<=0)=0; % refreezing or cooling
    
    melt = Emelt*1000/(rhow*Lf);    % mm melt/hour
    totalmelt = totalmelt+melt;
     
    % calculate melting or refreezing       
    if (T(n) >= 0)     % surface at 0 degC
       PDD(k) = T(n)/Nt;  % PDD from this hour 
    else
       % air is below 0 degC, could be refeeezing
       PDD(k) = 0;
    end
    cumPDD = cumPDD+PDD(k);
    
    % save AWS location
    aSWnet(k)=SWnet(1885,910);
    aLWin(k)=LWin(1885,910);
    aLWout(k)=LWout;
    aLWnet(k)=LWnet(1885,910);
    aQstar(k)=Qstar(1885,910);
    aQH(k)=QH(1885,910);
    aQE(k)=QE;
    aQnet(k)=Qnet(1885,910);
    aEnet(k)=Enet(1885,910);
    aEmelt(k)=Emelt(1885,910);
    amelt(k)=melt(1885,910);

end

%% mean July diagnostics
ndat_July=31*24;
TJul=mean(T(jul1:aug1))
meltJul=totalmelt(1885,910)
PDDJul=sum(PDD)
SWnetJul=mean(aSWnet)
LWinJulA=mean(aLWin)
LWoutJul=mean(aLWout)
QstarJul=mean(aQstar)
QHJul=mean(aQH)
QEJul=mean(aQE)
QnetJul=mean(aQnet)

%% Calculate 12-hour diagnostics from the month
np=31*2;       % number of periods, month of July
nfields=16;    % output fields, model results
Bylot_results=zeros(np,nfields);
ndat_period=12;    % 12 data points per period
for m=1:np
    startdat=1+(m-1)*ndat_period;   % start of 12-hr period
    enddat=m*ndat_period;           % end of 12-hr period
    startdatd=startdat+jul1-1;      % index, original dataset
    enddatd=startdat+jul1-1;        % index, original dataset
    % averages or sums for this period
    Bylot_results(m,1)=m;   % period (counter)
    Bylot_results(m,2)=day(enddat);   % day in July
    Bylot_results(m,3)=2-mod(m,2);      % 1/2 for am/pm
    Bylot_results(m,4)=mean(T(startdatd:enddatd));     % degC
    Bylot_results(m,5)=sum(PDD(startdat:enddat));    % W/m2
    Bylot_results(m,6)=mean(aSWnet(startdat:enddat));  % W/m2
    Bylot_results(m,8)=ice_albedo(1885,910);  % 
    Bylot_results(m,9)=mean(aLWin(startdat:enddat));    % W/m2
    Bylot_results(m,10)=mean(aLWout(startdat:enddat));  % W/m2
    Bylot_results(m,11)=mean(aQstar(startdat:enddat));  % W/m2
    Bylot_results(m,12)=mean(aQH(startdat:enddat));     % W/m2
    Bylot_results(m,13)=mean(aQE(startdat:enddat));     % W/m2
    Bylot_results(m,14)=mean(aQnet(startdat:enddat));   % W/m2
    Bylot_results(m,15)=sum(aEmelt(startdat:enddat)/1e3);  % kJ/m2
    Bylot_results(m,16)=sum(amelt(startdat:enddat));    % mm   
end
        
save 'Bylot_DistEbalance_July2016.dat' Bylot_results -ascii   

% save 'Distributed_Melt_July21.24.dat' totalmelt -ascii

%% Calculate daily diagnostics from the month
np=31;       % number of periods, month of July
nfields=15;    % output fields, model results
Bylot_results_daily=zeros(np,nfields);
ndat_period=24;    % 24 data points per period
for m=1:np
    startdat=1+(m-1)*ndat_period;   % start of 12-hr period
    enddat=m*ndat_period;           % end of 12-hr period
    startdatd=startdat+jul1-1;      % index, original dataset
    enddatd=startdat+jul1-1;        % index, original dataset
    % averages or sums for this period
    Bylot_results_daily(m,1)=m;   % period (counter)
    Bylot_results_daily(m,2)=day(enddat);   % day in July
    Bylot_results_daily(m,3)=mean(T(startdatd:enddatd));     % degC
    Bylot_results_daily(m,4)=sum(PDD(startdat:enddat));    % W/m2
    Bylot_results_daily(m,5)=mean(SWin(startdatd:enddatd));  % W/m2
    Bylot_results_daily(m,6)=mean(SWout(startdatd:enddatd)); % W/m2
    Bylot_results_daily(m,7)=mean(albedo(startdat:enddat));  % 
    Bylot_results_daily(m,8)=mean(LWin(startdat:enddat));    % W/m2
    Bylot_results_daily(m,9)=mean(LWout(startdat:enddat));  % W/m2
    Bylot_results_daily(m,10)=mean(Qstar(startdat:enddat));  % W/m2
    Bylot_results_daily(m,11)=mean(QH(startdat:enddat));     % W/m2
    Bylot_results_daily(m,12)=mean(QE(startdat:enddat));     % W/m2
    Bylot_results_daily(m,13)=mean(Qnet(startdat:enddat));   % W/m2
    Bylot_results_daily(m,14)=sum(Emelt(startdat:enddat)/1e3);  % kJ/m2
    Bylot_results_daily(m,15)=sum(melt(startdat:enddat));    % mm   
end
%         
% save 'Bylot_Ebalance_July2016_daily.dat' Bylot_results_daily -ascii        
% 
% subplot(3,1,1)
% plot(Bylot_results_daily(:,1),Bylot_results_daily(:,3))
% ylabel('T (\circC)')
% axis([1 31 -1 14])
% 
% subplot(3,1,2)
% hold off
% plot(Bylot_results_daily(:,1),Bylot_results_daily(:,10),'r')
% hold on
% plot(Bylot_results_daily(:,1),Bylot_results_daily(:,11)+Bylot_results_daily(:,12),'g')
% plot(Bylot_results_daily(:,1),Bylot_results_daily(:,13),'k')
% ylabel('Fluxes (W/m^2)')
% axis([1 31 -10 270])
% 
% subplot(3,1,3)
% plot(Bylot_results_daily(:,1),Bylot_results_daily(:,15))
% xlabel('period (July 2016)')
% ylabel('melt (mm)')
% axis([1 31 12 69])
% 
% print -djpeg90 'Bylot_daiymelt_July2016.jpg'
