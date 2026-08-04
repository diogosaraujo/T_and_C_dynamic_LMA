clear

fl = 'DailyData/';
fn = '45EM01';
cl = 1; % column - station - Just a for loop will get the job done for the rest

lon = [0.621601021; 0.763286489; 0.715583095; -0.018710311; -1.000866466; -1.470792551; ...
    -1.193064098; -1.283380229; -1.870021066; 0.084718318; -0.495372536; 0.148339227; ...
    -1.65194137; -1.754633382; -2.415435891; -2.338664342; -3.403767709; -2.780650893; ...
    -0.522415626; -0.488100921; -0.038555892; -2.676990592; -2.094256698; -5.21697241; ...
    -4.544712986; -4.287939068; -3.684625963; -1.611305601; -2.854473465; -3.116782759; ...
    -2.041749486; -4.349674597; -4.86967249];

lat = [51.61056055; 52.03904372; 52.93052051; 52.81915283; 52.70561195; 53.01431871; ...
    53.51620256; 54.41567571; 54.09415143; 51.55837072; 51.46831946; 51.37734342; ...
    54.57009638; 55.49602964; 55.10913606; 53.42855051; 54.6793072; 54.54111457; ...
    51.52262249; 51.2523841; 51.03884694; 50.99938727; 51.7926654; 50.21030892; ...
    52.80840545; 51.86936734; 51.66372352; 52.3045886; 52.87756796; 52.0931356; ...
    57.52687991; 58.34972028; 58.43797408];

% where the station is
Lon = lon(cl);
Lat = lat(cl);
Zbas = 100;

% hour definition
DeltaGMT = 0;
t_bef = 0;
t_aft = 1;

% representative atmospheric parameters - BETTER IF THEY WERE DEPENDENT ON
% MONTH!
beta_A= 0.04;
uo  = 0.35; %% 0.28;  %%% [cm] ozone amount in vertical column [0.22-0.3]
un=  0.0002 ; % [cm] total nitrogen dioxide amount
alpha_A= 1.3;%% Angstrom turbidity parameters
%%% Clear Sky diffuse
omega_A1     = 0.92;%%0.920;%%%omega_A1  %% aerosol single-scattering albedo band 1 [0.84-0.94]
omega_A2     = 0.84;%%0.833;%%omega_A2 %% aerosol single-scattering albedo band 2  [0.84-0.94]
rho_g = 0.15; %% [] spatial average regional albedo

So1 = 800;% [W/m^2]  Extraterrestrial  Radiation VIS band [0.29 um - 0.70 um ]
So2 = 300;% % [W/m^2] %%Extraterrestrial Radiation NIR band [0.70 um - 4.0 um ]

%
hmin = 5; % hour of minimum temperature
hmax = 15;% hour of maximum temperature

% Define Date (this needs to be adjusted. this is just an example)

Date_in = (datenum([1981,1,1]):1:datenum([2079,12,30]))';

Date_out = (datenum([1981,1,1]):(1/24):datenum([2079,12,30,23,0,0]))';

Date_360 = [];
for y = 1981:2079
    for m = 1:12

        if m < 3 || m > 5
            Date_360 = [Date_360; (datenum([y,m,1])+(0:29))'];
        elseif m == 3 % to avoid duplicate dates
            Date_360 = [Date_360; (datenum([y,m,3])+(0:29))'];
        elseif m == 4 % to avoid duplicate dates
            Date_360 = [Date_360; (datenum([y,m,2])+(0:29))'];
        elseif m == 5 % to avoid duplicate dates
            Date_360 = [Date_360; (datenum([y,m,2])+(0:29))'];
        end

    end
end

%% Load Daily and fill in problematic days

U_d = load([fl,'hurs_data_rcp',fn,'.mat'],'varData');%varData is a terrible name convention
U_d = U_d.varData(:, cl)/100;
Pr_d = load([fl,'pr_data_rcp',fn,'.mat'],'varData');%varData is a terrible name convention
Pr_d = Pr_d.varData(:, cl)*3600*24;% not sure - check units
S_d = load([fl,'rsds_data_rcp',fn,'.mat'],'varData');%varData is a terrible name convention
S_d = S_d.varData(:, cl);
Ta_d = load([fl,'tas_data_rcp',fn,'.mat'],'varData');%varData is a terrible name convention
Ta_d = Ta_d.varData(:, cl)-273.15;
Tmin_d = load([fl,'tasmin_data_rcp',fn,'.mat'],'varData');%varData is a terrible name convention
Tmin_d = Tmin_d.varData(:, cl)-273.15;
Tmax_d = load([fl,'tasmax_data_rcp',fn,'.mat'],'varData');%varData is a terrible name convention
Tmax_d = Tmax_d.varData(:, cl)-273.15;
Ws_d = load([fl,'sfcWind_data_rcp',fn,'.mat'],'varData');%varData is a terrible name convention
Ws_d = Ws_d.varData(:, cl);

% Interpolate daily - fill in the missing days because of the 360
% simulation calendar
U_d = interp1(Date_360, U_d, Date_in);
Pr_d = interp1(Date_360, Pr_d, Date_in);
S_d = interp1(Date_360, S_d, Date_in);
Ta_d = interp1(Date_360, Ta_d, Date_in);
Tmin_d = interp1(Date_360, Tmin_d, Date_in);
Tmax_d = interp1(Date_360, Tmax_d, Date_in);
Ws_d = interp1(Date_360, Ws_d, Date_in);


%% Interpolate hourly

Date_in_mid = Date_in + 0.5;

% Pressure? - Barometric correction
Pre_ref = 1000*exp(-9.80665*0.0289644*(Zbas)/8.3144598/(mean(Ta_d+273.15)));
Pre = Pre_ref*ones(length(Date_out), 1);%hPa

% Interpolate windspeed - assume constant during the day
Ws = interp1(Date_in_mid, Ws_d, Date_out, 'nearest');

% Interpolate Humidity - neglects diurnal variability
U = interp1(Date_in_mid, U_d, Date_out, 'nearest');

% Get temperature time series with min-max
Ta = NaN(size(Date_out));
Date_in_min = Date_in + hmin/24;
Date_in_max = Date_in + hmax/24;
ix_min = zeros(size(Date_in_min));
ix_max = zeros(size(Date_in_min));
for i = 1:length(ix_max)
    ix_min(i) = find(abs(Date_out-Date_in_min(i))<1e-3);
    ix_max(i) = find(abs(Date_out-Date_in_max(i))<1e-3);
end
Ta(ix_min) = Tmin_d;
Ta(ix_max) = Tmax_d;
% then interpolate for nans
ix = ~isnan(Ta);
Ta = interp1(Date_out(ix), Ta(ix), Date_out,'pchip');

%% Alternative for Ta - only works for DeltaGMT=0 currently

for i = 1:length(Tmin_d)

    Datam = datevec(Date_in(i));

    [delta_S,h_S,zeta_S,T_sunrise,T_sunset,L_day,E0,jDay,Delta_TSL] ...
        = SetSunVariables(Datam,DeltaGMT,Lon,Lat,t_bef,t_aft);

    Date_in_min(i) = Date_in(i) + T_sunrise/24;%coldest temperature at dawn
    Date_in_max(i) = Date_in(i) + 15/24;%(T_sunrise+L_day/2+2)/24;%warmest 2 hrs after midday

    ix_min(i) = find(abs(Date_out-Date_in_min(i))==min(abs(Date_out-Date_in_min(i))));
    ix_max(i) = find(abs(Date_out-Date_in_max(i))==min(abs(Date_out-Date_in_max(i))));

end
Ta = NaN(size(Date_out));

Ta(ix_min) = Tmin_d;
Ta(ix_max) = Tmax_d;
% then interpolate for nans
ix = ~isnan(Ta);
Ta = interp1(Date_out(ix), Ta(ix), Date_out,'pchip');

plot(Ta,'kx-')
hold on
plot(ix_min,Tmin_d,'ro')
plot(ix_max,Tmax_d,'bs')

%%
% Get dew point
Tdew = (243.12*(log(U) + (17.625.*Ta)./(Ta+243.04)) )./(17.62 - log(U)+(17.625.*Ta)./(Ta+243.04));
Tdew(Tdew>Ta) = Ta(Tdew>Ta);

% Precipitation

Pr = zeros(length(Date_out), 1);

for j = 1:length(Date_in)
    if Pr_d(j) > 0

        stp = ceil(Pr_d(j)/10);

        ix = find(abs(Date_out-Date_in(j)) == min(abs(Date_out-Date_in(j))));
        for k = 1:stp
            Pr(ix+k-1) = Pr_d(j)/stp;
        end
    end
end

% and finally radiation partitioning

Rsw_c = zeros(length(Date_out), 1); % clear sky radiation - for getting the diurnal shape right

for j = 1:length(Date_out)

    Datam = datevec(Date_out(j));
    [delta_S,h_S,zeta_S,T_sunrise,T_sunset,L_day,E0,jDay,Delta_TSL] = SetSunVariables(Datam(1:4),DeltaGMT,Lon,Lat,t_bef,t_aft);
    [Eb1,Eb2,Ed1,Ed2,Edp1,Edp2,rho_s1,rho_s2,Mb,Mg] = SetClearSkyRadiation(h_S,Zbas,Tdew(j),So1,So2,beta_A,alpha_A,omega_A1,omega_A2,uo,un,rho_g);
    Rsw_c(j) = Eb1 + Eb2 + Ed1 + Ed2;

end

Rsw = zeros(size(Rsw_c));
for j = 1:length(Date_in) % This will only work if we start at hour - midnight

    st = find(Date_out==Date_in(j));
    en = st + 23;

    sm = S_d(j);
    prfl = Rsw_c(st:en);
    Rsw(st:en) = prfl/mean(prfl)*sm;
end

% Partition radiation

% [~,~,SAD1,SAD2,SAB1,SAB2,PARB,PARD,N,~,~,~]=...
%     Automatic_Radiation_Partition(Date,Lat,Lon,Zbas,DeltaGMT,Pr,Tdew,Rsw,0);

[~,~,SAD1,SAD2,SAB1,SAB2,PARB,PARD,N,~,~,~] = ...
    C_Automatic_Radiation_Partition(Date_out,Lat,Lon,Zbas,DeltaGMT,Pr,Tdew,Rsw,0,t_bef,t_aft);
