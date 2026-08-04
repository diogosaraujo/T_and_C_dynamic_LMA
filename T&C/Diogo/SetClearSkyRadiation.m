function [Eb1,Eb2,Ed1,Ed2,Edp1,Edp2,rho_s1,rho_s2,Mb,Mg] = SetClearSkyRadiation(h_S,Zbas,Tdew,So1,So2,beta_A,alpha_A,omega_A1,omega_A2,uo,un,rho_g)
%%%% INPUT
%%h_S,[rad]  solar altitude
%Zbas [m a.s.l.]  watershed elevation
%Tdew [°C] dew point temperature
%So1 [W/m^2]  Extraterrestrial  Radiation VIS band [0.29 um - 0.70 um ]
%So2  % [W/m^2] %%Extraterrestrial Radiation NIR band [0.70 um - 4.0 um ]
%beta_A   0.05 [0-1] Angstrom turbidity parameters
%alpha_A  1.3 [0.5-2.5]Angstrom turbidity parameters
%omega_A1  %% 0.92 [0.74-0.94]aerosol single-scattering albedo band 1
%omega_A2 %% 0.84 [0.74-0.94] aerosol single-scattering albedo band 2
%%% uo 0.35 [0.22-0.35]  %% [cm] ozone amount in vertical column
%%% un 0.0002 [0.0001-0.046] %% [cm] total nitrogen dioxide amount
%%%  rho_g 0.15 [0.05-0.3] %% [] spatial average regional albedo
%%% OUTPUT
%%Eb1 [W/m^2]  clear sky beam irradiance VIS band [0.29 um - 0.70 um ]
%%Eb2  [W/m^2]  clear sky beam irradiance NIR band [0.70 um - 4.0 um ]
%%Ed1 [W/m^2]  clear sky total diffuse flux at the ground  VIS band [0.29 um - 0.70um ]
%%%Ed2   [W/m^2]  clear sky total diffuse flux at the ground NIR band [0.70 um - 4.0um ]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h_SD = 180/pi*h_S; %% [°]  solar altitude
Z=90-h_SD; %% sun zenit angle [°]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(h_S > 0.0) %% if there is sunshine
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DIRECT-BEAM TRANSMITTANCES %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Generally p/po = exp(-g*d_elevat/(Rd*Tm));  Rd =287.05; %% [J/kgK]  g= 9.81; %% [m/s^2] Tm = 15°C
    p=1013.25*exp(-Zbas/8434.5); %%% [mbar] Pressure after correction for differences in pressure between basin and seal level
    w=  exp(0.07*Tdew - 0.075); %% precipitable water [cm]
    %%% Gueymard 2003
    mR=(sin(h_S)+(0.48353*Z^(0.095846))/(96.741-Z)^1.1754)^-1; %%% Rayleigh scattering and uniformly mixd gas  Air mass
    mO=(sin(h_S) +(1.0651*Z^0.6379)/((101.8-Z)^2.2694))^-1; %% ozone absorption Air mass
    %% not necessary %% mn=(sin(h_S) +(1.1212*Z^1.6132)/((111.55-Z)^3.2629))^-1; %%
    mW=(sin(h_S) +(0.10648*Z^0.11423)/((93.781-Z)^1.9203))^-1;%% water vapor Air mass
    mA=(sin(h_S) +(0.16851*Z^0.18198)/((95.318-Z)^1.9542))^-1; %% aerosol extinction Air mass
    mRp=(p/1013.25)*mR;%% Rayleigh scattering and uniformely mixd gas  Air mass Correction pressure
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute ozone transmittances
    f1= uo*(10.979-8.5421*uo)/(1+2.0115*uo+40.189*uo^2);
    f2= uo*(-0.027589-0.005138*uo)/(1-2.4857*uo+13.942*uo^2);
    f3= uo*(10.995-5.5001*uo)/(1+1.6784*uo + 42.406*uo^2);
    TO1 = (1+f1*mO+f2*mO^2)/(1+f3*mO);
    TO2 = 1.0;
    % Compute nitrogen dioxide transmittances
    g1 = (0.17499 +41.654*un -2146.4*un^2)/(1+22295.0*un^2);
    g2= un*(-1.2134+ 59.324*un)/(1+ 8847.8*un^2);
    g3= (0.17499 +61.658*un + 9196.4*un^2)/(1+74109.0*un^2);
    TN1 = min(1,(1+g1*mW+g2*mW^2)/(1+g3*mW));
    TN2 = 1.0;
    % Compute Rayleigh scattering transmittances
    TR1    = (1+ 1.8169*mRp -0.033454*mRp^2)/(1+ 2.063*mRp +0.31978*mRp^2);
    TR2    = (1 -0.010394*mRp)/(1- 0.00011042*mRp^2);
    % Compute the uniformly mixd gas transmittances
    TG1    = (1+ 0.95885*mRp -0.012871*mRp^2)/(1+ 0.96321*mRp +0.015455*mRp^2);
    TG2    = (1+ 0.27284*mRp -0.00063699*mRp^2)/(1+0.30306*mRp);
    % Compute water vapor transmmitances
    h1=w*(0.065445+0.00029901*w)/(1+1.2728*w);
    h2=w*(0.065687+0.0013218*w)/(1+1.2008*w);
    c1=w*(19.566-1.6506*w+1.0672*w^2)/(1+ 5.4248*w+1.6005*w^2);
    c2=w*(0.50158-0.14732*w+0.047584*w^2)/(1+ 1.1811*w+1.0699*w^2);
    c3=w*(21.286-0.39232*w+1.2692*w^2)/(1+ 4.8318*w+1.412*w^2);
    c4=w*(0.70992-0.23155*w+0.096514*w^2)/(1+ 0.44907*w+0.75425*w^2);
    TW1    = (1+h1*mW)/(1+h2*mW);
    TW2    = (1+c1*mW+c2*mW^2)/(1+c3*mW+c4*mW^2);
    % Compute average aerosol transmittances
    alpha1 = alpha_A; %%Angstrom turbidity parameters exponent
    alpha2 = alpha_A; %%Angstrom turbidity parameters exponent
    beta1 = beta_A*0.7^(alpha1-alpha2); %Angstrom turbidity parameters
    beta2 = beta_A ; %Angstrom turbidity parameters
    uA     = log(1 + mA*beta2);%%% ---> Gueymard, (1989) Parameterization for effective wavelength computation
    d0    = 0.57664 - 0.024743*alpha1 ;
    d1    = (0.093942 -0.2269*alpha1 + 0.12848*alpha1^2)/(1+0.6418*alpha1);
    d2 = (-0.093819 + 0.36668*alpha1 - 0.12775*alpha1^2)/(1-0.11651*alpha1);
    d3 = alpha1*(0.15232-0.087214*alpha1+0.012664*alpha1^2)/(1-0.90454*alpha1+0.26167*alpha1^2);
    e0 =(1.183 -0.022989*alpha2 + 0.020829*alpha2^2)/(1+0.11133*alpha2);
    e1 =(-0.50003 -0.18329*alpha2 + 0.23835*alpha2^2)/(1+1.6756*alpha2);
    e2 =(-0.50001 +1.1414*alpha2 + 0.0083589*alpha2^2)/(1+11.168*alpha2);
    e3 =(-0.70003 -0.73587*alpha2 + 0.51509*alpha2^2)/(1+4.7665*alpha2);
    le1    = (d0 + d1*uA + d2*uA^2)/(1+d3*uA^2);  %% effective wavelength for band 1
    le2    = (e0 + e1*uA + e2*uA^2)/(1+e3*uA);  %% effective wavelength for band 2  %%% <--
    tauA1= beta1*le1^(-alpha1);
    tauA2= beta2*le2^(-alpha2);
    TA1    = exp(-mA*tauA1);
    TA2    = exp(-mA*tauA2);
    % Return the transmittances for each band
    T1(1)  = TO1;
    T1(2)  = TR1;
    T1(3)  = TG1;
    T1(4)  = TW1;
    T1(5)  = TA1;
    T1(6)  = TN1;
    %%%
    T2(1)  = TO2;
    T2(2)  = TR2;
    T2(3)  = TG2;
    T2(4)  = TW2;
    T2(5)  = TA2;
    T2(6)  = TN2;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute clear sky beam irradiance
    Eb1 = So1*prod(T1);
    Eb2 = So2*prod(T2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   DIFFUSE TRANSMITTANCES   %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    mWp= 1.66; % water vapor optical mass of reference
    TW1p    = (1+h1*mWp)/(1+h2*mWp);
    TW2p    = (1+c1*mWp+c2*mWp^2)/(1+c3*mWp+c4*mWp^2);
    TN1p = min(1,(1+g1*mWp+g2*mWp^2)/(1+g3*mWp));
    TN2p = 1.0;
    clear g1 g2 g3 h1 h2
    %%%%%%%%%%
    TAs1 = exp(-mA*omega_A1*tauA1);  %% aerosol transmittances due to scattering
    TAs2 = exp(-mA*omega_A2*tauA2);  %% aerosol transmittances due to scattering
    %%%%%%%%%%%%%%%%%%%%%%%
    BR1   = 0.5*(0.89013-0.0049558*mR+0.000045721*mR^2);  %% forward scattering fractions for Rayleigh extinction
    BR2  =  0.5 ; %% forward scattering fractions for Rayleigh extinction
    Ba   = 1 - exp(-0.6931 - 1.8326*sin(h_S));%% fractions aerosol scattered fluxes
    %%%%%%%%%%
    g0= (3.715 +0.368*mA +0.036294*mA^2)/(1+0.0009391*mA^2);
    g1= (-0.164-0.72567*mA +0.20701*mA^2)/(1+0.0019012*mA^2);
    g2= (-0.052288+0.31902*mA +0.17871*mA^2)/(1+0.0069592*mA^2);
    h0= (3.4352 +0.65267*mA +0.00034328*mA^2)/(1+0.034388*mA^1.5);
    h1= (1.231 -1.63853*mA +0.20667*mA^2)/(1+0.1451*mA^1.5);
    h2= (0.8889 -0.55063*mA +0.50152*mA^2)/(1+0.14865*mA^1.5);
    F1=(g0+g1*tauA1)/(1+g2*tauA1);
    F2=(h0+h1*tauA2)/(1+h2*tauA2);
    %%%% Compute the clear  sky albedo %%%%
    rho_s1=(0.13363 +0.00077358*alpha1 ...
        +beta1*(0.37567+0.22946*alpha1)/(1-0.10832*alpha1))/(1+beta1*(0.84057+0.68683*alpha1)/(1-0.08158*alpha1));
    rho_s2=(0.010191 +0.00085547*alpha2 ...
        +beta2*(0.14618+0.062758*alpha2)/(1-0.19402*alpha2))/(1+beta2*(0.58101+0.17426*alpha2)/(1-0.17586*alpha2));
    % Compute incidente irradiance a perfectly absorbing ground
    Edp1= TO1*TG1*TN1p*TW1p*(BR1*(1-TR1)*(TA1^0.25) +Ba*F1*TR1*(1-TAs1^0.25))*So1*sin(h_S);
    Edp2= TO2*TG2*TN2p*TW2p*(BR2*(1-TR2)*(TA2^0.25) +Ba*F2*TR2*(1-TAs2^0.25))*So2*sin(h_S);
    %% Backscattered diffuse component
    Edd1= rho_g*rho_s1*(Eb1*sin(h_S) +Edp1)/(1-rho_g*rho_s1);
    Edd2= rho_g*rho_s2*(Eb2*sin(h_S) +Edp2)/(1-rho_g*rho_s2);
    %%%% %%% Compute diffuse irradiances
    Ed1= Edp1 +Edd1;
    Ed2= Edp2 +Edd2;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Compute PAR Ratios
    m15 = min(mR,15);
    t0 = (0.90227 + 0.29*m15 + 0.22928*(m15^2)-0.0046842*(m15^3))/(1 +0.35474*m15 +0.19721*(m15^2));
    t1 = (-0.10591 + 0.15416*m15 - 0.048486*(m15^2)+0.0045932*(m15^3))/(1 -0.29044*m15 +0.026267*(m15^2));
    t2 = (0.47291 - 0.44639*m15 + 0.1414*(m15^2)-0.014978*(m15^3))/(1 -0.37798*m15 +0.052154*(m15^2));
    t3 = (0.077407+ 0.18897*m15 - 0.072869*(m15^2)+ 0.0068684*(m15^3))/(1 -0.25237*m15 +0.020566*(m15^2));
    v0 = (0.82725+ 0.86015*m15 + 0.00713*(m15^2)+ 0.00020289*(m15^3))/(1 +0.90358*m15 +0.015481*(m15^2));
    v1 = (-0.089088+ 0.089226*m15 - 0.021442*(m15^2)+ 0.0017054*(m15^3))/(1 -0.28573*m15 +0.024153*(m15^2));
    v2 = (-0.05342- 0.0034387*m15 + 0.0050661*(m15^2)- 0.00062569*(m15^3))/(1 -0.32663*m15 +0.029382*(m15^2));
    v3 = (-0.17797+ 0.13134*m15 - 0.030129*(m15^2)+ 0.0023343*(m15^3))/(1 -0.28211*m15 +0.023712*(m15^2));
    beta_e=beta1*(le1^(1.3-alpha1));
    Mb=(t0 +t1*beta_e +t2*(beta_e^2))/(1+t3*(beta_e^2)); %% PAR ratio for beam
    Mg=(v0 +v1*beta_e +v2*(beta_e^2))/(1+v3*(beta_e^2)); %% PAR ratio for global
else %%% not sunshine
    Eb1 = 0;
    Eb2 = 0;
    Ed1 = 0;
    Ed2 = 0;
    Edp1=0;
    Edp2=0;
    rho_s1=NaN;
    rho_s2=NaN;
    Mb=0;
    Mg=0;
end
end