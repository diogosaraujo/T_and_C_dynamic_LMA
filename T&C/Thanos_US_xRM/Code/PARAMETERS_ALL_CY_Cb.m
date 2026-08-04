%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% MAIN_FRAME SPATIAL TETHYS %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[aR,Zs,...
    EvL_Zs,Inf_Zs,Bio_Zs,Zinf,RfH_Zs,RfL_Zs,dz,Ks_Zs,Dz,...
    ms,Kbot,Krock,zatm,...
    Ccrown,Cbare,Crock,Curb,Cwat,...
    Color_Class,OM_H,OM_L,PFT_opt_H,PFT_opt_L,d_leaf_H,d_leaf_L,...
    SPAR,Phy,Soil_Param,Interc_Param,SnowIce_Param,VegH_Param,VegL_Param,fpr,...
    VegH_Param_Dyn,VegL_Param_Dyn,...
    Stoich_H,aSE_H,Stoich_L,aSE_L,fab_H,fbe_H,fab_L,fbe_L,...
    ZR95_H,ZR95_L,In_max_urb,In_max_rock,K_usle,...
    Urb_Par,Deb_Par,Zs_deb,... 
    Sllit,Kct,ExEM,ParEx_H,Mpar_H,ParEx_L,Mpar_L]=PARAMETERS_ALL_CY_Cb(code_dir,ANSWER,Psan,Pcla,Porg,md_max,BuildH,Fimp,Fgra,Ftree,FNonveg,UALB,Cbare_)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% PARAMETERS SOILS AND VEGETATION  %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ANSWER; 
%%%%
%%% Tropical Forest
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
ms = 8; %%% Soil Layer
if ANSWER == 0
    b_ix = true;
    ANSWER = 1;
else
    b_ix = false;
end
zatm =[25 25 25 25 25 25 25 25] ; %% Reference Height
zatm = zatm(ANSWER);
if zatm<BuildH
    zatm=BuildH+5;
end
fpr = 1;
%%%%
aR =50;
%Kh=Ks*aR;
Kbot = [0 0 0 0 0 0 0 0]; %% [mm/h] Conductivity at the bedrock layer
Krock = [0 0 0 0 0 0 0 0]; %% [mm/h] Conductivity of Fractured Rock
Kbot= Kbot(ANSWER);
Krock=Krock(ANSWER);
%%%%%%%%%% SOIL INPUT
Color_Class = 0;
%%%%
ExEM = [ 0 0 0 0 0 0 0 0];
ExEM =ExEM(ANSWER); 
%%%%%%%%%

for cc_ = 1:8
    Cwat = 0; 
    Curb = 0.0 ; 
    Crock = 0.0;
    if b_ix
        Cbare = 1;
        Ccrown = 0;
    else
        Cbare = Cbare_;% rough guess
        Ccrown = 1-Cbare_;%[0.7];
    end
    cc = length(Ccrown);%% Crown area
    II = false(1,8);
    II(ANSWER) = true;
end


%%%%%%%%%%%%%%%%%%%
SPAR=2; %%% SOIL PARAMETER TYPE
%%%%
[Osat,L,Pe,Ks,O33,rsd,lan_dry,lan_s,cv_s,K_usle]=Soil_parameters(Psan,Pcla,Porg);
%%%%%
rsd=rsd*ones(1,ms);
lan_dry=lan_dry*ones(1,ms);
lan_s =lan_s*ones(1,ms);
cv_s = cv_s*ones(1,ms);
%%%
%nVG=L+1;
%alpVG = 1/(-101.9368*Pe); %%[1/mm]%;
p=3+2/L;
m=2/(p-1); nVG= 1/(1-m);
alpVG=(((-101.9368*Pe)*(2*p*(p-1))/(p+3))*((55.6+7.4*p+p^2)/(147.8+8.1*p+0.092*p^2)))^-1; %%[1/mm]%;
%%%
Ks=Ks; 
L=L; 
%%%
Osat=Osat*ones(1,ms);
%Ohy = Ohy*ones(1,ms) ; %% [-]
L=L*ones(1,ms);
Pe = Pe*ones(1,ms);
O33 = O33*ones(1,ms);
alpVG= alpVG*ones(1,ms); %% [1/mm]
nVG= nVG*ones(1,ms); %% [-]
Ks_Zs= Ks*ones(1,ms); %%[mm/h]
%%%%%%%%%%%%%%%% Matric Potential
Kfc = 0.2; %% [mm/h]
Phy = 10000; %% [kPa]
%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
[Ofc,Oss_Lp,Owp_Lp,Ohy]=Soil_parametersII(ms,Osat,L,Pe,Ks_Zs,O33,nVG,alpVG,Kfc,1,1,Phy);
%[s_SVG,bVG]=Soil_parameters_VG(Phy,Osat,Ohy,nVG,alpVG,0);
%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
Zs = [0 10 50 100 200 400 800 1000 2500]; %%% [ms+1] 

if  not(length(Zs)==ms+1)
    disp('SOIL LAYER MESH INCONSISTENT')
    return
end
Zdes = 10; %%% Evaporation depth
Zinf = 10; %%% Infiltration depth
Zbio = 250; 
[EvL_Zs]=Evaporation_layers(Zs,Zdes); %%% Evaporation Layer fraction 
[Inf_Zs]=Evaporation_layers(Zs,Zinf); %%% Infiltration Depth Layer fraction
[Bio_Zs]=Evaporation_layers(Zs,Zbio); 
dz= diff(Zs); %%%% [mm]  Thickness of the Layers
Dz=zeros(1,ms); 
for ii = 1:ms
    if ii>1
        Dz(ii)= (dz(ii)+ dz(ii-1))/2; %%% Delta Depth Between Middle Layer  [mm]
    else
        Dz(ii)=dz(1)/2; %%% Delta Depth Between First Middle Layer and soil surface [mm]
    end
end
%%%%%%%%%%%%
% lVGM =0.5*ones(1,ms);
% nVGM = 3.0*ones(1,ms); %  %% [-]
% alpVGM= 33.38*alpVG; %% [1/mm]
% Omac = 0.00*Osat; % 0.021*ones(1,ms); 
% Ks_mac = 100*Ks_Zs; %%[mm/h]
% ZSS = 600; %%% Depth of soil-structural effects 
% %%%%
% etas=3./ZSS; %%
% ssw = exp(-etas*cumsum(Dz)); ssw(cumsum(Dz)>ZSS)=0;
% Ks_mac=Ks_mac.*ssw;
% Omac=Omac.*ssw;
%%%%%%%%%%%%%%%%% OTHER PARAMETER
In_max_urb=5;
In_max_rock=2; %% [mm]


%%%%%%%%%%%%% SNOW PARAMETER
TminS=-1.1;%% Threshold temperature snow
TmaxS= 2.3;%% Threshold temperature snow
ros_max1=580; %520 600; %%% [kg/m^3]
ros_max2=300; %320 450; %%% [kg/m^3]
Th_Pr_sno = 8; %%% [mm/day] Threshold Intensity of snow to consider a New SnowFall
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% ICE Parameter
Ice_wc_sp =0.01; %% [-] Specific Maximum water content ice
ros_Ice_thr = 500 ; %% [kg/m^3] Density Thrshold to transform snow into ice
Aice = 0.28; %% [-] Ice albedo
WatFreez_Th = -8; %% [°C] Threshold for freezing lake water
dz_ice = 0.45; %% [mm / h] Water Freezing Layer progression without snow-layer
%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% PARAMETERS VEGETATION
%%% cc -- number of crown area
%%% Root Depth
CASE_ROOT=1;  %%% Type of Root Profile 
%%%
ZR95_H = [700   800   200   200   400   900   900   800]; %% [mm]
ZR95_L = [0     0     0     0     0     0     0     0]; %% [mm]
ZR50_H = [NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN];
ZR50_L = [NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN];
ZRmax_H = [NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN];
ZRmax_L = [NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN];
ZR95_H =ZR95_H(II); ZR50_H =ZR50_H(II); ZRmax_H =ZRmax_H(II);
ZR95_L =ZR95_L(II); ZR50_L =ZR50_L(II); ZRmax_L =ZRmax_L(II);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Kct=0.75; %%% Factor Vegetation Cover --- for throughfall
%5 Interception Parameter
gcI=3.7; %%% [1/mm]
KcI=0.06; %%%% [mm] -- Mahfouf and Jacquemin 1989
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Interception Parameter
Sp_SN_In= 5.9; %% [mm/LAI]
Sp_LAI_H_In= [.1 .1 .1 .1 .1 .1 .1 .1]; %%[mm/LAI]
Sp_LAI_L_In= [0.2 0.2 .2 .2 .2 .2 .2 .2]; %%[mm/LAI]
Sp_LAI_H_In =Sp_LAI_H_In(II);
Sp_LAI_L_In =Sp_LAI_L_In(II);
%%%%%%%%%%% Leaf Dimension
d_leaf_H= [2     3     1     2     2     4     2     7]; %%[cm]
d_leaf_L= [0 0 0 0 0 0 0 0 ];  %% [cm]
d_leaf_H =d_leaf_H(II);
d_leaf_L =d_leaf_L(II);
%%%%%%%% Biochemical parameter
KnitH=[0.3000    0.5000    0.3000    0.3000    0.3000    0.3000    0.3000    0.4000]; %%% Canopy Nitrogen Decay
KnitL=[0 0 0 0 0 0 0 0];
mSl_H = [0 0 0 0 0 0 0 0];%% [m2 PFT /gC]  Linear increase in Sla with LAI 
mSl_L = [0 0 0 0 0 0 0 0];  % 
KnitH =KnitH(II);  mSl_H =mSl_H(II); 
KnitL =KnitL(II);  mSl_L =mSl_L(II);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%  Photosynthesis Parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FI_H=[0.0810    0.0810    0.0810    0.0810    0.0810    0.0810    0.0810    0.0810];% Intrinsec quantum Efficiency [umolCO2/umolPhotons]
Do_H=[600   900   900   600   900   600   600   600] ; %%[Pa]
a1_H=[4     4     7     4     6     5     4     6];
go_H=[0.0010    0.0050    0.0050    0.0010    0.0100    0.0010    0.0010    0.0010];% [mol / s m^2] minimum Stomatal Conductance
CT_H=[3 3 3 3 3 3 3 3]; %%--> 'CT' == 3  'CT' ==  4  %% Photosyntesis Typology for Plants
DSE_H =[0.6490    0.6490    0.6560    0.6490    0.6560    0.6490    0.6490    0.6490 ];  %% [kJ/mol] Activation Energy - Plant Dependent
Ha_H =[72    72    72    72    72    72    72    72]; %% [kJ / mol K]  entropy factor - Plant Dependent
gmes_H=[Inf   Inf   Inf   Inf   Inf   Inf   Inf   Inf]; %% [mol CO2 / s m^2 ];  mesophyll conductance
rjv_H=[2.1000    2.1000    2.4000    2.1000    2.2000    2.1000    2.1000    2.1000]; %%% Scaling Jmax - Vmax  [umol electrons / umolCO2 ]
%%%
FI_H=FI_H(II); Do_H=Do_H(II); a1_H=a1_H(II); go_H=go_H(II);
CT_H=CT_H(II); DSE_H=DSE_H(II); Ha_H=Ha_H(II); gmes_H=gmes_H(II);
rjv_H=rjv_H(II);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------
FI_L=[0.0810    0.0810    0.0810    0.0810    0.0810    0.0810    0.0810    0.0810];% Intrinsec quantum Efficiency [umolCO2/umolPhotons]
Do_L=[600   900   900   600   900   600   600   600] ; %%[Pa]
a1_L=[4     4     7     4     6     5     4     6];
go_L=[0.0010    0.0050    0.0050    0.0010    0.0100    0.0010    0.0010    0.0010];% [mol / s m^2] minimum Stomatal Conductance
CT_L=[3 3 3 3 3 3 3 3]; %%--> 'CT' == 3  'CT' ==  4  %% Photosyntesis Typology for Plants
DSE_L =[0.6490    0.6490    0.6560    0.6490    0.6560    0.6490    0.6490    0.6490 ];  %% [kJ/mol] Activation Energy - Plant Dependent
Ha_L =[72    72    72    72    72    72    72    72]; %% [kJ / mol K]  entropy factor - Plant Dependent
gmes_L=[Inf   Inf   Inf   Inf   Inf   Inf   Inf   Inf]; %% [mol CO2 / s m^2 ];  mesophyll conductance
rjv_L=[2.1000    2.1000    2.4000    2.1000    2.2000    2.1000    2.1000    2.1000]; %%% Scaling Jmax - Vmax  [umol electrons / umolCO2 ]
%%%

%%%
FI_L=FI_L(II); Do_L=Do_L(II); a1_L=a1_L(II); go_L=go_L(II);
CT_L=CT_L(II); DSE_L=DSE_L(II); Ha_L=Ha_L(II); gmes_L=gmes_L(II);
rjv_L=rjv_L(II);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Psi_sto_00_H = [ -0.7000   -0.5000   -0.5000   -0.7000   -0.5000   -0.7000   -0.7000   -0.7000]; %% [MPa]  Water Potential at 2% loss conductivity
Psi_sto_50_H = [-2.5000   -2.0000   -2.0000   -2.5000   -2.0000   -2.5000   -2.5000   -2.5000] ;%% [MPa]  Water Potential at 50% loss conductivity
%%% Leaf
PsiL00_H = -1-[ -1.4000   -1.4000   -1.0000   -1.4000   -1.4000   -1.4000   -1.4000   -1.4000 ] ;%%[MPa]  Water Potential at 50% loss conductivity
PsiL50_H =  [-5.5000   -4.000   -2.0000   -5.5000   -3.5000   -5.5000   -5.5000   -6.5000]; %% [MPa]  Water Potential at 2% loss conductivity
Kleaf_max_H = [10 10 10 10 10 10 10 10 ] ; %%  %%%  [mmolH20 m^2 leaf s /MPa]
Cl_H  = [ 1200        1200        1200        1200        1200        1200        1200        1200];  %%%  [500 - 3000]%  Leaf capacitance [mmolH20 / m^2 leaf MPa]
%%% Xylem
Axyl_H = [ 15 15 15 15 15 15 15 15] ; %% [cm^2 stem /m^2 PFT]
Kx_max_H = [80000 80000 80000 80000 80000 80000 80000 80000];  %%5550-555550 [mmolH20 /m s MPa]  Xylem Conductivity specific for water;
PsiX50_H = [-9    -9    -9    -9    -9    -9    -9    -9]; %%[MPa]  Water Potential at 50% loss conductivity
Cx_H= [-10 -10 -10 -10 -10 -10 -10 -10]; %%%% [kg / m^3 sapwood MPa]
%%------------------------
%
Psi_sto_00_L = [ -0.7000   -0.5000   -0.5000   -0.7000   -0.5000   -0.7000   -0.7000   -0.7000]; %% [MPa]  Water Potential at 2% loss conductivity
Psi_sto_50_L = [-2.5000   -2.0000   -2.0000   -2.5000   -2.0000   -2.5000   -2.5000   -2.5000] ;%% [MPa]  Water Potential at 50% loss conductivity
%%% Leaf
PsiL00_L = -1-[ -1.4000   -1.4000   -1.0000   -1.4000   -1.4000   -1.4000   -1.4000   -1.4000 ] ;%%[MPa]  Water Potential at 50% loss conductivity
PsiL50_L =  [-5.5000   -4.000   -2.0000   -5.5000   -3.5000   -5.5000   -5.5000   -6.5000]; %% [MPa]  Water Potential at 2% loss conductivity
Kleaf_max_L = [10 10 10 10 10 10 10 10 ] ; %%  %%%  [mmolH20 m^2 leaf s /MPa]
Cl_L = [ 1200        1200        1200        1200        1200        1200        1200        1200];  %%%  [500 - 3000]%  Leaf capacitance [mmolH20 / m^2 leaf MPa]
%%% Xylem
Axyl_L = [ 15 15 15 15 15 15 15 15] ; %% [cm^2 stem /m^2 PFT]
Kx_max_L = [80000 80000 80000 80000 80000 80000 80000 80000];  %%5550-555550 [mmolH20 /m s MPa]  Xylem Conductivity specific for water;
PsiX50_L = [-9    -9    -9    -9    -9    -9    -9    -9]; %%[MPa]  Water Potential at 50% loss conductivity
Cx_L= [-10 -10 -10 -10 -10 -10 -10 -10]; %%%% [kg / m^3 sapwood MPa]
%%%
Psi_sto_50_H =Psi_sto_50_H(II);  Psi_sto_00_H =Psi_sto_00_H(II); 
PsiL00_H = PsiL00_H(II); PsiL50_H=PsiL50_H(II);  Kleaf_max_H=Kleaf_max_H(II); 
Cl_H=Cl_H(II); Axyl_H=Axyl_H(II); Kx_max_H=Kx_max_H(II); PsiX50_H=PsiX50_H(II); Cx_H=Cx_H(II); 
Psi_sto_50_L =Psi_sto_50_L(II);  Psi_sto_00_L =Psi_sto_00_L(II); 
PsiL00_L = PsiL00_L(II); PsiL50_L=PsiL50_L(II);  Kleaf_max_L=Kleaf_max_L(II); 
Cl_L=Cl_L(II); Axyl_L=Axyl_L(II); Kx_max_L=Kx_max_L(II); PsiX50_L=PsiX50_L(II); Cx_L=Cx_L(II); 

%%%%%%%%%%%%%%%% Root Parameters
[RfH_Zs,RfL_Zs]=Root_Fraction_General(Zs,CASE_ROOT,ZR95_H,ZR50_H,ZR95_L,ZR50_L,ZRmax_H,ZRmax_L); 

%%%% Growth Parameters 
PsiG50_H= [-0.7000   -0.7000   -0.7000   -0.7000   -0.7000   -0.7000   -0.7000   -0.7000];  %%[MPa]
PsiG99_H= [-5.5000   -3.5000   -2.0000   -5.5000   -2.5000   -5.5000   -5.5000   -6.5000];  %%[MPa]
gcoef_H = [3.5000    3.5000    3.5000    3.5000    3.5000    3.5000    3.5000    3.5000]; % [gC/m2 day]
%%------  
PsiG50_L= [-0.7000   -0.7000   -0.7000   -0.7000   -0.7000   -0.7000   -0.7000   -0.7000];  %%[MPa]
PsiG99_L= [-5.5000   -3.5000   -2.0000   -5.5000   -2.5000   -5.5000   -5.5000   -6.5000];  %%[MPa]
gcoef_L = [3.5000    3.5000    3.5000    3.5000    3.5000    3.5000    3.5000    3.5000]; % [gC/m2 day]
%%%%%
PsiG50_H=PsiG50_H(II); PsiG99_H=PsiG99_H(II); gcoef_H=gcoef_H(II); 
PsiG50_L=PsiG50_L(II); PsiG99_L=PsiG99_L(II); gcoef_L=gcoef_L(II); 

OPT_PROP_H =[1 5 13 5 13 5 5 5 ];
OPT_PROP_L =[1 5 13 5 13 5 5 5 ];
OPT_PROP_H=OPT_PROP_H(II);
OPT_PROP_L=OPT_PROP_L(II);
for i=1:cc
    %%%%%%%% Vegetation Optical Parameter
    [PFT_opt_H(i)]=Veg_Optical_Parameter(OPT_PROP_H(i));
    [PFT_opt_L(i)]=Veg_Optical_Parameter(OPT_PROP_L(i));
end

OM_H=[1 1 1 1 1 1 1 1];
OM_L=[1 1 1 1 1 1 1 1];
Sllit = 2 ; %%% Litter Specific Leaf area [m2 Litter / kg DM]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% VEGETATION PART %%%%%%
%%% HIGH VEGETATION
%%%%%%%%%%%%%%%%%%%%%%%%%
Sl_H = [0.0150    0.0140    0.0300    0.0150    0.0300    0.0130    0.0140    0.0280 ]; % 0.05 -0.005 [m^2 gC] specific leaf area of  biomass [m^2 /gC]
Nl_H= [48    40    23    48    35    35    40    35]; %[gC/gN ] Leaf Carbon-Nitrogen ratio
r_H = [0.0400    0.0200    0.0450    0.0400    0.0300    0.0400    0.0400    0.0200];  %% [0.066 -0.011]respiration rate at 10° [gC/gN d ]
gR_H= [0.2500    0.2500    0.2500    0.2500    0.2200    0.2500    0.2500    0.2200]; % [0.22 - 0.28] growth respiration  [] -- [Rg/(GPP-Rm)]
aSE_H= [0     0     2     0     5     0     0     1]; %%% Plant Type -- 1 Seasonal Plant --  0 Evergreen  -- 2 Grass species -- 3 Crops
dd_max_H= [0.0001    0.0008    0.0200    0.0001    0.0050    0.0001    0.0001    0.0000]; %%%0.005  [1/d]  0.0250 -- 0.005-0.025 death maximum for drought
dc_C_H =  [0.2137    0.0274         0    0.2137    0.0274    0.2137    0.2137    0.2137]; %% [1/ d°C] -- [Factor of increasing mortality]
Tcold_H = [-5    -5   -10    -5   -10    -5    -5     5 ]; %% [°C] Cold Leaf Shed
drn_H=  [0.0020    0.0020    0.0022    0.0020    0.0020    0.0020    0.0020    0.0020]; %% turnover root  [1/d]
dsn_H= [0.0027    0.0027    0.0027    0.0027    0.0014    0.0027    0.0027    0.0027]; % normal transfer rate sapwood [1/d]
age_cr_H= [650   700    80   650    90   650   650   110]; %% [day] Critical Leaf Age
Bfac_lo_H= [ 0.9800    0.9000    0.9600    0.9800    0.9800    0.9800    0.9800    0.9000]; %% Leaf Onset Water Stress
Bfac_ls_H= [0.6000       NaN    0.6500    0.6000       NaN    0.6000    0.6000       NaN ]; %% Leaf Shed Water Stress [0-1]
Tlo_H = [10.3000    4.0000    5.0000   10.3000   10.0000   10.3000   10.3000   13.0000]; %% Mean Temperature for Leaf onset
Tls_H = [NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN    10]; %% Mean Temperature for Leaf Shed
PAR_th_H= [NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN];
dmg_H= [20    30    30    20    20    20    20    30]; %%% Tree 30 Grasses Day of Max Growth
LAI_min_H = [0.0010    0.0010    0.0001    0.0010    0.0100    0.0010    0.0010    0.0010];
Trr_H = [0.1500    1.5000    4.0000    0.1500    4.0000    0.2500    0.1500    3.0000]; %% Translocation rate [gC /m^2 d]
mjDay_H = [180   366  -320   180  -280   180   180   230]; %% Maximum Julian day for leaf onset
LDay_min_H =[12.7000    9.5000    7.5000   12.7000    5.5000   12.7000   12.7000   11.5000]; %% Minimum Day duration for leaf onset
LtR_H = [0.9000    0.7500    0.3500    1.2000    0.8000    0.7500    0.7500    1.5000]; %%% Leaf to Root ratio maximum
Mf_H= [0.0125    0.0125    0.0125    0.0125    0.0125    0.0125    0.0125    0.0125 ]; %% fruit maturation turnover [1/d]
Wm_H= [0 0 0 0 0 0 0 0 ] ; % wood turnover coefficient [1/d]
eps_ac_H = [0.4000    0.4000    0.2000    0.4000    0.2000    0.4000    0.4000    0.3000]; %% Allocation to reserve parameter [0-1]
LDay_cr_H = [11.8000    9.0000    5.0000   11.8000    8.0000   11.8000   11.8000   12.0000]; %%%  Threshold for senescence day light [h]
Klf_H =[0.0250    0.0250    0.0250    0.0250    0.0250    0.0250    0.0250    0.0250 ]; %% Dead Leaves fall turnover [1/d]
fab_H = [0.7400    0.7400         0    0.7400         0    0.7400    0.7400    0.7400]; %% fraction above-ground sapwood and reserve
fbe_H = [0.2600    0.2600    1.0000    0.2600    1.0000    0.2600    0.2600    0.2600]; %% fraction below-ground sapwood and reserve
ff_r_H= [0.1000    0.1000    0.1000    0.1000    0.1000    0.1000    0.1000    0.1000]; %% 
%%%%

soCrop_H = [0.02 0.02 0.02 0.02 0.02 0.02 0.02 0.02]; 
Sl_emecrop_H = [0.02 0.02 0.02 0.02 0.02 0.02 0.02 0.02];
MHcrop_H = [2 2 2 2 2 2 2 2];
soCrop_H=soCrop_H(II); Sl_emecrop_H=Sl_emecrop_H(II); MHcrop_H=MHcrop_H(II);

Sl_H =Sl_H(II); Nl_H=Nl_H(II); 
r_H=r_H(II); gR_H=gR_H(II); aSE_H=aSE_H(II); dd_max_H=dd_max_H(II);
dc_C_H=dc_C_H(II); Tcold_H=Tcold_H(II); drn_H=drn_H(II);
dsn_H=dsn_H(II);  age_cr_H=age_cr_H(II);
Bfac_lo_H=Bfac_lo_H(II); Bfac_ls_H=Bfac_ls_H(II);
Tlo_H = Tlo_H(II);  Tls_H=Tls_H(II);
dmg_H = dmg_H(II); LAI_min_H=LAI_min_H(II);
Trr_H = Trr_H(II);  mjDay_H=mjDay_H(II);
LDay_min_H= LDay_min_H(II); LtR_H =LtR_H(II);
Mf_H= Mf_H(II);  Wm_H= Wm_H(II);  eps_ac_H = eps_ac_H(II);
LDay_cr_H = LDay_cr_H(II);  Klf_H = Klf_H(II);
fab_H = fab_H(II); fbe_H = fbe_H(II); ff_r_H = ff_r_H(II);

for i=1:cc
    [Stoich_H(i)]=Veg_Stoichiometric_Parameter(Nl_H(i));
    [ParEx_H(i)]=Exudation_Parameter(0);
    [Mpar_H(i)]=Vegetation_Management_Parameter;
    if ANSWER == 5
        Mpar_H(i).Crop_B = [40 5];
        Mpar_H(i).Date_sowing = datenum(datetime([1970:2025], 11, 10));
        Mpar_H(i).Date_harvesting = datenum(datetime([1970:2025], 5, 7));
        Mpar_H(i).Crop_crown = ones(length(Mpar_H(1).Date_sowing ), 1);

    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% LOW VEGETATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sl_L = [0.0150    0.0140    0.0300    0.0150    0.0300    0.0130    0.0140    0.0280 ]; % 0.05 -0.005 [m^2 gC] specific leaf area of  biomass [m^2 /gC]
Nl_L = [48    40    23    48    35    35    40    35]; %[gC/gN ] Leaf Carbon-Nitrogen ratio
r_L = [0.0400    0.0200    0.0450    0.0400    0.0300    0.0400    0.0400    0.0200];  %% [0.066 -0.011]respiration rate at 10° [gC/gN d ]
gR_L= [0.2500    0.2500    0.2500    0.2500    0.2200    0.2500    0.2500    0.2200]; % [0.22 - 0.28] growth respiration  [] -- [Rg/(GPP-Rm)]
aSE_L= [0     0     2     0     5     0     0     1]; %%% Plant Type -- 1 Seasonal Plant --  0 Evergreen  -- 2 Grass species -- 3 Crops
dd_max_L= [0.0001    0.0008    0.0200    0.0001    0.0050    0.0001    0.0001    0.0000]; %%%0.005  [1/d]  0.0250 -- 0.005-0.025 death maximum for drought
dc_C_L =  [0.2137    0.0274         0    0.2137    0.0274    0.2137    0.2137    0.2137]; %% [1/ d°C] -- [Factor of increasing mortality]
Tcold_L = [-5    -5   -10    -5   -10    -5    -5     5 ]; %% [°C] Cold Leaf Shed
drn_L=  [0.0020    0.0020    0.0022    0.0020    0.0020    0.0020    0.0020    0.0020]; %% turnover root  [1/d]
dsn_L= [0.0027    0.0027    0.0027    0.0027    0.0014    0.0027    0.0027    0.0027]; % normal transfer rate sapwood [1/d]
age_cr_L= [650   700    80   650    90   650   650   110]; %% [day] Critical Leaf Age
Bfac_lo_L= [ 0.9800    0.9000    0.9600    0.9800    0.9800    0.9800    0.9800    0.9000]; %% Leaf Onset Water Stress
Bfac_ls_L= [0.6000       NaN    0.6500    0.6000       NaN    0.6000    0.6000       NaN ]; %% Leaf Shed Water Stress [0-1]
Tlo_L = [10.3000    4.0000    5.0000   10.3000   10.0000   10.3000   10.3000   13.0000]; %% Mean Temperature for Leaf onset
Tls_L = [NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN    10]; %% Mean Temperature for Leaf Shed
PAR_th_L= [NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN];
dmg_L= [20    30    30    20    20    20    20    30]; %%% Tree 30 Grasses Day of Max Growth
LAI_min_L = [0.0010    0.0010    0.0001    0.0010    0.0100    0.0010    0.0010    0.0010];
Trr_L= [0.1500    1.5000    4.0000    0.1500    4.0000    0.2500    0.1500    3.0000]; %% Translocation rate [gC /m^2 d]
mjDay_L = [180   366  -320   180  -280   180   180   230]; %% Maximum Julian day for leaf onset
LDay_min_L =[12.7000    9.5000    7.5000   12.7000    5.5000   12.7000   12.7000   11.5000]; %% Minimum Day duration for leaf onset
LtR_L = [0.9000    0.7500    0.3500    1.2000    0.8000    0.7500    0.7500    1.5000]; %%% Leaf to Root ratio maximum
Mf_L= [0.0125    0.0125    0.0125    0.0125    0.0125    0.0125    0.0125    0.0125 ]; %% fruit maturation turnover [1/d]
Wm_L= [0 0 0 0 0 0 0 0 ] ; % wood turnover coefficient [1/d]
eps_ac_L = [0.4000    0.4000    0.2000    0.4000    0.2000    0.4000    0.4000    0.3000]; %% Allocation to reserve parameter [0-1]
LDay_cr_L = [11.8000    9.0000    5.0000   11.8000    8.0000   11.8000   11.8000   12.0000]; %%%  Threshold for senescence day light [h]
Klf_L =[0.0250    0.0250    0.0250    0.0250    0.0250    0.0250    0.0250    0.0250 ]; %% Dead Leaves fall turnover [1/d]
fab_L = [0.7400    0.7400         0    0.7400         0    0.7400    0.7400    0.7400]; %% fraction above-ground sapwood and reserve
fbe_L = [0.2600    0.2600    1.0000    0.2600    1.0000    0.2600    0.2600    0.2600]; %% fraction below-ground sapwood and reserve
ff_r_L= [0.1000    0.1000    0.1000    0.1000    0.1000    0.1000    0.1000    0.1000]; %% 
%%%%
%%%

soCrop_L = [0.02 0.02 0.02 0.02 0.02 0.02 0.02 0.02]; 
Sl_emecrop_L = [0.02 0.02 0.02 0.02 0.02 0.02 0.02 0.02];
MHcrop_L = [2 2 2 2 2 2 2 2];
soCrop_L=soCrop_L(II); Sl_emecrop_L=Sl_emecrop_L(II); MHcrop_L=MHcrop_L(II);

Sl_L =Sl_L(II); Nl_L=Nl_L(II);
r_L=r_L(II); gR_L=gR_L(II); aSE_L=aSE_L(II); dd_max_L=dd_max_L(II);
dc_C_L=dc_C_L(II); Tcold_L=Tcold_L(II); drn_L=drn_L(II);
dsn_L=dsn_L(II);  age_cr_L=age_cr_L(II);
Bfac_lo_L=Bfac_lo_L(II); Bfac_ls_L=Bfac_ls_L(II);
Tlo_L = Tlo_L(II);  Tls_L=Tls_L(II);
dmg_L = dmg_L(II); LAI_min_L=LAI_min_L(II);
Trr_L = Trr_L(II);  mjDay_L=mjDay_L(II);
LDay_min_L= LDay_min_L(II); LtR_L =LtR_L(II);
Mf_L= Mf_L(II);  Wm_L= Wm_L(II);  eps_ac_L = eps_ac_L(II);
LDay_cr_L = LDay_cr_L(II);  Klf_L = Klf_L(II);
fab_L = fab_L(II); fbe_L = fbe_L(II); ff_r_L = ff_r_L(II);

for i=1:cc
    [Stoich_L(i)]=Veg_Stoichiometric_Parameter(Nl_L(i));
    [ParEx_L(i)]=Exudation_Parameter(0);
    [Mpar_L(i)]=Vegetation_Management_Parameter;
end

%%%%
Vmax_H = [ 40    60    90    40    95    50    50    50]; %
Vmax_L = [0     0     0     0     0     0     0     0 ]; %
Vmax_H =Vmax_H(II); Vmax_L =Vmax_L(II);
%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%
Restating_parameters;
if dbThick>5
    nst=md_max;
    k=(dbThick/5)^(1/(nst-1));
    Zs_deb = [0 5*k.^(1:nst-1)]; %% [mm]
    clear nst k
else
end

end