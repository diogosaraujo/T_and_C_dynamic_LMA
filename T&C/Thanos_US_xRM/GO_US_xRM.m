clear
addpath('Code')

%%
dt=3600; %%[s] %%% 
dth=1; %%[h]

ms = 13; %%% Soil Layer 
cc = 1; %% Crown area

SLA_ex = load('LMA_US_xRM.mat');

id_location = 'US_xRM';
load('Meteo_US_xRM_1985_2020.mat')
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
x1=1;
x2=length(Ta);
NN = x2-x1+1;% %%% time Step

%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%5
Date=Date(x1:x2);
Pr=Pr(x1:x2);
Ta=Ta(x1:x2);
Ws=Ws(x1:x2); ea=ea(x1:x2);  SAD1=SAD1(x1:x2);
SAD2=SAD2(x1:x2); SAB1=SAB1(x1:x2); Pre=Pre(x1:x2);
SAB2=SAB2(x1:x2); N=N(x1:x2); Tdew=Tdew(x1:x2);esat=esat(x1:x2);
PARB=PARB(x1:x2); PARD = PARD(x1:x2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%
% t_bef= -0.67; t_aft= 1.67;
%%%%%%%%%%%%%%%%%%%%
Ds=esat-ea; %% [Pa] Vapor Pressure Deficit
Ds(Ds<0)=0;

Oa= 210000;% Intercellular Partial Pressure Oxygen [umolO2/mol] -

Ws(Ws<=0)=0.01;
%%dt,Pr(i),Ta(i),Ws(i),ea(i),Pre(i),Rdir(i),Rdif(i),N(i),z,Tdew(i),esat(i),.
[YE,MO,DA,HO,MI,SE] = datevec(Date);
Datam(:,1) = YE; Datam(:,2)= MO; Datam(:,3)= DA; Datam(:,4)= HO;
clear YE MO DA HO MI SE

PARAM_IC = strcat('MOD_PARAM_',id_location);

MAIN_FRAME_SLA ;

rmpath('Code')

save(['RES_', id_location])

load('RES_US_xRM.mat')
GRAPH_MOD 