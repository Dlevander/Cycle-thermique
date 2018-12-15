function [ETA,DATEN,DATEX,DAT,MASSFLOW,COMBUSTION,FIG] = GT(P_e,options,display)
% GT Gas turbine modelisation
% GT(P_e,options,display) compute the thermodynamics states for a Gas
% turbine based on several inputs (given in OPTION) and based on a given 
% electricity production P_e. It returns the main results. It can as well
% plots graphs if input argument DISPLAY = true (<=> DISPLAY=1)
%
% INPUTS (some inputs can be dependent on others => only one of these 2 can
%         be activated)
% P_E = electrical power output target [kW]
% OPTIONS is a structure containing :
%   -options.k_mec [-] : Shaft losses 
%   -options.T_0   [C] : Reference temperature
%   -options.T_ext [C] : External temperature
%   -options.r     [-] : Comperssion ratio
%   -options.k_cc  [-] : Coefficient of pressure losses due to combustion
%                        chamber
%   -options.T_3   [C] : Temperature after combustion (before turbine)
%   -option.eta_PiC[-] : Intern polytropic efficiency (Rendement
%                        polytropique interne) for compression
%   -option.eta_PiT[-] : Intern polytropic efficiency (Rendement
%                        polytropique interne) for expansion
%DISPLAY = 1 or 0. If 1, then the code should plot graphics. If 0, then the
%          do not plot.
%
%OUPUTS : 
% ETA is a vector with :
%   -eta(1) : eta_cyclen, cycle energy efficiency
%   -eta(2) : eta_toten, overall energy efficiency
%   -eta(3) : eta_cyclex, cycle exegy efficiency
%   -eta(4) : eta_totex, overall exergie efficiency
%   -eta(5) : eta_rotex, compressor-turbine exergy efficiency
%   -eta(6) : eta_combex, Combustion exergy efficiency
%   FYI : eta(i) \in [0;1] [-]
% DATEN is a vector with : 
%   -daten(1) : perte_mec [kW]
%   -daten(2) : perte_ech [kW]
% DATEX is a vector with :
%   -datex(1) : perte_mec [kW]
%   -datex(2) : perte_rotex [kW]
%   -datex(3) : perte_combex [kW]
%   -datex(4) : perte_echex  [kW]
% DAT is a matrix containing :
% dat = {T_1       , T_2       , T_3       , T_4; [°C]
%        p_1       , p_2       , p_3       , p_4; [bar]
%        h_1       , h_2       , h_3       , h_4; [kJ/kg]
%        s_1       , s_2       , s_3       , s_4; [kJ/kg/K]
%        e_1       , e_2       , e_3       , e_4;};[kJ/kg]
% MASSFLOW is a vector containing : 
%   -massflow(1) = m_a, air massflow [kg/s]
%   -massflow(2) = m_c, combustible massflow [kg/s] 
%   -massflow(3) = m_f, exhaust gas massflow [kg/s]
% 
% COMBUSTION is a structure with :
%   -combustion.LHV    : the Lower Heat Value of the fuel [kJ/kg]
%   -combustion.e_c    : the combuistible exergie         [kJ/kg]
%   -combustion.lambda : the air excess                   [-]
%   -combustion.Cp_g   : heat capacity of exhaust gas at 400 K [kJ/kg/K]
%   -combustion.fum  : is a vector of the exhaust gas composition :
%       -fum(1) = m_O2f  : massflow of O2 in exhaust gas [kg/s]
%       -fum(2) = m_N2f  : massflow of N2 in exhaust gas [kg/s]
%       -fum(3) = m_CO2f : massflow of CO2 in exhaust gas [kg/s]
%       -fum(4) = m_H2Of : massflow of H2O in exhaust gas [kg/s] 
%
% FIG is a vector of all the figure you plot. Before each figure, define a
% figure environment such as:  
%  "FIG(1) = figure;
%  plot(x,y1);
%  [...]
%   FIG(2) = figure;
%  plot(x,y2);
%  [...]"
%  Your vector FIG will contain all the figure plot during the run of this
%  code (whatever the size of FIG).
%


%% Your Work

if nargin<3
    display=1;
   if nargin<2
       options=struct();
       if nargin<1
           P_e=230e3; %[kW]
       end
   end
end

if isfield(options,'T_0')
    T_0 = options.T_0;
else
    T_0 = 15;%C   
end

if isfield(options,'k_mec')
    k_mec = options.k_mec;
else
    k_mec = 0.015;   
end

if isfield(options,'T_ext')
    T_ext = options.T_ext;
else
    T_ext = 15;   %C
end

if isfield(options,'r')
    r = options.r;
else
    r = 18;   
end

if isfield(options,'k_cc')
    k_cc = options.k_cc;
else
    k_cc = 0.95;   
end

if isfield(options,'T_3')
    T_3 = options.T_3;
else
    T_3 = 1400;   %[C]
end

if isfield(options,'eta_PiC')
    eta_PiC = options.eta_PiC;
else
    eta_PiC = 0.9;   
end

if isfield(options,'eta_PiT')
    eta_PiT = options.eta_PiT;
else
    eta_PiT = 0.9;   
end

x = 0;
y = 4; % Hypothese : CH4


%%%%%%%%%%%%%%% OUTPUT %%%%%%%%%%%%%%%%%%%
ETA = zeros(6,1);

DATEN = zeros(2,1);

DATEX = zeros(4,1);

DAT= zeros(5,4);

MASSFLOW = zeros(3,1); 

FIG = 0; %A MODIFIER

%%%%% Calcul des états %%%%%%
R_air = 287.1; %J/kg.K 
x_a_O2 = 0.21*32/28.96; %fraction massique de O2 dans l'air
x_a_N2 = 0.79*28/28.96; %fraction massique de N2 dans l'air
% %calcul point 1 : air atmosphérique

p_1 = 1,01325;  %bar
T_1 = T_ext ;
Cp_air_27 = 1000*(x_a_N2*janaf('c','N2',300)+x_a_O2*janaf('c','O2',300)) ; %J/kg*K
h_1 = T_ext * Cp_air_27 ;
s_1 = Cp_air_27*log((T_ext+273.15)/273.15); %J/kg*K
e_1 = 0; % point de reference

%%calcul point 2 : apres la pompe

p_2 = r*p_1 ;
%T_2=429;
T_2 = transf_poly('compression',T_ext,r,eta_PiC,R_air,0,0,0,0); %renvoie T en sortie de compresseur en C, pas besoin des compo fumee pour 1 compr
TK_2=T_2+273.15;
Cp_12 = 1000*(x_a_N2*mean(janaf('c','N2',linspace(300,TK_2)))+x_a_O2*mean(janaf('c','O2',linspace(300,TK_2)))) ; %J/kg*K faire cp moyen entre T1 et T2 ?
h_2 = h_1 + Cp_12*(T_2-T_1);
s_2 = s_1 + (1-eta_PiC)*Cp_12*log((TK_2)/(T_1+273.15)); %eq 3.15
%h_2=443400;
%s_2=142;
e_2 = (h_2-h_1) - 273.15*(s_2-s_1);

%%calcul point 3 : après la combustion

T_3 = T_3 ; % really ?
p_3 = p_2*k_cc; %pertes de charges dans chambre combustion
TK_3=T_3+273.15;
[x_N2 x_O2 x_CO2 x_H2O R_fum lambda ma1 LHV e_c] = combustion(x,y,T_2,T_3,0);
Cp_23 = CP(x_O2,x_CO2,x_H2O,x_N2,[TK_2,TK_3])*1000;
%CpMoy_23 =CPmoy(x_O2,x_CO2,x_H2O,x_N2,[TK_2,TK_3])*1000;
%Cp_23_2 =(x_N2*mean(janaf('c','N2',linspace(TK_2,TK_3))) + x_O2*mean(janaf('c','O2',linspace(TK_2,TK_3))) + x_CO2*mean(janaf('c','CO2',linspace(TK_2,TK_3))) + x_H2O*mean(janaf('c','H2O',linspace(TK_2,TK_3))))*1000; %faire cp moyen entre T2 et T3 ?
h_3 = CP(x_O2,x_CO2,x_H2O,x_N2,[300,TK_3])*1000*(T_3-15)+h_1;
%h_3 =Cp_23*(T_3-T_2) + h_2; %un peu trop petit, h2 et t3 et t2 sont bons 
s_3 = s_2+ Cp_23*log(TK_3/TK_2) - R_fum*log(p_3/p_2);
e_3 = (h_3-h_1) - 273.15*(s_3-s_1);

%%Calcul point 4 : aprs la turbine
p_4 = p_1 ;% atm
T_4 = transf_poly('detente',T_3,p_4/p_3,eta_PiT,R_fum,x_CO2,x_H2O,x_O2,x_N2);
TK_4=T_4+273.15;
Cp_34 = CP(x_O2,x_CO2,x_H2O,x_N2,[TK_3,TK_4])*1000;
%Cp_34_2 =(x_N2*mean(janaf('c','N2',linspace(TK_3,TK_4))) + x_O2*mean(janaf('c','O2',linspace(TK_3,TK_4))) + x_CO2*mean(janaf('c','CO2',linspace(TK_3,TK_4))) + x_H2O*mean(janaf('c','H2O',linspace(TK_3,TK_4))))*1000;
h_4 = h_3 + Cp_34*(T_4-T_3);
s_4 = s_3 - Cp_34*log(TK_4/TK_3)* ((1-eta_PiT)/eta_PiT);%eq3.16
e_4 = (h_4-h_1) - 273.15*(s_4-s_1);

%%Remplissage etats

DAT(:,1) = [T_1 p_1 h_1 s_1 e_1]';
DAT(:,2) = [T_2 p_2 h_2 s_2 e_2]';
DAT(:,3) = [T_3 p_3 h_3 s_3 e_3]';
DAT(:,4) = [T_4 p_4 h_4 s_4 e_4]';
%%Rendements énergétiques %%
P_e=P_e*10^3;
eta_cyclen = 1-(((1+1/(lambda*ma1))*h_4-h_1)/((1+1/(lambda*ma1))*h_3-h_2)); %eq3.12
%eta_toten = eta_mec*eta_cyclen; %eta_mec=P_e/(P_T-P_C)
W_T = h_3-h_4; %J/kg
W_C = h_2-h_1;

m_a = P_e/(W_T+W_T/(lambda*ma1)-W_C-k_mec*W_T-k_mec*W_T/(lambda*ma1)-W_C*k_mec); %eq3.29,3.1,3.7 et 3.8
m_c = m_a/(lambda*ma1);
m_g = m_a+m_c;
eta_toten = P_e/(LHV*1000*m_c);
eta_mec = eta_toten/eta_cyclen;

% m_c=P_e/(LHV*1000*eta_toten);%debit comb
% m_a=lambda*ma1*m_c;%debit air
% m_g=m_a+m_c;%debit fumees

%% rendement exergétique

eta_cyclex=(m_g*(h_3-h_4)-m_a*(h_2-h_1))/(m_g*e_3-m_a*e_2);
eta_rotex=(m_g*(h_3-h_4)-m_a*(h_2-h_1))/(m_g*(e_3-e_4)-m_a*(e_2-e_1));
%eta_combex=(m_g*e_3-m_a*e_2)/(m_g*h_3-m_a*h_2);%eq3.36 sans prendre en compte f
eta_combex = (m_g*e_3-m_a*e_2)/(m_c*e_c*1000); %eq3.36 avec e_c = 52215 pour CH4
eta_totex=eta_cyclex*eta_mec*eta_combex; %p124

%% Pertes
Pm=P_e/k_mec;%puissance motrice
Pprim=m_c*LHV*1000;%puissance totale
pertes_mec=k_mec*(W_T*m_g+W_C*m_a);%eq3.29
puissance_effect=Pprim*eta_totex;
pertesTotex= Pprim*(1-eta_totex);
pertes_cyclex= (m_g*e_3-m_a*e_2)-Pm;
pertes_combu = (1-eta_combex)*e_c*m_c*1000;
pertes_cycle=Pprim*(1-eta_cyclex);
pertes_echapex=(m_g*e_3-m_a*e_2)-(m_g*(e_3-e_4)-m_a*(e_2-e_1));
pertes_ExC =  ((h_2-h_1) - (e_2-e_1))*m_a; %Watts/sec
pertes_ExT =  ((e_3-e_4) - (h_3-h_4))*m_g;
pertes_rotex = pertes_ExC + pertes_ExT;
pertes_echap = h_4*m_g;

%% Remplissage vecteurs
MASSFLOW(1) = m_a;
MASSFLOW(2) = m_c;
MASSFLOW(3) = m_g;

ETA(1) = eta_cyclen;
ETA(2) = eta_toten;
ETA(3) = eta_cyclex;
ETA(4) = eta_totex;
ETA(5) = eta_rotex;
ETA(6) = eta_combex;

DATEN(1) = pertes_mec/1000;
DATEN(2) = pertes_echap/1000; 

DATEX(1) = pertes_mec/1000;
DATEX(2) = pertes_rotex/1000;
DATEX(3) = pertes_combu/1000;
DATEX(4) = pertes_echapex/1000 ; 

COMBUSTION.LHV = LHV ;
COMBUSTION.e_c = e_c ;
COMBUSTION.lambda = lambda ;
COMBUSTION.Cp_g = CP(x_O2,x_CO2,x_H2O,x_N2,400+273.15);
%COMBUSTION.Cp_g  = (x_N2*janaf('c','N2',400) + x_O2*janaf('c','O2',400) + x_CO2*janaf('c','CO2',400) + x_H2O*janaf('c','H2O',400))*1000;
COMBUSTION.fum(1) =  x_O2*m_g ;
COMBUSTION.fum(2) =  x_N2*m_g ;
COMBUSTION.fum(3) =  x_CO2*m_g;
COMBUSTION.fum(4) =  x_H2O*m_g;
for i=1:4 %DAT en kJ/kg
    DAT(3,i)=DAT(3,i)/1000;
    DAT(4,i)=DAT(4,i)/1000;
    DAT(5,i)=DAT(5,i)/1000;
end
if display ==1
    %FIG(1)=plot_GT(DAT,eta_PiT,eta_PiC,R_air);
end

end
