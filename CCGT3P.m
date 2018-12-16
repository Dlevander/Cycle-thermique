function [ETA MASSFLOW] = CCGT3P(P_eg,options,display)
% CCGT3P is a Combine cycle Gas Turbine with 2 pressure level
% CCGT3P(P_e,options,display) compute the thermodynamics states for a CCGT
% with 3 pressure level (cfr p166 english reference book) including
% combustion, exchanger and cycles. This is done based on several inputs
% (given in OPTION) and based on a given electricity production P_e.
% It returns the main results. It can as well plots graphs if input 
% argument DISPLAY = true (<=> DISPLAY=1)
%
% INPUTS (some inputs can be dependent on others => only one of these 2 can
%         be activated)
% P_EG = electrical power output target for gas turbine [kW]
% OPTIONS is a structure containing :
%   -options.T0       [°C] : Reference temperature
%   -options.T_ext    [°C] : External temperature
%   -options.T_STmax  [°C] : maximum temperature on ST cycle
%   -options.eta_mec  [-] : mecanic efficiency of shafts bearings
%   -options.pdrum   [bar]: Drum pressure
%   -options.pmid    [bar]: Intermediary pressure level
%   -options.x7       [-] : Vapor ratio [gaseous/liquid] (titre)
%   -option.eta_SiC   [-] : Isotrenpic efficiency for compression
%   -option.eta_SiT   [-] : Isotrenpic efficiency for compression
%   -options.GT    [struct] : options for Gas turbine (see GT function) 
% DISPLAY = 1 or 0. If 1, then the code should plot graphics. If 0, then 
%          do not plot.
%
%OUPUTS : 
% ETA is a vector with :
%   -eta(1)  : eta_STcyclen, cycle energy efficiency
%   -eta(2)  : eta_GTcyclen, cycle energy efficiency
%   -eta(3)  : eta_toten, overall energy efficiency
%   -eta(4)  : eta_STcyclex, cycle exegy efficiency
%   -eta(5)  : eta_GTcyclex, cycle exegy efficiency
%   -eta(6)  : eta_totex, overall exergie efficiency
%   -eta(7)  : eta_gen, Steam generator energy efficiency
%   -eta(8)  : eta_gex, Steam generator exergy efficiency
%   -eta(9)  : eta_combex, Combustion exergy efficiency
%   -eta(10) : eta_chemex, Chimney exergy efficiency (losses)
%   -eta(11) : eta_transex, Heat exchanger overall exergy efficiency
%   FYI : eta(i) \in [0;1] [-]
% MASSFLOW is a vector containing : 
%   -massflow(1) [kg/s]: water massflow at high pressure turbine inlet
%   -massflow(2) [kg/s]: water massflow at medium pressure turbine inlet
%   -massflow(3) [kg/s]: water massflow at low pressure turbine inlet
%   -massflow(4) [kg/s]: air massflow at gas turbine inlet 
%   -massflow(5) [kg/s]: combustible massflow
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

%% Your work
if nargin<3
    display=1;
   if nargin<2
       options=struct();
       if nargin<1
           P_eg=100e3;%100MW
       end
   end
end

if isfield(options,'T_0')
    T_0 = options.T_0;
else
    T_0 = 15;%C   
end

if isfield(options,'T_ext')
    T_ext = options.T_ext;
else
    T_ext = 15;   %C
end

if isfield(options,'T_STmax')
    T_STmax = options.T_STmax;
else
    T_STmax = 565;   
end

if isfield(options,'eta_mec')
    eta_mec = options.eta_mec;
else
    eta_mec = 0.98;   
end


if isfield(options,'pdrum')
    pdrum = options.pdrum;
else
    pdrum = 4;   %bar
end

if isfield(options,'pmid')
    pmid = options.pmid;
else
    pmid = 28;   %bar
end


if isfield(options,'eta_SiC')
    eta_SiC = options.eta_SiC;
else
    eta_SiC = 0.9;   
end

if isfield(options,'eta_SiT')
    eta_SiT = options.eta_SiT;
else
    eta_SiT = 0.9;   
end

%% Outputs
ETA= zeros(11,1);
MASSFLOW=zeros(5,1);

%% Etats et rendements cycle GT
[ETA DATEN DATEX DAT MASSFLOW COMBUSTION] = GT(P_eg,options.GT,display);
eta_cycleng = ETA(1);
eta_cyclexg = ETA(3);
p_1g = DAT(2,1);
p_2g = DAT(2,2);
p_3g = DAT(2,3);
p_4g = DAT(2,4);
T_1g = DAT(1,1);
T_2g = DAT(1,2);
T_3g = DAT(1,3);
T_4g = DAT(1,4);
h_1g = DAT(3,1);
h_2g = DAT(3,2);
h_3g = DAT(3,3);
h_4g = DAT(3,4);
s_1g = DAT(4,1);
s_2g = DAT(4,2);
s_3g = DAT(4,3);
s_4g = DAT(4,4);
e_1g = DAT(5,1);
e_2g = DAT(5,2);
e_3g = DAT(5,3);
e_4g = DAT(5,4);
m_ag = MASSFLOW(1);
m_cg = MASSFLOW(2) ;
m_gg = MASSFLOW(3);
LHV = COMBUSTION.LHV;
e_c = COMBUSTION.e_c;
x_O2 = COMBUSTION.fum(1)/m_gg;
x_N2 = COMBUSTION.fum(2)/m_gg;
x_H2O = COMBUSTION.fum(4)/m_gg;
x_CO2 = COMBUSTION.fum(3)/m_gg;

eta_SiT = 0.9; %hypothèse
eta_SiC = 0.9;
T_STmax = 565;
x7 =0.95;
T_7 = 36; %comment faire autrement ? hypothèse
%% Etats ST
p_HP=122.8 ; %bar
p_IP=27.3;
p_LP=3.8;
P_es=153800; %kW

% Etat de référence
% T0=15;
% p0 = 1.01325;
% Cp_air_27 = 1000*(x_a_N2*janaf('c','N2',300)+x_a_O2*janaf('c','O2',300)) ; %J/kg*K
% h0 = T0 * Cp_air_27 ;
% s0 = Cp_air_27*log((T_ext+273.15)/273.15); %J/kg*K
% e0 = 0; % point de reference

%Etat 3 : entrée turbine HP
p_3 = p_HP;
T_3 = T_STmax;
h_3 = XSteam('h_pT',p_3,T_3);
s_3 = XSteam('s_pT',p_3,T_3);
x_3 = XSteam('x_ph',p_3,h_3);
e_3 = exergie(h_3,s_3);

%Etat 4 : sortie turbine HP (=entree turb IP ?)
p_4 = p_IP;
s_4s = s_3;
h_4s = XSteam('h_ps',p_4,s_4s);
x_4s = XSteam('x_ps',p_4,s_4s);

h_4 = h_3 - eta_SiT * (h_3-h_4s);
s_4 = XSteam('s_ph',p_4,h_4);
x_4 = XSteam('x_ps',p_4,s_4);
T_4 = XSteam('T_ph',p_4,h_4);
e_4 = exergie(h_4,s_4);

%Etat 5 : sortie du reheater entree turb IP
T_5 = T_STmax; %=T_3
p_5 = p_4 ; %resurchauffe isobare
h_5 = XSteam('h_pT',p_5,T_5);
s_5 = XSteam('s_pT',p_5,T_5);
e_5 = exergie(h_5,s_5);

%Etat 1 : sortie condenseur, liq sature
T_1 = T_7;
p_1 = XSteam('psat_T',T_1);
h_1 = XSteam('hL_T',T_1);
s_1 = XSteam('sL_T',T_1);
e_1 = exergie(h_1,s_1);
x_1 = 0;

%Etat 2 : après pompe alimentaire
p_2 = p_LP ; 
h_2s = XSteam('h_ps',p_2,s_1); %h dans le cas isentropique
h_2=h_1+(h_2s-h_1)/eta_SiC;%h en prenant compte rendement is
s_2=XSteam('s_ph',p_2,h_2);
T_2=XSteam('T_ph',p_2,h_2);
x_2=XSteam('x_ph',p_2,h_2);
e_2=exergie(h_2,s_2);

% Etat 81 entrée evaporisateur LP
p_81 = p_LP ;
T_81 = XSteam('Tsat_p',p_81);
h_81 = XSteam('hL_p',p_81);
s_81 = XSteam('sL_p',p_81);
x_81 = 0; %liquide saturé
e_81 = exergie(h_81,s_81);

%Etat 82 : sortie evap LP, vap saturée
p_82 = p_81;
T_82 = T_81;
h_82 = XSteam('hV_p',p_82);
s_82 = XSteam('sV_p',p_82);
x_82 = 1;
e_82 = exergie(h_82,s_82);

%Etat 91 : entree évap IP, liq saturee
p_91 = p_IP ;
T_91 = XSteam('Tsat_p',p_91);
h_91 = XSteam('hL_p',p_91);
s_91 = XSteam('sL_p',p_91);
x_91 = 0;
e_91 = exergie(h_91,s_91);

%Etat 92 : sortie evap IP, vap saturee
p_92 = p_91;
T_92 = T_91;
h_92 = XSteam('hV_p',p_92);
s_92 = XSteam('sV_p',p_92);
x_92 = 1;
e_92 = exergie(h_92,s_92);

%Etat 80 : sortie surchauffeur LP
p_80 = p_82;
T_80 = T_91;
h_80 = XSteam('h_pT',p_80,T_80);
s_80 = XSteam('s_pT',p_80,T_80);
e_80 = exergie(h_80,s_80);

%Etat 101 : entre evap HP
p_101 = p_HP ;
T_101 = XSteam('Tsat_p',p_101);
h_101 = XSteam('hL_p',p_101);
s_101 = XSteam('sL_p',p_101);
x_101 = 0; %liq sature
e_101 = exergie(h_101,s_101);

%Etat 102 : sortie evap HP
p_102 = p_101;
T_102 = T_101;
h_102 = XSteam('hV_p',p_102);
s_102 = XSteam('sV_p',p_102);
x_102 = 1; %vap sature
e_102 = exergie(h_102,s_102);

%Etat 90 : sortie surchauffeur IP
p_90 = p_92;
T_90 = T_101;
h_90 = XSteam('h_pT',p_90,T_90);
s_90 = XSteam('s_pT',p_90,T_90);
e_90 = exergie(h_90,s_90);

%Etat 6 : sortie turb IP et entree turbine LP
p_6 = p_LP; 
T_6 = T_102;
h_6 = XSteam('h_pT',p_6,T_6);
s_6 = XSteam('s_pT',p_6,T_6);
e_6 = exergie(h_6,s_6);

%Etat 7 : sortie turbine LP, entree condenseur
x_7 = x7;
T_7; %hypothese....si jms le temps, faire le truc juste en dessous
p_7 = XSteam('psat_T',T_7);
h_7 = XSteam('h_Tx',T_7,x_7);
s_7 = XSteam('s_ph',p_7,h_7);
e_7 = exergie(h_7,s_7);
% n=0;
% delta = 10;
% p_7 = 0.1;
% while (delta>0.01 && n<10)
%     p_7old = p_7;
%     s_7s = s_6;
%     h_7s = XSteam('h_ps',p_7old,s_7s);
%     h_7 = h_6 - eta_SiT * (h_6-h_7s);
%     s_7 = XSteam('s_ph',p_7old,h_7);
%     p_7 = XSteam('p_hs',h_7,s_7);
%     x_7 = XSteam('x_ph',p_7,h_7);
%     delta = abs(x_7-x7);
%     n=n+1;
% end


%% Calcul debit vapeur
T_pinch = 10; % hypothese, p153 livre en

T_P = T_pinch + T_101; %temp des fumees dans cheminee a chaque niv de pression, C
T_K = T_pinch + T_91;
T_M = T_pinch + T_81;
TK_P = T_pinch + T_101 +273.15; %temp des fumees dans cheminee a chaque niv de pression, K
TK_K = T_pinch + T_91 + 273.15;
TK_M = T_pinch + T_81 + 273.15;

Cp_fum_4gP = CP(x_O2,x_CO2,x_H2O,x_N2,[T_4g+273.15,TK_P]);
h_P = h_4g - Cp_fum_4gP*(T_4g-T_P); %kJ

Cp_fum_PK = CP(x_O2,x_CO2,x_H2O,x_N2,[TK_P,TK_K]);
h_K = h_P - (T_P-T_K)*Cp_fum_PK;

Cp_fum_KM = CP(x_O2,x_CO2,x_H2O,x_N2,[TK_K,TK_M]);
h_M = h_K -(T_K-T_M)*Cp_fum_KM;

solution = fsolve(@solver,[0,0,0]);
mv_HP = solution(1) ;
mv_IP =solution(2);
mv_LP = solution(3);
mv_tot = mv_HP+mv_LP+mv_IP;
    function Solve = solver(s)
        Solve(1) = s(1)*(h_3-h_101+h_5-h_4)+s(2)*(h_5-h_90)-m_gg*(h_4g-h_P);
        Solve(2) = s(1)*(h_101-h_91)+s(2)*(h_90-h_91)+s(3)*(h_6-h_80)-m_gg*(h_P-h_K);
        Solve(3) =s(1)*(h_91-h_81)+s(2)*(h_91-h_81)+s(3)*(h_80-h_81)-m_gg*(h_K-h_M);
    end
%% Etat 5g : sortie de cheminee
h_5g = h_M - mv_tot*(h_81-h_2)/m_gg;
TK_5g_approx = TK_M - 50; 
Cp_fum_M5g = CP(x_O2,x_CO2,x_H2O,x_N2,[TK_M,TK_5g_approx]); %J/kg*C
deltah = h_M - h_5g; %positif
T_5g = T_M - deltah/Cp_fum_M5g;
p_5g = p_4g;
Cp_fum_45 = CP(x_O2,x_CO2,x_H2O,x_N2,[T_4g + 273.15, T_5g + 273.15]);
s_5g = s_4g + Cp_fum_45*log((T_5g+273.15)/(T_4g +273.15));
e_5g=(h_5g-h_1g) - 288.15*(s_5g-s_1g);

%% Rendements 

%% Gas cycle %%%%%
W_TG = h_3g-h_4g;
W_CG = h_2g-h_1g;
P_T_G = W_TG*m_gg;
P_C_G = W_CG*m_ag;
W_mcyG = (m_gg/m_ag)*W_TG - W_CG;
P_mcyG = W_TG*m_gg - W_CG*m_ag;

P_primv_tot = LHV*m_cg; %kJ/kg

eta_mecG = P_eg/P_mcyG;
%P_fmecG = (1-k_mec)*(P_T_G + P_C_G);
pertesEn_mecG = P_mcyG - P_eg;

eta_totenG = 1 - ((m_gg*h_4g) - (m_ag*h_1g)) / (m_cg*LHV);
pertesEn_tot = (m_cg*LHV)* (1-eta_totenG);
pertesEx_C_G =  ((h_2g - h_1g) - (e_2g - e_1g))*m_ag; %Watts/sec
pertesEx_T_G =  ((e_3g - e_4g) - (h_3g- h_4g))*m_gg;
pertesEx_rotG = pertesEx_C_G + pertesEx_T_G; 
fluxEx_T_G =  ((e_3g - e_4g))*m_gg;
fluxEx_C_G =  ((e_2g - e_1g))*m_ag;


%Q_combTGV = Q_combG + Q_combR;
Q_HRSG = (h_4g - h_5g)*m_gg; % kJ/sec
Q_exh = h_5g*m_gg - h_1g*m_ag; %kJ/sec

delta_e_fum = (e_4g - e_5g) * m_gg; %kJ
delta_e_vap = (e_3 - e_81)*mv_HP + (e_5-e_4)*mv_HP + (e_5 - e_90)*(mv_IP)...
    + (e_90 - e_81)*mv_IP ...
    + (e_6 - e_81)*mv_LP + (e_81 - e_2)*mv_tot;
delta_h_vap = (h_3 - h_81)*mv_HP + (h_5-h_4)*mv_HP + (h_5 - h_90)*(mv_IP)...
    + (h_90 - h_81)*mv_IP ...
    + (h_6 - h_81)*mv_LP + (h_81 - h_2)*mv_tot;
zero_bilanH = Q_HRSG - delta_h_vap; % donne 3.4 - 3.2

eta_transex = delta_e_vap/ delta_e_fum;
pertesEx_trans = delta_e_fum - delta_e_vap;
eta_chemex = (e_4g - e_5g)/(e_4g-0.04);
eta_totexG = P_eg / (m_cg*e_c); % not including reheating

%% Vapour cycle %%%%%

W_T_HP = h_3-h_4; %en kJ/kg
P_T_HP = W_T_HP * mv_HP; %kJ/sec =kW
pertesEx_T_HP = - mv_HP * ( h_3-h_4 - e_3 + e_4);
fluxEx_T_HP = mv_HP * ( e_3-e_4);

W_T_IP = h_5-h_6;
P_T_IP = W_T_IP * (mv_IP + mv_HP);
pertesEx_T_IP = - (mv_IP + mv_HP) * ( h_5-h_6 - e_5 + e_6 );
fluxEx_T_IP = (mv_IP + mv_HP) * ( e_5-e_6 );

W_T_LP = h_6 - h_7;
P_T_LP = W_T_LP * mv_tot;
pertesEx_T_LP = - mv_tot * ( h_6 - h_7 - e_6 + e_7);
fluxEx_T_LP = mv_tot * ( e_6-e_7);

P_T_totV= P_T_HP + P_T_IP + P_T_LP;
pertesEx_T = pertesEx_T_LP + pertesEx_T_IP + pertesEx_T_HP;
fluxEx_T_V = fluxEx_T_LP + fluxEx_T_IP + fluxEx_T_HP;

W_C_IP = h_91 - h_81;
P_C_IP = W_C_IP * (mv_HP + mv_IP);
pertesEx_C_IP = - (mv_HP + mv_IP)* (e_91 - e_81 - h_91 + h_81 );
fluxEx_C_IP = (mv_HP + mv_IP)* (e_91 - e_81 );

W_C_HP = h_101 - h_91;
P_C_HP = W_C_HP * mv_HP;
pertesEx_C_HP = - mv_HP * (e_101 - e_91 - h_101 + h_91 );
fluxEx_C_HP = mv_HP * (e_101 - e_91 );

W_C_Pe = h_2 - h_1;
P_C_Pe = W_C_Pe * mv_tot;
pertesEx_C_Pe = -  mv_tot* (e_2 - e_1 - h_2 + h_1 );
fluxEx_C_Pe = mv_tot* (e_2 - e_1 );

P_C_totV = P_C_IP + P_C_HP + P_C_Pe;
pertesEx_C = pertesEx_C_Pe + pertesEx_C_HP + pertesEx_C_IP;
fluxEx_C_V = fluxEx_C_Pe + fluxEx_C_HP + fluxEx_C_IP;

P_mcyV = (P_T_totV - P_C_totV); %kJ/sec 
P_e_V = P_mcyV * eta_mec; %kW 166MW

P_primTGV= m_cg*LHV;
fluxEx_prim = m_cg*e_c; %kJ/kg
PmG=P_eg/eta_mec;
PmV=P_es/eta_mec;
PmTGV=PmG+PmV;

eta_cyclexV=P_mcyV/delta_e_fum;
eta_totex = PmTGV/(m_cg*e_c);
eta_cyclenV = P_mcyV / (m_gg* (h_4g - h_5g));
eta_totenV = P_es / Q_HRSG;

pertesEn_mecV = P_mcyV - P_es;
pertesEn_mec_tot = pertesEn_mecG + pertesEn_mecV;

pertesEn_cond = mv_tot * (h_7 - h_1);
pertesEn_chem = Q_exh;

%eta_totenTGV = eta_totenG + eta_totenV - eta_totenG*eta_totenV
eta_totenTGV = (P_es + P_eg) / (m_cg*LHV) ;
eta_gen = delta_h_vap/(m_cg*LHV);
eta_gex = delta_e_vap/(m_cg*e_c);

flux_Ex_HRSG = (e_4g - e_5g) * m_gg;
eta_combex = (m_gg*e_3g - m_ag*e_2g ) /...
    (m_cg*e_c); 
%eta_totexV = P_e_V / flux_Ex_HRSG;  % faux !! ????
%eta_totexTGV = eta_totexG + eta_totexV - eta_totexG*eta_totexV
eta_totexTGV = (P_es + P_eg) / fluxEx_prim;

pertesEx_exh = m_gg*e_5g - e_1g*m_ag;
%pertesEx_exh2 = m_fumV* (1+1/(3.21*17.1))*(data(5).e-data(1).e)
pertesEx_cond = mv_tot*(e_7-e_1);
pertesEx_comb = P_primTGV * (1-eta_combex);
pertesEx_rotV = pertesEx_C + pertesEx_T;
pertesEx_mecV = (1-eta_mec)* (fluxEx_T_V - fluxEx_C_V);
pertesEx_mecG = (1-eta_mec)* (fluxEx_T_G - fluxEx_C_G);
pertesEx_mec = pertesEx_mecV + pertesEx_mecG;
pertesEx_rot_tot = pertesEx_rotV + pertesEx_rotG ; 


ETA(1) = eta_cyclenV;
ETA(2) = eta_cycleng;
ETA(3) = eta_totenTGV;
ETA(4)  = eta_cyclexV;
ETA(5) = eta_cyclexg;
ETA(6)  = eta_totex;
ETA(7)  = eta_gen;
ETA(8)  = eta_gex;
ETA(9)  = eta_combex;
ETA(10) = eta_chemex;
ETA(11) = eta_transex;

MASSFLOW(1) =mv_HP;
MASSFLOW(2) =mv_IP;
MASSFLOW(3) =mv_LP;
MASSFLOW(4) =m_ag;
MASSFLOW(5) =m_cg;


end
