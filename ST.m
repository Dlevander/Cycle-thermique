function [ETA XMASSFLOW DATEN DATEX DAT MASSFLOW COMBUSTION FIG] = ST(P_e,options,display)
% ST Steam power plants modelisation
% ST(P_e,options,display) compute the thermodynamics states for a Steam
% power plant (combustion, exchanger, cycle) turbine based on several 
% inputs (given in OPTION) and based on a given electricity production P_e.
% It returns the main results. It can as well plots graphs if input 
% argument DISPLAY = true (<=> DISPLAY=1)
%
% INPUTS (some inputs can be dependent on others => only one of these 2 can
%         be activated)
% P_E = electrical power output target [kW]
% OPTIONS is a structure containing :
%   -options.nsout     [-] : Number of feed-heating 
%   -options.reheat    [-] : Number of reheating
%   -options.T_max     [°C] : Maximum steam temperature
%   -options.T_cond_out[°C] : Condenseur cold outlet temperature
%   -options.p3_hp     [bar] : Maximum pressure
%   -options.drumFlag  [-] : if =1 then drum if =0 => no drum. 
%   -options.eta_mec   [-] : mecanic efficiency of shafts bearings
%   -options.comb is a structure containing combustion data : 
%       -comb.Tmax     [°C] : maximum combustion temperature
%       -comb.lambda   [-] : air excess
%       -comb.x        [-] : the ratio O_x/C. Example 0.05 in CH_1.2O_0.05
%       -comb.y        [-] : the ratio H_y/C. Example 1.2 in CH_1.2O_0.05
%   -options.T_exhaust [°C] : Temperature of exhaust gas out of the chimney
%   -options.p_3       [-] : High pressure after last reheating
%   -options.x4        [-] : Vapor ratio [gaseous/liquid] (in french : titre)
%   -options.T_0       [°C] : Reference temperature
%   -options.TpinchSub [°C] : Temperature pinch at the subcooler
%   -options.TpinchEx  [°C] : Temperature pinch at a heat exchanger
%   -options.TpinchCond[°C] : Temperature pinch at condenser 
%   -options.Tdrum     [°C] : minimal drum temperature
%   -option.eta_SiC    [-] : Isotrenpic efficiency for compression
%   -option.eta_SiT    [-] : Isotrenpic efficiency for Turbine. It can be a vector of 2 values :
%             	             eta_SiT(1)=eta_SiT_HP,eta_SiT(2)=eta_SiT_others
% DISPLAY = 1 or 0. If 1, then the code should plot graphics. If 0, then 
%          do not plot.
%
%OUPUTS : 
% ETA is a vector with :
%   -eta(1) : eta_cyclen, cycle energy efficiency
%   -eta(2) : eta_toten, overall energy efficiency
%   -eta(3) : eta_cyclex, cycle exegy efficiency
%   -eta(4) : eta_totex, overall exergie efficiency
%   -eta(5) : eta_gen, Steam generator energy efficiency
%   -eta(6) : eta_gex, Steam generator exergy efficiency
%   -eta(7) : eta_combex, Combustion exergy efficiency
%   -eta(8) : eta_chemex, Chimney exergy efficiency (losses)
%   -eta(9) : eta_transex, Heat exchanger overall exergy efficiency
%   FYI : eta(i) \in [0;1] [-]
% Xmassflow is a vector with each feedheating massflow [kg/s] (respect to figure 
%           2.33, page 91 "Thermal Power Plants" English version).
%           Xmassflow(1) = mass flow at 6_1 etc...
% DATEN is a vector with : 
%   -daten(1) : perte_gen [kW]
%   -daten(2) : perte_mec [kW]
%   -daten(3) : perte_cond [kW]
% DATEX is a vector with :
%   -datex(1) : perte_mec    [kW]
%   -datex(2) : perte_totex  [kW]
%   -datex(3) : perte_rotex  [kW]
%   -datex(4) : perte_combex [kW]
%   -datex(5) : perte_condex [kW]
%   -datex(6) : perte_chemex [kW]
%   -datex(7) : perte_transex[kW]
% DAT is a matrix containing :
% dat = {T_1       , T_2       , ...       , T_6_I,     T_6_II, ... ;  [°C]
%        p_1       , p_2       , ...       , p_6_I,     p_6_II, ... ;  [bar]
%        h_1       , h_2       , ...       , h_6_I,     h_6_II, ... ;  [kJ/kg]
%        s_1       , s_2       , ...       , s_6_I,     s_6_II, ... ;  [kJ/kg/K]
%        e_1       , e_2       , ...       , e_6_I,     e_6_II, ... ;  [kJ/kg]
%        x_1       , x_2       , ...       , x_6_I,     x_6_II, ... ;   };[-]
% MASSFLOW is a vector containing : 
%   -massflow(1) = m_a, air massflow [kg/s]
%   -massflow(2) = m_v, water massflow at 2 [kg/s]
%   -massflow(3) = m_c, combustible massflow [kg/s] 
%   -massflow(4) = m_f, exhaust gas massflow [kg/s]
% 
% COMBUSTION is a structure with :
%   -combustion.LHV    : the Lower Heat Value of the fuel [kJ/kg]
%   -combustion.e_c    : the combustible exergie         [kJ/kg]
%   -combustion.lambda : the air excess                   [-]
%   -combustion.Cp_g   : heat capacity of exhaust gas     [kJ/kg/K]
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

%% YOUR WORK

% Exemple of how to use 'nargin' to check your number of inputs
if nargin<3
    display = 1;
    if nargin<2
        options = struct();
        if nargin<1
            P_e = 250e3; % [kW] Puissance énergétique de l'installation
        end
    end
end


% Exemple of how to use (isfield' to check if an option has been given (or
% not)
if isfield(options,'nsout')
     nsout = options.nsout;
else
     nsout = 0;  % [-]
end

if isfield(options,'reheat')
     reheat = options.reheat;
else
    reheat = 0;  % [-]
end

if isfield(options,'T_max')
    T_max = options.T_max;
else
    T_max = 525.0;  % [éC]
end

if isfield(options,'T_cond_out')
     T_cond_out = options.T_cond_out;
else
     T_cond_out = 30.0 ;  % [éC]
end

if isfield(options,'p3_hp')
    p3_hp = options.p3_hp;
else
    p3_hp = 200;  % [bar]
end
%présence ou non d'un tiroir pour le superheater
if isfield(options,'drumFlag')
     drumFlag = options.drumFlag;
else
     drumFlag = 1 ;  % [-]
end

if isfield(options,'eta_mec')
     eta_mec = options.eta_mec;
else
     eta_mec = 0.98;  % [-]
end

if isfield(options,'comb')
    if isfield(options.comb,'Tmax')
        Tmax = options.comb.Tmax;
    else
        Tmax = 500;  % [éC] A MODIFIER
    end
    if isfield(options.comb,'lambda')
        lambda = options.comb.lambda;
    else
        lambda = 1.05;  % [-] 
    end
    if isfield(options.comb,'x')
        x = options.comb.x;
    else
        x = 0;  % [-] CH4
    end
    if isfield(options.comb,'y')
        y = options.comb.y;
    else
        y = 4;  % [-] CH4
    end
else
    Tmax = 500;  % [éC] A MODIFIER
    lambda = 1.05;  % [-]
    x = 0;  % [-] CH4
    y = 4;  % [-] CH4
end

if isfield(options,'T_exhaust')
     T_exhaust = options.T_exhaust;
else
     T_exhaust = 120.0;  % [éC]
end

if isfield(options,'p_3')
     p_3 = options.p_3;
else
     p_3 = 62;  % [bar]
end

if isfield(options,'x4')
     x4 = options.x4;
else
     x4 = 0.88;  % [-]
end

if isfield(options,'T_0')
    T_60 = options.T_0;
else
    T_60 = 15.0;  % [éC]
end

if isfield(options,'TpinchSub')
     TpinchSub= options.TpinchSub;
else
     TpinchSub = 115.0;  % [éC]
end

if isfield(options,'TpinchEx')
     TpinchEx = options.TpinchEx;
else
     TpinchEx = 490.0;  % [éC]
end

if isfield(options,'TpinchCond')
     TpinchCond= options.TpinchCond;
else
     TpinchCond = 18.0;  % [éC]
end

if isfield(options,'Tdrum')
     Tdrum = options.Tdrum;
else
     Tdrum = 30.0;  % [-]
end

if isfield(options,'eta_SiC')
     eta_SiC = options.eta_SiC;
else
     eta_SiC = 1;  % [-]
end

if isfield(options,'eta_SiT')
     eta_SiT_HP = options.eta_SiT(1);
     eta_SiT_others = options.eta_SiT(2);
     
else
     eta_SiT_HP = 0.9;  % [-] A MODIFIER , PLUSIEURS VALEURS ? DIFFERENTS PAR TURBINE ?
     eta_SiT_others = 0.92; % on considère rendement égale pour les turbines MP et BP
end

if P_e == null 
    P_e = 35e3; %[kW]
end

%%%%%%%%%%%%%%% OUTPUT %%%%%%%%%%%%%%%%%%%
ETA = zeros(9,1);
eta_cyclen;
eta_toten;
eta_cyclex;
eta_totex;
eta_gen;
eta_gex;
eta_combex;
eta_chemex;
eta_transex;


Xmassflow = zeros(nsout,1);

DATEN = zeros(3,1);

DATEX = zeros(7,1);

numEtat = 7+reheat+nsout;

DAT= zeros(6,20);

%%%%%%%%%%%%%%% ETAT 30 %%%%%%%%%%%%%%% 

% Sortie chaudière
T_30 = T_max;
p_30 = p3_hp;
h_30 = XSteam('h_pT',p_30,T_30);
s_30 = XSteam('s_pT',p_30,T_30);
%x_30 = XSteam('x_ph',p_30,h_30); = 1 car vapeur surchauffée
e_30 = exergie(h_30,s_30);



%%%%%%%%%%%%%%%% Si resurchauffe %%%%%%%%%%%%%%%
if reheat == 1
    %%%%%%%%%%%%%%% ETAT 40 %%%%%%%%%%%%%%%
    % Sortie de la turbine HP dans cas d'une resurchauffe
    p_40 = p_3; % car surchauffe isobare
    s_40s = s_30;
    x_40s = XSteam('x_ps',p_40,s_40s);
    h_40s = XSteam('h_ps',p_40,s_40s);
    
    h_40 = h_30 - eta_SiT_HP * (h_30-h_40s);
    s_40 = XSteam('s_ph',p_40,h_40);
    x_40 = XSteam('x_ps',p_40,s_40);
    T_40 = XSteam('T_ph',p_40,h_40);
    e_40 = exergie(h_40,s_40);
    
    %%%%%%%%%%%%%%% ETAT 50 %%%%%%%%%%%%%%%
    % Entrée de la turbine MP dans le cas d'une resurchauffe
    
    p_50 = p_3; %la pression apres la resurchauffe
    T_50 = T_30;
    h_50 = XSteam('h_pT',p_50,T_50);
    s_50 = XSteam('s_pT',p_50,T_50);
    %x_50 = XSteam('x_ph',p_50,h_50); = 1 car vapeur surchauffée
    e_50 = exergie(h_50,s_50);
    
    %%%%%%%%%%%%%%% ETAT 60 %%%%%%%%%%%%%%%
    % Sortie de la turbine BP dans cas isentropique (60s) et réel (60)
    s_60s = s_50;
    T_60s = T_cond_out + TpinchCond;
    T_60 = T_cond_out + TpinchCond;
    p_60 = XSteam('psat_T',T_60);
    x_60s = XSsteam('x_ps',p_60,s_60s);
    h_60s = XSteam('h_Tx',T_60s,x_60s);
    
    h_60 = h_50 - eta_SiT_others * (h_50-h_60s) ;
    x_60 = XSteam('x_ph',p_60,h_60);
    s_60 = XSteam('s_ph',p_60,h_60);
    e_60 = exergie(h_60,s_60);

%%%%%%%%%%%%%%%% Si pas resurchauffe %%%%%%%%%%%%%%%
elseif reheat == 0
    %%%%%%%%%%%%%%% ETAT 60S + 60  %%%%%%%%%%%%%%%
    % Sortie de la turbine BP (60) dans cas isentropique (60s) et réel (60)
    
    s_60s = s_30;
    T_60s = T_cond_out + TpinchCond;
    T_60 = T_cond_out + TpinchCond;
    p_60 = XSteam('psat_T',T_60);
    x_60s = XSsteam('x_ps',p_60,s_60s);
    %x_40s = (s_40s - XSteam('sL_T',T_40s)) / (XSteam('sV_T',T_40s) - XSteam('sL_T',T_40s)) ;
    h_60s = XSteam('h_Tx',T_60s,x_60s);
    %h_40s = (x_40s*XSteam('hV_T',T_40s)) + ((1-x_40s)*XSteam('hL_T',T_40s));
    h_60 = h_30 - eta_SiT_HP * eta_SiT_others * (h_30-h_60s) ;
    x_60 = XSteam('x_ph',p_60,h_60);
    %x_40 = (h_40 - XSteam('hL_T',T_40s)) / (XSteam('hV_T',T_40s) - XSteam('hL_T',T_40s));
    s_60 = XSteam('s_ph',p_60,h_60);
    %s_40 = (x_40*XSteam('sV_T',T_40)) + ((1-x_40)*XSteam('sL_T',T_40));
    e_60 = exergie(h_60,s_60);
end




end