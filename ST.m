function [ETA,XMASSFLOW,DATEN,DATEX,DAT,MASSFLOW,COMBUSTION,FIG] = ST(P_e,options,display)
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
%   -options.p_3       [] : High pressure after last reheating
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
% dat = {T_1       , T_2       , ...       , T_6_I,     T_6_II, ... ;  [Â°C]
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
    display == 1;
    if nargin<2
        options = struct();
        if nargin<1
            
            P_e = 250e3; % [kW] Puissance energetique de l'installation
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
    
    T_max = 525.0;  % [C]
    
end

if isfield(options,'T_cond_out')
    T_cond_out = options.T_cond_out;
else
    T_cond_out = 30.0 ;  % [C]
end

if isfield(options,'p3_hp')
    p3_hp = options.p3_hp;
else
    p3_hp = 200;  % [bar]
end
%presence ou non d'un tiroir pour le superheater
if isfield(options,'drumFlag')
    drumFlag = options.drumFlag;
else
    drumFlag = 1 ;  % [-] % cas de base sans degazificateur non ?
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
        Tmax = 500;  % [C] A MODIFIER
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
    Tmax = 500;  % [C] A MODIFIER
    lambda = 1.05;  % [-]
    x = 0;  % [-] CH4
    y = 4;  % [-] CH4
end

if isfield(options,'T_exhaust')
    T_exhaust = options.T_exhaust;
else
    T_exhaust = 120.0;  % [C]
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
    T_0 = options.T_0;
else
    T_0 = 15.0;  % [C]
end

if isfield(options,'TpinchSub')
    TpinchSub= options.TpinchSub;
else
    TpinchSub = 115.0;  % [C]
end

if isfield(options,'TpinchEx')
    TpinchEx = options.TpinchEx;
else
    TpinchEx = 490.0;  % [C]
end

if isfield(options,'TpinchCond')
    TpinchCond= options.TpinchCond;
else
    TpinchCond = 18.0;  % [C]
end

if isfield(options,'Tdrum')
    Tdrum = options.Tdrum;
else
    Tdrum = 120.0;  % [�C]
end

if isfield(options,'eta_SiC')
    eta_SiC = options.eta_SiC;
else
    eta_SiC = 1;  % [-]
end

if isfield(options,'eta_SiT')
    if length(options.eta_SiT) == 2
        eta_SiT_HP = options.eta_SiT(1);
        eta_SiT_others = options.eta_SiT(2);        
    elseif length(options.eta_SiT) == 1
        eta_SiT_HP = options.eta_SiT; % A MODIFIER SUREMENT
        eta_SiT_others = eta_SiT_HP;
    end
else
    eta_SiT_HP = sqrt(0.88);  % [-]
    eta_SiT_others = sqrt(0.88); % on considere rendement egale pour les turbines MP et BP
end

if P_e <= 0
    P_e = 35e3; %[kW]
end

%%%%%%%%%%%%%%% OUTPUT %%%%%%%%%%%%%%%%%%%
ETA = zeros(9,1);
XMASSFLOW = zeros(nsout,1);
DATEN = zeros(3,1);
DATEX = zeros(7,1);
%numEtat = 7+reheat+nsout;

DAT= zeros(6,35);

MASSFLOW = 0; %A MODIFIER
COMBUSTION = 0; %A MODIFIER

FIG = 0; %A MODIFIER
%% Calcul

%%%%%%%%%%%%%%% ETAT 30 %%%%%%%%%%%%%%%

% Sortie chaudiere
T_30 = T_max;
p_30 = p3_hp;
h_30 = XSteam('h_pT',p_30,T_30);
s_30 = XSteam('s_pT',p_30,T_30);
x_30 = NaN;% XSteam('x_ph',p_30,h_30); = NaN car vapeur surchauffee
e_30 = exergie(h_30,s_30);

%%%%%%%%%%%%%%% ETAT 22 %%%%%%%%%%%%%%%
p_22= p_30;% surchauffe isobare
T_22= XSteam('Tsat_p',p_22);% Tsat
h_22= XSteam('hV_p',p_22);
s_22= XSteam('sV_p',p_22);
e_22= exergie(h_22,s_22);
x_22= 1;

%%%%%%%%%%%%%%% ETAT 21 %%%%%%%%%%%%%%%
T_21 = T_22; % evaporation isoterme
p_21 = p_22; %XSteam('psat_T',T_21)
h_21 = XSteam('hL_T',T_21);
s_21= XSteam('sL_p',p_21);
e_21= exergie(h_22,s_22);
x_21= 0;
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
    % Entree de la turbine MP dans le cas d'une resurchauffe
    
    p_50 = p_3; %la pression apres la resurchauffe
    T_50 = T_30;
    h_50 = XSteam('h_pT',p_50,T_50);
    s_50 = XSteam('s_pT',p_50,T_50);
    x_50 = XSteam('x_ph',p_50,h_50);% = NaN car vapeur surchauffee
    e_50 = exergie(h_50,s_50);
    
    %%%%%%%%%%%%%%% ETAT 60 %%%%%%%%%%%%%%%
    % Sortie de la turbine BP dans cas isentropique (60s) et reel (60)
    s_60s = s_50;
    T_60 = T_cond_out;%-3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% pour RH1
    p_60 = XSteam('psat_T',T_60);
    x_60s = XSteam('x_ps',p_60,s_60s);
    h_60s = XSteam('h_Tx',T_60,x_60s);
    
    h_60 = h_50 - eta_SiT_others * (h_50-h_60s) ;
    x_60 = XSteam('x_ph',p_60,h_60);
    s_60 = XSteam('s_ph',p_60,h_60);
    e_60 = exergie(h_60,s_60);
    
    %%%%%%%%%%%%%%%% Si pas resurchauffe %%%%%%%%%%%%%%%
elseif reheat == 0
    
    %%%%%%%%%%%%%%% ETAT 60S + 60  %%%%%%%%%%%%%%%
    % Sortie de la turbine BP (60) dans cas isentropique (60s) et reel (60)
    s_60s = s_30;
    
    % En se basant sur T_cond
    T_60 = T_cond_out;
    p_60 = XSteam('psat_T',T_60);
    x_60s = XSteam('x_ps',p_60,s_60s);
    h_60s = XSteam('h_Tx',T_60,x_60s);
    %     h_60b = h_30 - eta_SiT_HP * eta_SiT_others * (h_30-h_60s);
    %     h_N = (h_60s*eta_SiT_HP * eta_SiT_others - h_60b)/(eta_SiT_HP * eta_SiT_others - 1); %determination du point N en vue de la regle de Bauman
    %     eta_vap_humide = x4*eta_SiT_HP * eta_SiT_others;
    %     h_60 = h_N - eta_vap_humide * (h_N-h_60s);
    h_60 = h_30 - eta_SiT_HP * eta_SiT_others * (h_30-h_60s);
    x_60 = XSteam('x_ph',p_60,h_60);
    s_60 = XSteam('s_ph',p_60,h_60);
    e_60 = exergie(h_60,s_60);
    
end

%%%%%%%%%%%%%%% ETAT 70 %%%%%%%%%%%%%%%
%Sortie du condenseur, liquide sature
T_70 = T_60; %car passe simplement de vapeur a liquide saturee -3 pour RH1, -1 pour 1reheat8nsout
h_70 = XSteam('hL_T',T_70);
p_70 = XSteam('psat_T',T_70);
s_70 = XSteam('sL_T',T_70);
e_70 = exergie(h_70,s_70);
x_70=XSteam('x_ph',p_70,h_70); %=0 car liquide saturee en sortie de condenseur

if nsout==0
    %%%%%%%%% Calcul etat  10 %%%%%%%%%%
    %S'il n'y a pas de soutirage le point 1 = point 7 car on coupe avant Pe et on revient au point 1 avant Pa
    T_10 = T_70;
    h_10 = h_70;
    p_10 = p_70;
    s_10 = s_70;
    e_10 = e_70;
    x_10 = x_70;
    
    %%%%%%%%% Calcul etat  20 %%%%%%%%%%
    p_20 = p_21+4; %hypothese %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% pour RH1
    h_20s = XSteam('h_ps',p_20,s_10); %h dans le cas isentropique
    h_20=h_10+(h_20s-h_10)/eta_SiC;%h en prenant compte rendement is
    s_20=XSteam('s_ph',p_20,h_20);
    T_20=XSteam('T_ph',p_20,h_20);
    x_20=XSteam('x_ph',p_20,h_20);
    e_20=exergie(h_20,s_20);
    
elseif nsout>0
    %%%%%%%%% Calcul etat 61s et 61 %%%%%%%%
    % nsout -1 pour retirer le sout en sortie de HP, +2 pour le linspace
    
    h6s_temp = linspace(h_60s,h_50,(nsout+1)); %h_6x, x plus bas => enthalpie plus basse
    h6s = h6s_temp(2:length(h6s_temp)-1);
    s_60 = s_50;
    p6s = zeros(1,length(h6s));
    
    for i=1:length(h6s)
        p6s(i) = XSteam('p_hs',h6s(i),s_60); % XSteam, pas de vecteur en arg
    end
    p6 = p6s; %hypothese
    h6=h_50+eta_SiC*(h6s-h_50);
    
    % l'echangeur en sortie de HP a son etat deja defini
    h6 = [h6 h_40];
    p6 = [p6 p_40];
    %x_60 = x_40;
    %preallocation
    T6 = zeros(1,length(h6)-1);
    x6 = zeros(1,length(T6));
    e6 = zeros(1,length(T6));
    for i=1:length(h6)-1 % l'etat 40 est deja defini
        T6(i) = XSteam('T_hs',h6(i),s_60);
        x6(i) = XSteam('x_ph',p6(i),h6(i));
        e6(i) = exergie(h6(i),s_60);
    end
    %% Degazification
    %%%%%%%%% Calcul etat 7x %%%%%%%%%%
    p7 = p6; % on considere les vannes de detente apres les etats 7
    %preallocation
    T7 = zeros(1,length(T6));
    h7 = zeros(1,length(T6));
    % les temperature de condensation sont les temp de sortie des echangeurs
    flag = 0;
    p_degaz = 0;
    d = 0;
    for i=1:length(p7)
        T7(i) = XSteam('Tsat_p',p7(i));
        if T7(i) >= 393.15 && flag==0 %flag pour s'assurer que ce ne soit pas letat avec la temp max qui soit retenu mais bien celui juste au dessus de 120Â°C
            % calcul de la pression dans le degazificateur
            T_degaz = T7(i);
            p_degaz = p7(i); %temperature avant la pompe pb
            flag = 1;
            d = i ; % endroit du degazificateur
        end
        h7(i) = XSteam('h_pT',p7(i),T7(i));
    end
    
    %%%%%%%%%%%%%%% ETAT 80 %%%%%%%%%%%%%%%
    %Apres la pompe Pe, de rendement isentropique option.eta_SiC
    p_80 = p_degaz ;
    h_80s = XSteam('h_ps',p_80,s_70); %h dans le cas isentropique
    h_80=h_70+eta_SiC*(h_80s-h_70);%h en prenant compte rendement is
    s_80=XSteam('s_ph',p_80,h_80);
    T_80=XSteam('T_ph',p_80,h_80);
    x_80=XSteam('x_ph',p_80,h_80);
    e_80=exergie(h_80,s_80);
    
    %%%%%%%%% Calcul etat  10 %%%%%%%%%%
    T_10 = T7(length(T7)-1)-TpinchEx;
    p_10 = p_degaz + (p_21 - p_degaz)/3; % hypothese
    h_10 = XSteam('hL_T',T_10);
    s_10 = XSteam('sL_T',T_10);
    x_10 = NaN;
    e_10 = exergie(h_10,h_10);
    
    %%%%%%%%% Calcul etat  20 %%%%%%%%%%
    p_20 = p_21; %hypothese
    h_20s = XSteam('h_ps',p_20,s_10); %h dans le cas isentropique
    h_20=h_10+eta_SiC*(h_20s-h_10);%h en prenant compte rendement is
    s_20=XSteam('s_ph',p_20,h_20);
    T_20=XSteam('T_ph',p_20,h_20);
    x_20=XSteam('x_ph',p_20,h_20);
    e_20=exergie(h_20,s_20);
    
    %%%%%%%%% Calcul des etat 9x %%%%%%%%
    %preallocation
    p9 = zeros(1,nsout);
    s9 = zeros(1,nsout);
    e9 = zeros(1,nsout);
    
    T9 = T7-TpinchEx; %si on considere des echangeurs parfaits
    for i = 1:nsout
        if i<d
            p9(i) = p_degaz;
        else
            p9(i) = p_10;
        end
        h9(i) = XSteam('h_pT',p9(i),T9(i));
        s9(i) = XSteam('s_pT',p9(i),T9(i));
        e9(i) = exergie(h9(i),s9(i));
    end
    
    %%%%%%%%% Calcul etat 100 %%%%%%%%%%
    %etat avant la vanne et aprÃ¨s le subcooler
    T_100 = T_80 + TpinchSub;
    p_100 = p7(1) ;
    h_100 = XSteam('h_pT',p_100,T_100);
    s_100 = XSteam('s_pT',p_100,T_100);
    x_100 = XSteam('x_ph',p_100,T_100);
    e_100 = exergie(h_100,s_100);
    
    %%%%%%%%% Calcul etat  110 %%%%%%%%%%
    %etat apres la vanne (isenthalpique), juste avant le condenseur
    h_110 = h_100;
    T_110 = T_70;
    p_110 = p_70;
    s_110 = XSteam('s_pT',p_110,T_110);
    x_110 = XSteam('x_ph',p_110,T_110);
    e_100 = exergie(h_110,s_110);
    
    %%%%%%%%% Calcul des Soutirages X %%%%%%%%%%
    A = zeros(nsout,nsout);
    B = zeros(nsout,1);
    %remplissage de la matrice A et B
    %premiere ligne particuliere
    A(1,1:d-1) = A(1,1:d-1) + ( h9(1) - h_80 - h_100 + h7(1) );
    A(1,1) = A(1,1) - ( h7(1) -h6(1) );
    A(1,2:d-1) = A(1,2:d-1) - ( h7(1) - h7(2) );
    
    B(1) = h9(1)-h_80-h_100+h7(1);
    for i=2:nsout % i les lignes de la matrice
        deltaH9  = h9(i) - h9(i-1);
        deltaH7  = h7(i) - h7(i+1);
        deltaH76 = h7(i) - h6(i);
        if i < d
            A(i,1:d-1) = A(i,1:d-1) + deltaH9;
            A(i,i) = A(i,i) - deltaH76; % remplissage diagonale
            A(i,i+1:d-1) = A(i,i+1:d-1) - deltaH7;
            
            B(i) = deltaH9;
        end
        if i > d
            A(i,:) = A(i,:) + deltaH9;
            A(i,i) = A(i,i) - deltaH76; % remplissage diagonale
            if i < nsout
                A(i,i+1:nsout) = A(i,i+1:nsout) - deltaH7;
            end
        end
    end
    % remplissage ligne d particuliere
    A(d,1:d-1) = A(d,1:d-1) + h7(d) - h9(d-1);
    A(d,d) = h7(d) - h6(d);
    A(d,d+1:nsout) = A(d,d+1:nsout) + h7(d) - h7(d+1);
    
    B(d) = h7(d) - h9(d-1);
    % Resolution
    XMASSFLOW = A\B;
    
    %%%%%%%%% Calcul etat 90 %%%%%%%%%%
    h_90 = h_80 + (h_100 - h7(1)) * sum( XMASSFLOW(1:d-1) )/( 1+ sum( XMASSFLOW(1:d-1) ));
end
%% remplissage output

DAT(:,1) = [T_10 p_10 h_10 s_10 e_10 x_10]';
DAT(:,2) = [T_20 p_20 h_20 s_20 e_20 x_20]';
DAT(:,3) = [T_21 p_21 h_21 s_21 e_21 x_21]';
DAT(:,4) = [T_22 p_22 h_22 s_22 e_22 x_22]';
DAT(:,5) = [T_30 p_30 h_30 s_30 e_30 x_30]';
if reheat == 0
    DAT(:,6) = [T_60 p_60 h_60 s_60 e_60 x_60]';
    if nsout > 0
        
    end    
elseif reheat == 1
    DAT(:,6) = [T_40 p_40 h_40 s_40 e_40 x_40]';
    DAT(:,7) = [T_50 p_50 h_50 s_50 e_50 x_50]';
    DAT(:,8) = [T_60 p_60 h_60 s_60 e_60 x_60]';
    if nsout > 0
        
    end  
end
%% Rendements
X_tot = sum(XMASSFLOW);
    %%%%%%%%% Calcul des rendements %%%%%%%%%%
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
    
    %  ex_mT    exergie du travail moteur de la turbine [kJ/kg]
    %  W_mT     Travail moteur de la turbine
    %  Q_I      Action calorifique a la chaudiere [kJ/kg]
    %  ex_I     delta d' Exergie�a la chaudiere [kJ/kg]
    ex_mT = 0;
    if reheat == 0
        W_mT = h_30-h_60;
        Q_I = (X_tot+1)*(h_30 - h_20);
        ex_I = (X_tot+1)*(e_30 - e_20);
        if nsout > 0
            for i = 1:nsout
                W_mT = W_mT + XMASSFLOW(i)*(h_30 - h6(i));
                ex_mT = ex_mT + XMASSFLOW(i)*(e_30 - e6(i));
            end
            ex_mT = ex_mT + e_30 - e_60;
        elseif reheat == 1
            
            W_mT = (X_tot+1)*(h_30 - h_40) + (h_50 - h_60);
            Q_I = (X_tot+1)*(h_30 - h_20) + (sum(XMASSFLOW(1:nsout-1))+1)*(h_50 - h_40);% par kg
            ex_I =(X_tot+1)*(e_30 - e_20) + (sum(XMASSFLOW(1:nsout-1))+1)*(e_50-e_40);
            
            if nsout > 1
                for i = 1:(nsout-1)
                    W_mT = W_mT + XMASSFLOW(i)*(h_50 - h6(i));
                    ex_mT = ex_mT + (XMASSFLOW(i)*(e_50 - e6(i)));
                end
            end
            
            ex_mT = ex_mT  +(X_tot+1)*(e_30-e_40) + e_50-e_60;
        end
        ex_mT = ex_mT + e_30 - e_60;
        
        % W_mP     travail fourni par les pompes
        %  ex_mP    exergie des pompes
        if  nsout ~= 0
            if drumFlag ==0 % travail fourni aux 2 pompes sans bache et avec soutirages
                W_mP = (X_tot+1)*(h_80-h_70+h_20-h_10);
                ex_mP = (X_tot+1)*(e_80-e_70+e_20-e_10);%W_mP - T_riv*(data(80).s - data(70).s + data(20).s - data(10).s);
            else % travail fourni aux 3 pompes avec bache et avec soutirages
                W_mP = (X_tot+1)*(h_80-h_70 + h9(d) - h7(d) + h_20-h_10);
                ex_mP = W_mP - T_0*(s_80-s_70 + s9(d) - s7(d) + s_20-s_10);
            end
            
        else % travail fourni a� la pompe sans bache et sans soutirages
            W_mP = (X_tot+1)*(h_20-h_10);
            ex_mP = W_mP - T_0*(s_20-s_10);
        end
        
        % m_vap    debit massique de vapeur
        m_vap = P_e/(eta_mec*(W_mT-W_mP)); %page 60
        
        %%%%%%%% Combustion  %%%%%%
        %%determination des fractions massiques des fumees
        
        if x == 0 && y == 4 % CH4 + 2*lambda*(O2+3.76N2) => 2*(lambda-1)*O2 + CO2 + 2*H2O + 2*lambda*3.76*N2
            LHV = 5.020625*10^4 ; % kJ/kg
            PCI_comb = 8.033*10^5; %J/mole
            e_c = 52215; %kJ/kg fuel exergy. tableau pg 25 du livre
            MmComb = 16 ; %(kg/kmol)
            Cp_comb = 1000*35.8/MmComb; %J/(kg_comb*K)
            % Fractions massiques des reactifs (kg/kg_comb)
            Comb_R = 1;
            O2_R = 2*lambda*32;
            N2_R = 2*lambda*3.76*28;
            % Fractions massiques des produits (g/mol_comb)
            x_O2 = 2*(lambda-1)*32;
            x_CO2 = 44;
            x_H2O = 2*18;
            x_N2 = 2*lambda*3.76*28;
            Sum = x_O2 + x_CO2 + x_H2O + x_N2 ; %TEST
            
            
        elseif x == 0 && y == 0 % C + lambda*(O2+3.76N2) => (lambda-1)*O2 + CO2 + lambda*3.76*N2  p131
            LHV = 32780; % [kJ/kg]
            MmComb = 12; %[kg/kmol]
            Cp_comb = 1000*10.4/MmComb; %[J/kg_comb*K]
            e_c = 32400; %[kJ/kg]
            % Fractions massiques des reactifs [kg/kg_comb]
            Comb_R = 1;
            O2_R = lambda*32/MmComb;
            N2_R = lambda*3.76*28/MmComb;
            % Fractions massiques des produits [kg/kg_comb]
            x_O2 = (lambda-1)*32/MmComb;
            x_CO2 = 44/MmComb;
            x_H2O = 0;
            x_N2 = lambda*3.76*28/MmComb;
            
        elseif x == 0 && y == 8/3 % Propane : C3H8 + 5*lambda*(O2+3.76N2) => 5*(lambda-1)*O2 + 3*CO2 + 4*H2O + 5*lambda*3.76*N2
            
            LHV = 46465; % [kJ/kg]
            MmComb = 44; %[kg/kmol]
            CpComb = 1000*70.9/MmComb; %[J/kg_comb*K]
            e_c = 49045; %[kJ/kg]
            % Fractions massiques des reactifs (kg/kg_comb)
            Comb_R = 1;
            O2_R = 5*lambda*32/MmComb;
            N2_R = 5*lambda*3.76*28/MmComb;
            
            % Fractions massiques des produits (kg/kg_comb)
            x_O2 = 5*(lambda-1)*32/MmComb;
            x_CO2 = 3*44/MmComb;
            x_H2O = 4*18/MmComb;
            x_N2 = 5*lambda*3.76*28/MmComb;
            
        end
        %Soit on donne lambda soit Tf, liee par p.207 eq 10.6 dans l'ancien
        %bouquin epsilon = 0 et hc =0
        ma1 = (32+3.76*28)*(1+(y/4))/(12+y); % pouvoir comburivore [kg_air_stoech/kg_comb] % x toujours egal a 0
        %ma1 = (1+(y-2*x)/4)*(32+3.76*28.15)/(12.01+1.008*y+16*x);
        combustion.LHV = LHV;
        combustion.e_c = e_c;
        combustion.lambda = lambda;
        %%Debit massique combustible et fumees
        T_exh = T_exhaust;%T en sortie de cheminee
        T_f = Tmax; %T fumees juste en sortie de combustion
        TK_exh = T_exh+273.15;
        TK_f = T_f+273.15;
        CpMoyO2_f = 1000*mean(janaf('c','O2',linspace(TK_exh,TK_f,50)));
        CpMoyCO2_f = 1000*mean(janaf('c','CO2',linspace(TK_exh,TK_f,50)));
        CpMoyN2_f = 1000*mean(janaf('c','N2',linspace(TK_exh,TK_f,50)));
        CpMoyH2O_f = 1000*mean(janaf('c','H2O',linspace(TK_exh,TK_f,50)));
        CpMoy_f = (CpMoyH2O_f*x_H2O + CpMoyCO2_f*x_CO2 + CpMoyO2_f*x_O2+CpMoyN2_f*x_N2)/Sum; % somme pondere non ?
        delta_h=CpMoy_f*(T_f-T_exh);%kj/kg %inverse par rapport a avant
        m_fum = m_vap*Q_I*1e3/delta_h; %debit fumees
        %m_comb = m_fum/(1+ma1*lambda); %debit combustible
        m_comb = m_vap*Q_I/(LHV*0.945); %[kg/s]
        m_a = lambda*ma1*m_comb; %[kg/s]
        %%Calcul des enthalpies, entropies et exergies des fumees et du fuel+air
        %on change les fractions massiques de (kg/kg_comb) (kg/kg_fum)

        x_O2 = x_O2/Sum;
        x_CO2 = x_CO2/Sum;
        x_N2 = x_N2/Sum;
        x_H2O = x_H2O/Sum;
        % Enthalpie, entropie et exergie des fumees juste apres la combustion
        h_f = x_O2*janaf('h','O2',T_f+273.15) + x_CO2*janaf('h','CO2',T_f+273.15) + x_H2O*janaf('h','H2O',T_f+273.15) + x_N2*janaf('h','N2',T_f+273.15);
        s_f = x_O2*janaf('s','O2',T_f+273.15) + x_CO2*janaf('s','CO2',T_f+273.15) + x_H2O*janaf('s','H2O',T_f+273.15) + x_N2*janaf('s','N2',T_f+273.15);
        e_f = exergie(h_f,s_f);
        % Enthalpie, entropie et exergie des fumees en sortie de cheminee
        h_exh = x_O2*janaf('h','O2',T_exh+273.15) + x_CO2*janaf('h','CO2',T_exh+273.15) + x_H2O*janaf('h','H2O',T_exh+273.15) + x_N2*janaf('h','N2',T_exh+273.15);
        s_exh = x_O2*janaf('s','O2',T_exh+273.15) + x_CO2*janaf('s','CO2',T_exh+273.15) + x_H2O*janaf('s','H2O',T_exh+273.15) + x_N2*janaf('s','N2',T_exh+273.15);
        e_exh = exergie(h_exh,s_exh);
        % Enthalpie, entropie et exergie du mÃ©lange fuel+air p32
        e_r = 0; %Pris Ã  letat de reference
        
        %%Rendement du boiler
        if reheat ==0
            rend_boiler = m_vap*(X_tot+1)*(h_30-h_20)/(m_comb*LHV*10^3);
        elseif reheat ==1
            if n_sout>0
                rend_boiler = (m_vap*((X_tot+1)*(h_30-h_20)+ (X_tot-X(nsout)+1)*(h_50-h_40)))/(m_comb*LHV*10^3);
            else
                rend_boiler = (m_vap*((X_tot+1)*(h_30-h_20)+ (X_tot+1)*(h_50-h_40)))/(m_comb*LHV*10^3);
            end
        end
        %%Rendement Ã©nergÃ©tique du cycle
        rend_cyclen = (W_mT - W_mP)/Q_I;
        rend_toten = rend_boiler*eta_SiT_HP*rend_cyclen;
        %% pertes et puissances
        % Pu_tot 	Puissance Totale
        % Per_boiler 	pertes Ã  la chaudiÃ¨re
        % Pu_turb    Puissance Ã  la turbine
        % Pe_meca 	pertes mecanique
        % Pe_cond    pertes au condenseur
        Pu_tot = m_comb*LHV*10^3;
        Per_boiler = Pu_tot*(1-rend_boiler);
        Pu_turb = P_e/eta_SiT_HP;
        Per_meca = Pu_turb*(1-eta_SiT_HP) + W_mP;
        Per_cond =  Pu_tot - P_e - Per_boiler - Per_meca;
        
        %% Rendement exergÃ©tique du cycle
        % P_totex    Flux exergetique total
        P_totex = m_comb*e_c*10^3;
        % rend_totex    rendement exergetique total
        % rend_I_ex     rendement exergetique du generateur
        % rend_comb_ex  rendement exergetique de la combustion
        % rend_chim_ex  rendement exergetique de la cheminÃ©e; P.65
        % rend_trans_ex rendement exergetique du transfert de chaleur; P.65
        % rend_rot_ex   rendement exergetique de la turbine; P.65
        % rend_cycl_ex  rendement exergetique du cycle; P.65
        rend_totex = P_e/P_totex;
        rend_gex = m_vap*ex_I/(m_comb*e_c*10^3);
        rend_combex =( m_fum*(e_f-e_r))/(P_totex);%unitÃ©?
        rend_chemex = (e_f-e_exh)/(e_f-e_r);
        rend_transex = (m_vap*ex_I)/(m_fum*(e_f-e_exh));
        rend_rotex = (W_mT-W_mP)/(ex_mT-ex_mP);
        rend_cyclex = rend_rotex*(ex_mT-ex_mP)/ex_I;
        
        %% Pertes exergÃ©tiques
        % perteex_pompe   pertes exergetiques du aux  pompes
        % perteex_comb    pertes exergetiques du Ã  la combustion
        % perteex_chimn   pertes exergetiques du Ã  la cheminÃ©e
        % pertex_trans    pertes exergetiques du au transfert de chaleur
        perte_turbex = P_e/eta_SiT_HP*(1/rend_rotex-1);
        perte_combex = P_totex*(1-rend_combex);
        perte_chemex = (e_exh-e_f)-(h_exh-h_f);
        
        %if n_bache_sout>0
        %perteEx_cond = abs((( sum(X(1:n_bache_sout))*data(140).e+data(60).e) -data(70).e)-(( sum(X(1:n_bache_sout))*data(140).h+data(60).h) -data(70).h)) ;
        %  else
        %  perteEx_cond =( (data(70).e -data(60).e)-(data(70).h -data(60).h));
        %end
        % perteEx_pompe = ((data(20).h - data(10).h) - data(20).e-data(10).e);
        %if n_sout> 0
        %   perteEx_pompe = ( perteEx_pompe+(data(80).h - data(70).h) - (data(80).e-data(70).e ));
        %  if n_bache_sout>0
        %     perte_pompex = perte_pompex + ((data(90+n_bache_sout).h-data(70+n_bache_sout).h)-(data(90+n_bache_sout).e-data(70+n_bache_sout).e));
        %end
        % end
        %       perte_pompex = abs(perte_pompex)*1000;
        %      perte_condex = perte_condex *1000;
        %
        %           perte_transex = P_totex - P_eff - perte_combex - perte_chemex - perte_turbex - perte_condex - Per_meca;% m_f*(ex_f - ex_exh) - m_vap*ex_I;
        
        %%%%%% Remplissage des vecteurs %%%%%
        ETA(1) = rend_cyclen;
        ETA(2)= rend_toten;
        ETA(3) = rend_cyclex;
        ETA(4)  = rend_totex;
        ETA(5) = rend_boiler;
        ETA(6)= rend_gex;
        ETA(7) =rend_combex;
        ETA(8) =rend_chemex;
        ETA(9) =rend_transex;
        
        DATEN(1) = Per_boiler;
        DATEN(2) = Per_meca;
        DATEN(3) = Per_cond;
        
        DATEX(1) = 0;
        DATEX(2) = 0;
        DATEX(3) = 0;
        DATEX(4) = perte_combex;
        DATEX(5) = 0;
        DATEX(6) = perte_chemex;
        DATEX(7) = 0;
        
        MASSFLOW(1) = m_a; %ma1; % faux m_air pas ma1
        MASSFLOW(2) = m_vap;
        MASSFLOW(3) = m_comb;
        MASSFLOW(4) = m_fum;
        
        combustion.LHV = LHV;
        combustion.e_c = e_c;
        combustion.lambda = lambda;
        combustion.Cp_g = CpMoy_f;
        combustion.fum(1) = x_O2*m_fum;
        combustion.fum(1) = x_N2*m_fum;
        combustion.fum(1) = x_CO2*m_fum;
        combustion.fum(1) = x_H2O*m_fum;
        
        
        
    end
    %% Display part
    if display == 1
        T = linspace(0.001,400,400);
        SL = zeros(1,400);
        SV = zeros(1,400);
        
        for i=1:length(T)
            SL(i) = XSteam('sL_T',T(i));
            SV(i) = XSteam('sV_T',T(i));
        end
        FIG(1) = figure;
        hold on
        %cloche de base
        plot(SL,T,'-b',SV,T,'-b')
        linS = [linspace(DAT(4,1),DAT(4,5),1000)];% linspace(DAT(4,2),DAT(4,3),1000) linspace(DAT(4,3),DAT(4,4),1000) linspace(DAT(4,4),DAT(4,5),1000) linspace(DAT(4,5),DAT(4,6),100)];
        linT = arrayfun( @(s) XSteam('T_ps',DAT(2,2),s),linS);
        plot(DAT(4,:),DAT(1,:),'x')
        plot(linS,linT)
        grid on
    end

end


