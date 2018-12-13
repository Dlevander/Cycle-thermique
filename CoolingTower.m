function [DAT_WATER DAT_AIR MASSFLOW] = CoolingTower(P_w,options)
% COOLINGTOWER is a cooling tower 0D modelisation
% COOLINGTOWER(P_w,options) compute the thermodynamics states for a Cooling
% tower based on several inputs (given in OPTION) and based on a given
% water power to dissipate P_w.
% It returns the main results.
%
% INPUTS :
% P_W = Heat power output at the condenser [kW]
% OPTIONS is a structure containing :
%   -options.Tcond  [C]: Temperature in the condenser
%   -options.Tpinch [C]: Minimum tempearture pinch between Tw_out and the
%                         condenser temperature.
%   -options.Tw_out [C]: Cooling water temperature at the condenser outlet
%   -options.Tw_in  [C]: Cooling water temperature at the condenser inlet
%   -options.Triver [C]: River temperature
%   -options.Ta_in  [C]: Atmospheric air temperature
%   -options.Ta_out [C]: Air outlet temperature of cooling tower
%   -options.Phi_atm [-]: Relative humidity of atmospheric air
%   -options.Phi_out [-]: Maximum relative humidity of air at the cooling
%                         tower outlet.
%
% OUTPUT :
% MassFlow [kg/s]: Vector containing the different massflow :
%   -massflow(1) : water massflow at the condenser
%   -massflow(2) : additionnal water massflow = water flow evaporated
%   -massflow(3) : air massflow at the cooling tower
%
%  dat_water = [T_e1       , T_e2       , T_e3       , T_e4;  %[°C]
%               h_e1       , h_e2       , h_e3       , h_e4;  %[kJ/kg]
%               m_e1       , m_e2       , m_e3       , m_e4]; %[kg/s]
%
%  dat_air   = [Ta_in       , Ta_out  ;  %[°C]
%               ha_in       , ha_out  ;  %[kJ/kg]
%               xa_in       , xa_out  ;  %[kg_water/kg_dry_air]
%               Phia_in     , Phia_out]; %[-] relative humidity
%
%
% ADDITIONNAL INFORMATIONS
% Water points :
%       1 : water outlet of cooling tower
%       2 : water just before condenser
%       3 : water just after  condenser
%       4 : water from the river (coming between 1 & 2)
%
% Air points :
%       a_in : air at the cooling tower inlet
%       a_out : air at the cooling tower outlet
%

%% YOUR WORK
if nargin<2
    options = struct();
    if nargin<1
        P_w = 444e3; %[kW]
    end
end

if isfield(options,'Tcond')
    Tcond = options.Tcond;
else
    Tcond = 30; %[C]
end

if isfield(options,'Tpinch')
    Tpinch = options.Tpinch;
else
    Tpinch = 4; %[C]
end

if isfield(options,'Tw_out')
    Tw_out = options.Tw_out;
else
    Tw_out = 42.8; %[C]
end

if isfield(options,'Tw_in')
    Tw_in = options.Tw_in;
else
    Tw_in = 4; %[C]
end

if isfield(options,'Triver')
    Triver = options.Triver;
else
    Triver = 15; %[C]
end

if isfield(options,'Ta_in')
    Ta_in = options.Ta_in;
else
    Ta_in = 15; %[C]
end

if isfield(options,'Ta_out')
    Ta_out = options.Ta_out;
else
    Ta_out = 25; %[C]
end

if isfield(options,'Phi_atm')
    Phi_atm = options.Phi_atm;
else
    Phi_atm = 0.8; %[K]
end

if isfield(options,'Phi_out')
    Phi_out = options.Phi_out;
else
    % minimisation du debit d'air necessaire avec
    Phi_out = 1; %[K]
end
Cp_el = 4.186; % [kj/kg/K] Cp eau liquide suppose constant
%% etat a_in : air at the cooling tower inlet

psat_a_in = XSteam('psat_T',Ta_in);
p_a_in = 1.01325; %pression atm
phia_in = Phi_atm;
xa_in = 0.622*phia_in*psat_a_in/(p_a_in-phia_in*psat_a_in);
ha_in = 1.009*(Ta_in+273.15) + xa_in*(2501.6 + 1.854*(Ta_in+273.15));

%% etat a_out : air at the cooling tower outlet

phia_out = Phi_out;
psat_a_out = XSteam('psat_T',Ta_out);
p_a_out = 1.01325; %pression atm
xa_out = 0.622*phia_out*psat_a_out/(p_a_out-phia_out*psat_a_out);
ha_out = 1.009*(Ta_out+273.15) + xa_out*(2501.6 + 1.854*(Ta_out+273.15));

%% etat e3 : water just after  condenser

T_e3 = Tw_out;
h_e3 = Cp_el*(T_e3+273.15); % [kj/kg]
%% etat e2 : water just before condenser

T_e2 = Tw_in;
h_e2 = Cp_el*(T_e2+273.15); % [kj/kg]
%% etat e4 : water from the river (coming between 1 & 2)

T_e4 = Triver;
h_e4 = Cp_el*(T_e4+273.15); % [kj/kg]

%% debits

m_as = P_w/(ha_out-ha_in);
m_ev = m_as*(xa_out-xa_in);
m_e = P_w/(Cp_el*(T_e3-T_e2));

m_e1 = m_e-m_ev;
m_e2 = m_e;
m_e3 = m_e;
m_e4 = m_ev;
%% etat e1 : water outlet of cooling tower

T_e1 = (m_e*Cp_el*(T_e3+273.15)-m_as*(ha_out-ha_in))/(m_e1*Cp_el);
h_e1 = Cp_el*(T_e1+273.15);

MASSFLOW = [m_e,m_ev,m_as];

DAT_WATER = [T_e1       , T_e2       , T_e3       , T_e4;  %[C]
             h_e1       , h_e2       , h_e3       , h_e4;  %[kJ/kg]
             m_e1       , m_e2       , m_e3       , m_e4]; %[kg/s]
DAT_AIR = [Ta_in,   Ta_out;     %[C]
           ha_in,   ha_out;     %[kJ/kg]
           xa_in,   xa_out;     %[-]
           phia_in, phia_out];  %[-]
       
% DAT_AIR(:,1) = [Ta_in ;ha_in;xa_in;phia_in];
% DAT_AIR(:,2) = [Ta_out;ha_out;xa_out;phia_out];

end