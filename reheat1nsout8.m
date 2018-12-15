clear all 
close all
%% Rankine Hirn reheat 1
% Toujours donnee :
% - nsout
% - reheat
% - T_cond_out
% - eta_mec
% - eta_SiC
% - eta_SiT

optionsST.nsout = 8;   %  [-] : Number of feed-heating
optionsST.reheat  = 1; %  [-] : Number of reheating
optionsST.T_max = 565; %  [°C] : Maximum steam temperature
optionsST.T_cond_out = 33; %  [°C] : Condenseur cold outlet temperature
optionsST.p3_hp = 310;    %  [bar] : Maximum pressure
optionsST.drumFlag = 1; %  [-] : if =1 then drum if =0 => no drum.
options.eta_mec = 0.98;  %[-] : mecanic efficiency of shafts bearings
%optionsST.comb        % is a structure containing combustion data :
%optionsST.comb.Tmax = 1900;%    [Â°C] : maximum combustion temperature
optionsST.comb.lambda = 1.05;%  [-] : air excess
%       -comb.x        [-] : the ratio O_x/C. Example 0.05 in CH_1.2O_0.05
%       -comb.y        [-] : the ratio H_y/C. Example 1.2 in CH_1.2O_0.05
optionsST.T_exhaust = 120; % [C] : Temperature of exhaust gas out of the chimney
optionsST.p_3 = 62; %      [-] : High pressure after last reheating
%optionsST.x4 = 0.89; %        [-] : Vapor ratio [gaseous/liquid] (in french : titre)
%   -options.T_0       [C] : Reference temperature
optionsST.TpinchSub = 5;%[C] : Temperature pinch at the subcooler
optionsST.TpinchEx = 15; %  [C] : Temperature pinch at a heat exchanger
optionsST.TpinchCond = 15; % [C] : Temperature pinch at condenser
optionsST.Tdrum = 140; %    [C] : minimal drum temperature
optionsST.eta_SiC = 0.85;%   [-] : Isotrenpic efficiency for compression
optionsST.eta_SiT = 0.89;   %[-] : Isotrenpic efficiency for Turbine. It can be a vector of 2 values :
%             	             eta_SiT(1)=eta_SiT_HP,eta_SiT(2)=eta_SiT_others

P_e = 288e3; %[kW]

[ETA,XMASSFLOW,DATEN,DATEX,DAT,MASSFLOW,COMBUSTION,FIG] = ST(P_e,optionsST,1)

