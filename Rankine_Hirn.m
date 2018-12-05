clear all 
close all
%% Rankine Hirn
optionsST.nsout = 0;   %  [-] : Number of feed-heating
optionsST.reheat  = 0; %  [-] : Number of reheating
optionsST.T_max = 520; %  [°C] : Maximum steam temperature
optionsST.T_cond_out = 33; %  [°C] : Condenseur cold outlet temperature
optionsST.p3_hp = 40;    %  [bar] : Maximum pressure
%optionsST.drumFlag =1; %[-] :if =1 then drum if =0 => no drum.
optionsST.eta_mec = 0.98; %   [-] : mecanic efficiency of shafts bearings
%is a structure containing combustion data :
optionsST.comb.Tmax = 1900;  % [°C] : maximum combustion temperature  
optionsST.comb.lambda = 1.05;% [-] : air excess
optionsST.comb.x = 0;        % [-] : the ratio O_x/C. Example 0.05 in CH_1.2O_0.05
optionsST.comb.y = 4;        % [-] : the ratio H_y/C. Example 1.2 in CH_1.2O_0.05
optionsST.T_exhaust = 120;   % [°C] : Temperature of exhaust gas out of the chimney
%optionsST.p_3 = 40;          % [-] : High pressure after last reheating
%optionsST.x4                 % [-] : Vapor ratio [gaseous/liquid] (in french : titre)
optionsST.T_0 = 15;          % [°C] : Reference temperature
%optionsST.TpinchSub          % [°C] : Temperature pinch at the subcooler
%optionsST.TpinchEx           % [°C] : Temperature pinch at a heat exchanger
optionsST.TpinchCond = 3;    % [°C] : Temperature pinch at condenser
%   -options.Tdrum            % [°C] : minimal drum temperature
optionsST.eta_SiC = 0.85;    % [-] : Isotrenpic efficiency for compression
%   -option.eta_SiT    [-] : Isotrenpic efficiency for Turbine. It can be a vector of 2 values :
%             	             eta_SiT(1)=eta_SiT_HP,eta_SiT(2)=eta_SiT_others

P_e = 35e3;

[ETA,XMASSFLOW,DATEN,DATEX,DAT,MASSFLOW,COMBUSTION,FIG] = ST(P_e,optionsST,1)
%[ETA,XMASSFLOW,DATEN,DATEX,DAT,MASSFLOW,COMBUSTION,FIG] = ST2(P_e,optionsST,1);

% 
% T = linspace(0.001,400,400);
% SL = zeros(1,400);
% SV = zeros(1,400);
% 
% for i=1:length(T)
%     SL(i) = XSteam('sL_T',T(i));
%     SV(i) = XSteam('sV_T',T(i));
% end
% S = [SL SV];
% T = [T T];
% plot(S,T)
% grid on