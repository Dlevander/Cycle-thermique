%% Rankine Hirn
optionsST.nsout = 0;   %  [-] : Number of feed-heating
optionsST.reheat  = 0; %  [-] : Number of reheating
optionsST.T_max = 520; %  [�C] : Maximum steam temperature
optionsST.T_cond_out = 33; %  [�C] : Condenseur cold outlet temperature
optionsST.p3_hp = 40;    %  [bar] : Maximum pressure
%   -options.drumFlag  [-] : if =1 then drum if =0 => no drum.
%   -options.eta_mec   [-] : mecanic efficiency of shafts bearings
%   -options.comb is a structure containing combustion data :
%       -comb.Tmax     [�C] : maximum combustion temperature
%       -comb.lambda   [-] : air excess
%       -comb.x        [-] : the ratio O_x/C. Example 0.05 in CH_1.2O_0.05
%       -comb.y        [-] : the ratio H_y/C. Example 1.2 in CH_1.2O_0.05
%   -options.T_exhaust [�C] : Temperature of exhaust gas out of the chimney
%optionsST.p_3 = 40; %      [-] : High pressure after last reheating
%   -options.x4        [-] : Vapor ratio [gaseous/liquid] (in french : titre)
%   -options.T_0       [�C] : Reference temperature
%   -options.TpinchSub [�C] : Temperature pinch at the subcooler
%   -options.TpinchEx  [�C] : Temperature pinch at a heat exchanger
%   -options.TpinchCond[�C] : Temperature pinch at condenser
%   -options.Tdrum     [�C] : minimal drum temperature
optionsST.eta_SiC = 0.85;%   [-] : Isotrenpic efficiency for compression
%   -option.eta_SiT    [-] : Isotrenpic efficiency for Turbine. It can be a vector of 2 values :
%             	             eta_SiT(1)=eta_SiT_HP,eta_SiT(2)=eta_SiT_others

P_e = 35e3;

[~,~,~,~,DAT,~,~,~] = ST(P_e,optionsST,1)

T = linspace(0.001,400,1000);
S = zeros(1,2000);
v = 0;
i = 1;
while v ==0
    s = XSteam('sL_T',T(i));
    if isnan(s)
        v = i;
        s = XSteam('sV_T',T(i));
    end
    S(i) = s;
    i = i+1;
end
for j=v-1:-1:1
    s = XSteam('sV_T',T(j));
    S(i) = s;
    i = i+1;
end
S = S(1:2*v);
T = [linspace(0.001,T(v),v) linspace(T(v),0.001,v)];
plot(S,T)
grid on