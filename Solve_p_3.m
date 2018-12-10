function [F] = Solve_p_3(p_3,x4,T_cond_out,T_50,eta_SiT_others)
h_50 = XSteam('h_pT',p_3,T_50);
s_50 = XSteam('s_pT',p_3,T_50);
p_60 = XSteam('psat_T',T_cond_out);
h_60s = XSteam('h_ps',p_60,s_50);
h_60 = h_50 - eta_SiT_others*(h_50-h_60s);
x_60 = XSteam('x_ph',p_60,h_60);
F = x4 - x_60;
end