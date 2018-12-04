function [T_f] = TempComb(LHV, lambda, ma1, h_air, x_O2 , x_CO2 , x_H2O , x_N2)
T_f_init = 1000 ; %pour commencer l'itération
delta= 10;
n=0;
epsilon = 0.5;
CpMoyO2_f = 1000*mean(janaf('c','O2',linspace(300,T_f_init+273,50)));
CpMoyCO2_f = 1000*mean(janaf('c','CO2',linspace(300,T_f_init+273,50)));
CpMoyN2_f = 1000*mean(janaf('c','N2',linspace(300,T_f_init+273,50)));
CpMoyH2O_f = 1000*mean(janaf('c','H2O',linspace(300,T_f_init+273,50)));
CpMoy_f = x_H2O*CpMoyH2O_f + x_CO2*CpMoyCO2_f + x_O2*CpMoyO2_f+ x_N2*CpMoyN2_f;
while (delta>1 && n<50)
    T_f_old = T_f_init;
    CpMoyO2_f = 1000*mean(janaf('c','O2',linspace(300,T_f_init+273,50)));
    CpMoyCO2_f = 1000*mean(janaf('c','CO2',linspace(300,T_f_init+273,50)));
    CpMoyN2_f = 1000*mean(janaf('c','N2',linspace(300,T_f_init+273,50)));
    CpMoyH2O_f = 1000*mean(janaf('c','H2O',linspace(300,T_f_init+273,50)));
    CpMoy_f = x_H2O*CpMoyH2O_f + x_CO2*CpMoyCO2_f + x_O2*CpMoyO2_f+ x_N2*CpMoyN2_f;
    T_f_init = ((1-epsilon)*LHV + lambda*ma1*h_air)/((lambda*ma1+1)*CpMoy_f);
    delta = abs(T_f_init-T_f_old);
    n=n+1;
end
T_f = T_f_init;
