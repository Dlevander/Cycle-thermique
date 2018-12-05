function [T_f] = TempComb(LHV, lambda, ma1, h_air, x_O2 , x_CO2 , x_H2O , x_N2)
T_f_init = 1000 ; %pour commencer l'itération
delta= 10;
n=0;
while (delta>1 && n<50)
    T_f_old = T_f_init;
    Cp_f = CP(x_O2,x_CO2,x_H2O,x_N2,linspace(300,T_f_old,100)); %[kj/kg/K]
    T_f_init = (LHV + lambda*ma1*h_air)/((lambda*ma1+1)*Cp_f);
    delta = abs(T_f_init-T_f_old);
    n=n+1;
end
T_f = T_f_init;
end
