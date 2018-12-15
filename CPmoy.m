function CPmoy = CPmoy(x_O2,x_CO2,x_H2O,x_N2,T)
if length(T) == 2
    linT = linspace(T(1),T(2),1000);
elseif T ~=300
    linT = linspace(300,T,1000); 
end
    CpO2 = janaf('c','O2',linT);
    CpCO2 = janaf('c','CO2',linT);
    CpH2O = janaf('c','H2O',linT);
    CpN2 = janaf('c','N2',linT);
    Cp_f = x_H2O*CpH2O + x_CO2*CpCO2 + x_O2*CpO2+ x_N2*CpN2;
    %CP = mean(Cp_f); % [kj/kg*K]
    CPmoy = sum((Cp_f./linT)*(linT(2)-linT(1)))/log(linT(length(linT))/linT(1));
end