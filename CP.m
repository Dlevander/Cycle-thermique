function CP = CP(x_O2,x_CO2,x_H2O,x_N2,T)
        CpO2 = janaf('c','O2',T);
        CpCO2 = janaf('c','CO2',T);
        CpH2O = janaf('c','H2O',T);
        CpN2 = janaf('c','N2',T);
        Cp_f = x_H2O*CpH2O + x_CO2*CpCO2 + x_O2*CpO2+ x_N2*CpN2;
        CP = mean(Cp_f); % [kj/kg*K]
end