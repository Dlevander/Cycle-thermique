function [x_N2 x_O2 x_CO2 x_H2O R_fum lambda ma1 LHV e_c] = combustion(x,y,T2,T3,lambda);

Mm_O2 = 32*1e-3; %kg/mol
Mm_CO2 = 44*1e-3;
Mm_H2O = 18*1e-3;
Mm_N2 = 28*1e-3;
%% On a besoin du LHV pour trouver lambda
if x == 0 && y == 4 %CH4
    LHV = 50150; % [kj/kg_comb]
    LHV_mol = 803.3 ;% kJ/mol
elseif x == 0 && y == 0 %C
    LHV = 32780; % [kJ/kg]  % Propane
elseif x == 0 && y == 8/3 
    LHV = 46465; % [kJ/kg]
end
ma1 = (32+3.76*28.15)*(1+(y/4))/(12.01+1.008*y); % pouvoir comburivore [kg_air_stoech/kg_comb]
%% Determination de lambda
if lambda==0
cp_CH4 = 2.44*1000; %cp methane a 15C, trouve sur internet J/kg*K
b2 = y/2; %pg 26 cours combu
a2 = 1; % 
b1 = 0; %hypothese
a1 = 0;

hO2_2 = 1000*mean(janaf('c','O2',linspace(300,T2+273.15)))*T2;
hN2_2 = 1000*mean(janaf('c','N2',linspace(300,T2+273.15)))*T2;
hO2_3 = 1000*mean(janaf('c','O2',linspace(300,T3+273.15)))*T3;
hN2_3 = 1000*mean(janaf('c','N2',linspace(300,T3+273.15)))*T3;
hCO2_3 = 1000*mean(janaf('c','CO2',linspace(300,T3+273.15)))*T3;
hH2O_3 = 1000*mean(janaf('c','H2O',linspace(300,T3+273.15)))*T3;

Num = (b2/2+a2)*hO2_3*Mm_O2 - b2*hH2O_3*Mm_H2O - a2*hCO2_3*Mm_CO2 + LHV_mol*1000 + cp_CH4*15*0.016; 
Denom = (hO2_3*Mm_O2+ 3.76*hN2_3*Mm_N2- hO2_2*Mm_O2- 3.76*hN2_2*Mm_N2);
w = Num/Denom; %w, coefficient stoech de l'air
lambda = w/(1+y/4); %X=0 pour nous
end

%%determination des fractions massiques des fumees

if x == 0 && y == 4 % CH4 + 2*lambda*(O2+3.76N2) => 2*(lambda-1)*O2 + CO2 + 2*H2O + 2*lambda*3.76*N2
    %Livre
    LHV = 50150; % [kj/kg_comb]
    e_c = 52215; %kJ/kg fuel exergy. tableau pg 25 du livre
    MmComb = 16 ; %(kg/kmol)
    % Fractions massiques des produits (kg/kg_comb)
    x_O2 = 2*(lambda-1)*32/MmComb;
    x_CO2 = 44/MmComb;
    x_H2O = 2*18/MmComb;
    x_N2 = 2*lambda*3.76*28/MmComb; 
    sum_mol = 2*(lambda-1)+1+2+2*lambda*3.76;
    
elseif x == 0 && y == 0 % C + lambda*(O2+3.76N2) => (lambda-1)*O2 + CO2 + lambda*3.76*N2  p131
    LHV = 32780; % [kJ/kg]
    MmComb = 12; %[kg/kmol]
    e_c = 32400; %[kJ/kg]
    % Fractions massiques des produits [kg/kg_comb]
    x_O2 = (lambda-1)*32/MmComb;
    x_CO2 = 44/MmComb;
    x_H2O = 0;
    x_N2 = lambda*3.76*28/MmComb;
    
elseif x == 0 && y == 8/3 % Propane : C3H8 + 5*lambda*(O2+3.76N2) => 5*(lambda-1)*O2 + 3*CO2 + 4*H2O + 5*lambda*3.76*N2
    
    LHV = 46465; % [kJ/kg]
    MmComb = 44; %[kg/kmol]
    e_c = 49045; %[kJ/kg]
    % Fractions massiques des produits (kg/kg_comb)
    x_O2 = 5*(lambda-1)*32/MmComb;
    x_CO2 = 3*44/MmComb;
    x_H2O = 4*18/MmComb;
    x_N2 = 5*lambda*3.76*28/MmComb;
end
Mm_fum = 0.001*(x_O2*MmComb + x_CO2*MmComb + x_N2*MmComb + x_H2O*MmComb)/sum_mol;
%on change les fractions massiques en (kg/kg_fumee)
Sum = x_O2 + x_CO2 + x_H2O + x_N2 ;
x_O2 = x_O2/Sum;
x_CO2 = x_CO2/Sum;
x_N2 = x_N2/Sum;
x_H2O = x_H2O/Sum;
ma1 = (32+3.76*28.15)*(1+(y/4))/(12.01+1.008*y); % pouvoir comburivore [kg_air_stoech/kg_comb]
R_fum = 8.314/Mm_fum;
end
