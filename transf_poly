function [T2]=transf_poly(transf,T1,r,rend,Rstar,x_CO2,x_H2O,x_O2,x_N2)
% T1 en °C
% T2 en °C
% type pour savoir si c'est une compression ou une detente
% r = rapport de compression
% rend= rendement polytropique interne
% Rstar= 8.314/Mm

x_a_O2 = 0.21*32/28.96; %fraction massique de O2 dans l'air
x_a_N2 = 0.79*28/28.96; %fraction massique de N2 dans l'air

T1=T1+273.15;
if T1<300
    Tjanaf = 300;%janaf va pas en dessous de 300K
else
    Tjanaf = T1;
end
T2=Tjanaf;
T2old=0;
n=0;
if strcmp(transf,'compression')
    while(abs(T2-T2old)>0.1 && n<30)
    T2old=T2; 
    Cp = 1000*(x_a_N2*mean(janaf('c','N2',linspace(Tjanaf,T2old,50)))+x_a_O2*mean(janaf('c','O2',linspace(Tjanaf,T2old,50)))) ; %J/kg*K
    %exp = Rstar/(rend*Cp);
    %T2=T1*r^(exp); %T en kelvin p.118 eq 3.22 et 3.19
    %T2=T1*exp(log(r)*Rstar/(rend*Cp)); %eq 3.17
    Cv=Cp-Rstar;
    gamma=Cp/Cv;
    T2=T1*r^((gamma-1)/(gamma*rend));
    n=n+1;
    end
    
elseif strcmp(transf,'detente')
    while(abs(T2-T2old)>0.1 && n<30)
    T2old=T2;
    Cp=(x_N2*mean(janaf('c','N2',linspace(Tjanaf,T2old))) + x_O2*mean(janaf('c','O2',linspace(Tjanaf,T2old))) + x_CO2*mean(janaf('c','CO2',linspace(T1,T2old))) + x_H2O*mean(janaf('c','H2O',linspace(T1,T2old))))*1000; %faire cp moyen entre T2 et T3 ?
    %Cv=Cp-Rstar;
    %gamma=Cp/Cv;
    %T2=T1*r^(rend*(gamma-1)/gamma); %r=p4/p3 = 1/(r_comp*k_cc) eq3.23 et 3.25
    exp=rend*Rstar/Cp;
    T2=T1*r^exp;
    n=n+1;
    end
else
    display('error')
    return  
end
T2=T2-273.15; %sortie en °C

end
