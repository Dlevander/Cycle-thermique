function [FIG] = plot_GT(DAT,eta_PiT,eta_PiC,x,y,k_cc)

    R_air = 287.1; %[kJ/kg.K] 
    x_a_O2 = 0.21*32/28.96; %fraction massique de O2 dans l'air
    x_a_N2 = 0.79*28/28.96; %fraction massique de N2 dans l'air

    %Compression 1-2
    linP12 = linspace(DAT(2,1),DAT(2,2),100);
    linr12 = linP12./DAT(2,1);
    linT12 = arrayfun(@(r) transf_poly('compression',DAT(1,1),r,eta_PiC,R_air,0,0,0,0),linr12);
    linT12K =linT12+273.15;
%     for i=1:length(linT12K)
%         if linT12K(i)<300
%             linT12Kjanaf(i) = 300;
%         else
%             linT12Kjanaf(i) = linT12K(i);
%         end
%     end
    Cp_12 = 1000*arrayfun(@(t) (x_a_N2*mean(janaf('c','N2',linspace(300,t)))+x_a_O2*mean(janaf('c','O2',linspace(300,t)))),linT12K) ;
    %Cp_12 =(arrayfun(@(t) CP(x_a_O2,0,0,x_a_N2,[300 t]),linT12K)); %[kj/kg/K]
    linS12 = DAT(4,1) + (1-eta_PiC).*Cp_12.*log((linT12K)./(DAT(1,1)+273.15));%[kj/kg/K]
    %Combustion isobare 2-3
    TK_2 = DAT(1,2)+273.15;
    linT23 = linspace(DAT(1,2),DAT(1,3),1000);
    linT23K = linT23+273.15;
    p_3 = DAT(2,2)*k_cc; %pertes de charges dans chambre combustion
    %preallocation
    x_N2 = zeros(1,length(linT23));
    x_O2 = zeros(1,length(linT23));
    x_CO2 = zeros(1,length(linT23));
    x_H2O = zeros(1,length(linT23));
    R_fum = zeros(1,length(linT23));
    Cp_23 = zeros(1,length(linT23));
    linT34= zeros(1,length(linT23));
    Cp_34 = zeros(1,length(linT23));
    %Detente 3-4
    linP34 = linspace(DAT(2,3),DAT(2,4),length(linT23));
    linr34 = linP34./DAT(2,3);
    TK_3 = DAT(1,3)+273.15;
    
    for i=1:length(linT23)
        %Combustion isobare 2-3
        [x_N2(i),x_O2(i),x_CO2(i),x_H2O(i),R_fum(i),~,~,~,~] = combustion(x,y,DAT(1,2),linT23(i),0);
        Cp_23(i) = CP(x_O2(i),x_CO2(i),x_H2O(i),x_N2(i),[TK_2 linT23K(i)]);
        %Detente 3-4
        linT34(i) = transf_poly('detente',DAT(1,3),linr34(i),eta_PiT,R_fum(i),x_CO2(i),x_H2O(i),x_O2(i),x_N2(i));
        linT34K = linT34+273.15;
        Cp_34(i) = CP(x_O2(i),x_CO2(i),x_H2O(i),x_N2(i),[TK_3 linT34K(i)]);
    end
    
    %Combustion isobare 2-3
    linS23 = DAT(4,2)+Cp_23.*log(linT23K./TK_2) - R_fum.*log(p_3/DAT(2,2))*1e-3;
    %Detente 3-4
    linS34 = DAT(4,3) - Cp_34.*log(linT34K./TK_3)* ((1-eta_PiT)/eta_PiT);
    
    %Plot
    close all
    FIG = figure;
    hold on 
    grid on
    plot(DAT(4,:),DAT(1,:),'*');
    plot(linS12,linT12,'r')
    plot(linS23,linT23,'r')
    plot(linS34,linT34,'r')
    title('Diagramme T-s')
    xlabel('s [kj/kg K]')
    ylabel('T [°C]')
    hold off
end