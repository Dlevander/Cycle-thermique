function [FIG] = plotRankineHirn(DAT,eta_SiT_HP,eta_SiT_others)
        T = linspace(0,400,400);
        SL = arrayfun( @(t) XSteam('sL_T',t),T);
        SV = arrayfun( @(t) XSteam('sV_T',t),T);
        HL = arrayfun( @(t) XSteam('hL_T',t),T);
        HV = arrayfun( @(t) XSteam('hV_T',t),T);
        %% Diagramme T-s
        FIG(1) = figure;
        hold on
        grid on
        title('Diagramme T-s')
        xlabel('s [kj/kg K]')
        ylabel('T [°C]')
        %cloche de base
        plot(SL,T,'-b',SV,T,'-b')
        %plot isobare 1 a 3
        linS13 = linspace(DAT(4,1),DAT(4,5),1000);% linspace(DAT(4,2),DAT(4,3),1000) linspace(DAT(4,3),DAT(4,4),1000) linspace(DAT(4,4),DAT(4,5),1000) linspace(DAT(4,5),DAT(4,6),100)];
        linT13 = arrayfun( @(s) XSteam('T_ps',DAT(2,2),s),linS13);
        %plot detente
        linP36 = linspace(DAT(2,5),DAT(2,6),1000);
        linS36s = DAT(4,5)*ones(1,1000); %linspace(DAT(4,5),DAT(4,6),1000);
        linH36s = arrayfun( @(p,s) XSteam('h_ps',p,s),linP36,linS36s);
        linH36 = DAT(3,5) - eta_SiT_HP * eta_SiT_others * (DAT(3,5)-linH36s);
        linS36 = arrayfun( @(p,h) XSteam('s_ph',p,h),linP36,linH36);
        linT36 = arrayfun( @(p,s) XSteam('T_ps',p,s),linP36,linS36s);
        %plot condens
        linS61 = linspace(DAT(4,1),DAT(4,6),1000);
        linT61 = DAT(1,6)*ones(1,1000);
 
        %plot
        plot([linS13 linS36 linS61],[linT13 linT36 linT61])
        plot(DAT(4,:),DAT(1,:),'x')
        %% Diagramme h-s
        FIG(2) = figure;
        hold on
        grid on
        title('Diagramme h-s')
        xlabel('h [kj/kg]')
        ylabel('T [°C]')
        %cloche de base
        plot(SL,HL,'-b',SV,HV,'-b')
        %plot isobare 1 a 3
        linH13 = arrayfun( @(s) XSteam('h_ps',DAT(2,2),s),linS13);
        %plot condens
        linH61 = arrayfun( @(s) XSteam('h_ps',DAT(2,6),s),linS61);
        
        %plot
        plot(DAT(4,:),DAT(3,:),'x')
        plot([linS13 linS36 linS61],[linH13 linH36 linH61],'r')
end