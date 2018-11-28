function [FIG] = plotRankineHirn(DAT,eta_SiT_HP,eta_SiT_others)
        T = linspace(0,400,400);
        SL = arrayfun( @(t) XSteam('sL_T',t),T);
        SV = arrayfun( @(t) XSteam('sV_T',t),T);
        HL = arrayfun( @(t) XSteam('hL_T',t),T);
        HV = arrayfun( @(t) XSteam('hV_T',t),T);
        %% Diagramme T-s
        FIG(1) = figure;
        hold on
        %cloche de base
        plot(SL,T,'-b',SV,T,'-b')
        %plot isobare 1 a 5
        linS15 = linspace(DAT(4,1),DAT(4,5),1000);% linspace(DAT(4,2),DAT(4,3),1000) linspace(DAT(4,3),DAT(4,4),1000) linspace(DAT(4,4),DAT(4,5),1000) linspace(DAT(4,5),DAT(4,6),100)];
        linT15 = arrayfun( @(s) XSteam('T_ps',DAT(2,2),s),linS15);
        %plot detente
        linP56 = linspace(DAT(2,5),DAT(2,6),1000);
        linS56s = DAT(4,5)*ones(1,1000); %linspace(DAT(4,5),DAT(4,6),1000);
        linH56s = arrayfun( @(p,s) XSteam('h_ps',p,s),linP56,linS56s);
        linH56 = DAT(3,5) - eta_SiT_HP * eta_SiT_others * (DAT(3,5)-linH56s);
        linS56 = arrayfun( @(p,h) XSteam('s_ph',p,h),linP56,linH56);
        linT56 = arrayfun( @(p,s) XSteam('T_ps',p,s),linP56,linS56s);
        %plot condens
        linS61 = linspace(DAT(4,1),DAT(4,6),1000);
        linT61 = DAT(1,6)*ones(1,1000);
        plot(DAT(4,:),DAT(1,:),'x')
        plot([linS15 linS56 linS61],[linT15 linT56 linT61])
        grid on
        %% Diagramme h-s
        FIG(2) = figure;
        hold on
        grid on
        %cloche de base
        plot(SL,HL,'-b',SV,HV,'-b')
        %plot isobare 1 a 5
        linH15 = arrayfun( @(s) XSteam('h_ps',DAT(2,2),s),linS15);
        %plot condens
        linH61 = arrayfun( @(s) XSteam('h_ps',DAT(2,6),s),linS61);
        plot([linS15 linS56 linS61],[linH15 linH56 linH61],'r')
end