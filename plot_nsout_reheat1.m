function [FIG] = plot_nsout_reheat1(DAT,dat7,eta_SiT_HP,eta_SiT_others,nsout)
        close all
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
        
        %plot isobare 1 - 3
        linS13 = linspace(DAT(4,1),DAT(4,5),1000);% linspace(DAT(4,2),DAT(4,3),1000) linspace(DAT(4,3),DAT(4,4),1000) linspace(DAT(4,4),DAT(4,5),1000) linspace(DAT(4,5),DAT(4,6),100)];
        linT13 = arrayfun( @(s) XSteam('T_ps',DAT(2,2),s),linS13);
        
        %plot detente 3 - 4
        linP34 = linspace(DAT(2,5),DAT(2,6),1000);
        linS34s = DAT(4,5)*ones(1,1000); %linspace(DAT(4,5),DAT(4,6),1000);
        linH34s = arrayfun( @(p,s) XSteam('h_ps',p,s),linP34,linS34s);
        linH34 = DAT(3,5) - eta_SiT_HP * (DAT(3,5)-linH34s);
        linS34 = arrayfun( @(p,h) XSteam('s_ph',p,h),linP34,linH34);
        linT34 = arrayfun( @(p,s) XSteam('T_ps',p,s),linP34,linS34);
        
        %plot reheat 4-5
        linS45 = linspace(DAT(4,6),DAT(4,7),1000);
        linT45 = arrayfun( @(s) XSteam('T_ps',DAT(2,6),s),linS45);
        
        %plot detente 5-6
        linP56 = linspace(DAT(2,7),DAT(2,8),1000);
        linS56s = DAT(4,7)*ones(1,1000);
        linH56s = arrayfun( @(p,s) XSteam('h_ps',p,s),linP56,linS56s);
        linH56 = DAT(3,7) - eta_SiT_others * (DAT(3,7)-linH56s);
        linS56 = arrayfun( @(p,h) XSteam('s_ph',p,h),linP56,linH56);
        linT56 = arrayfun( @(p,s) XSteam('T_ps',p,s),linP56,linS56);
        
        %plot refroidissement isobare et condens soutirage 6-7
        linS67 = zeros(nsout,1000);
        linT67 = zeros(nsout,1000);
        for i=1:nsout
            linS67(i,:) = linspace(DAT(4,8+i),dat7(4,i),1000);
            linT67(i,:) = arrayfun( @(s) XSteam('T_ps',dat7(2,i),s),linS67(i,:));
        end
       
        
        %detente isenthalpique 7-10
        
        %plot condens
        linS61 = linspace(DAT(4,1),DAT(4,8),1000);
        linT61 = DAT(1,8)*ones(1,1000);
        
        %plot 
        plot([linS13 linS34 linS45 linS56 linS61],[linT13 linT34 linT45 linT56 linT61],'r')
        plot(DAT(4,:),DAT(1,:),'x')
        %plot Soutirage
        plot(linS67',linT67)
        hold off
%% Diagramme h-s
        FIG(2) = figure;
        hold on
        grid on
        title('Diagramme h-s')
        xlabel('s [kj/kg K]')
        ylabel('h [kj/kg]')
        
        %cloche de base
        plot(SL,HL,'-b',SV,HV,'-b')
        
        %plot isobare 1 a 3
        linH13 = arrayfun( @(s) XSteam('h_ps',DAT(2,5),s),linS13);
        
        %plot reheat 4-5
        linH45 = arrayfun( @(s) XSteam('h_ps',DAT(2,6),s),linS45);
        
        %plot condens 6-1
        linH61 = arrayfun( @(s) XSteam('h_ps',DAT(2,1),s),linS61);
        
        %plot
        plot(DAT(4,:),DAT(3,:),'x')
        plot([linS13 linS34 linS45 linS56 linS61],[linH13 linH34 linH45 linH56 linH61],'r')
end