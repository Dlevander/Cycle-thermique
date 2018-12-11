function [FIG] = plot_nsout_reheat1(DAT,dat7,d,eta_SiT_HP,eta_SiT_others,nsout)
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
        
        %plot isobare 2 - 3
        linS23 = linspace(DAT(4,2),DAT(4,5),1000);
        linT23 = arrayfun( @(s) XSteam('T_ps',DAT(2,2),s),linS23);
        
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
        
        %plot condens 60-70
        linS6070 = [DAT(4,8) DAT(4,9+nsout)];
        linT6070 = [DAT(1,8) DAT(1,9+nsout)];
        
        %plot passage pompe Pe 7-8
        linS7080 = [DAT(4,9+nsout) DAT(4,10+nsout)];
        linT7080 = [DAT(1,9+nsout) DAT(1,10+nsout)];
        
        %plot chauffage flux principale 8-9d-1
        linS89d_1 = linspace(DAT(4,10+nsout),dat7(4,d),1000);
        linT89d_1 = arrayfun( @(s) XSteam('T_ps',DAT(2,12+nsout+d-1),s),linS89d_1);
        
        %plot bache degaz
        
        %plot Pompe Pb 7d-9d
        linS7d9d = [dat7(4,d) DAT(4,12+nsout+d+1)];
        linT7d9d = [dat7(1,d) DAT(1,12+nsout+d+1)];
        
        %plot chauffage flux principale 9d-1
        linS9d1 = linspace(DAT(4,12+nsout+d-1),DAT(4,1),1000);
        linT9d1 = arrayfun( @(s) XSteam('T_ps',DAT(2,1),s),linS9d1);
        
        %plot Pompe Pa 1-2
        linS12 = [DAT(4,1) DAT(4,2)];
        linT12 = [DAT(1,1) DAT(1,2)]; 
        
        %plot refroidissement isobare et condens soutirage 6-7
        linS67 = zeros(nsout,1000);
        linT67 = zeros(nsout,1000);
        for i=1:nsout
            linS67(i,:) = linspace(DAT(4,8+i),dat7(4,i),1000);
            linT67(i,:) = arrayfun( @(s) XSteam('T_ps',dat7(2,i),s),linS67(i,:));
        end
       
        
        %detente isenthalpique 7-10
        
        
        %plot 
        plot([linS23 linS34 linS45 linS56],[linT23 linT34 linT45 linT56],'r')
        %plot([linS6070 linS78 linS7d9d linS12],[linT6070 linT78 linT7d9d linT12],'k');
        plot([linS6070 linS7080],[linT6070 linT7080],'r');
        
        plot(linS89d_1,linT89d_1,'c')
        plot(linS7d9d,linT7d9d,'g');
        plot(linS9d1,linT9d1,'k');
        plot(linS12,linT12,'m');
     
        plot(DAT(4,:),DAT(1,:),'*')
        %plot(linS67,linT67,'linewidth',2)
        %plot Soutirage
        %plot(linS67,linT67)
        hold off
%% Diagramme h-s
%         FIG(2) = figure;
%         hold on
%         grid on
%         title('Diagramme h-s')
%         xlabel('s [kj/kg K]')
%         ylabel('h [kj/kg]')
%         
%         %cloche de base
%         plot(SL,HL,'-b',SV,HV,'-b')
%         
%         %plot isobare 1 a 3
%         linH13 = arrayfun( @(s) XSteam('h_ps',DAT(2,5),s),linS13);
%         
%         %plot reheat 4-5
%         linH45 = arrayfun( @(s) XSteam('h_ps',DAT(2,6),s),linS45);
%         
%         %plot condens 6-1
%         linH61 = arrayfun( @(s) XSteam('h_ps',DAT(2,1),s),linS67);
%         
%         %plot
%         plot(DAT(4,:),DAT(3,:),'x')
%         plot([linS13 linS34 linS45 linS56 linS61],[linH13 linH34 linH45 linH56 linH61],'r')
end