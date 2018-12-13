function [FIG] = plot_nsout_reheat1(DAT,dat7,dat10,dat100,dat110,d,eta_SiT_HP,eta_SiT_others,nsout)
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
        cloche = [SL fliplr(SV)];
        Tcloche = [T T];
        plot(cloche,Tcloche,'-b');
        %plot(SL,T,'-b',SV,T,'-b')
        
        %plot Pompe Pa 1-2
        linS12 = [DAT(4,1) DAT(4,2)];
        linT12 = [DAT(1,1) DAT(1,2)];
        
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
        
        %plot detente 5-60
        linP560 = linspace(DAT(2,7),DAT(2,8),1000);
        linS560s = DAT(4,7)*ones(1,1000);
        linH560s = arrayfun( @(p,s) XSteam('h_ps',p,s),linP560,linS560s);
        linH560 = DAT(3,7) - eta_SiT_others * (DAT(3,7)-linH560s);
        linS560 = arrayfun( @(p,h) XSteam('s_ph',p,h),linP560,linH560);
        linT560 = arrayfun( @(p,s) XSteam('T_ps',p,s),linP560,linS560);
        
        %plot condens 60-70
        linS6070 = [DAT(4,8) DAT(4,9+nsout)];
        linT6070 = [DAT(1,8) DAT(1,9+nsout)];
        
        %plot sortie condenseur a entree pompe alim 70-1
        linS701 = DAT(4,9+nsout:11+2*nsout);
        linT701 = DAT(1,9+nsout:11+2*nsout); 
        
        %plot refroidissement isobare et condens Soutirage 6-7
        linS67 = zeros(nsout,1000);
        linT67 = zeros(nsout,1000);
        linH67 = zeros(nsout,1000);
        linS710 = zeros(nsout,2);
        linT710 = zeros(nsout,2);
        linH710 = zeros(nsout,2);
        for i=1:nsout
            linS67(i,:) = linspace(DAT(4,8+i),dat7(4,i),1000);
            linT67(i,:) = arrayfun( @(s) XSteam('T_ps',dat7(2,i),s),linS67(i,:));
            linH67(i,:) = arrayfun( @(s) XSteam('h_ps',dat7(2,i),s),linS67(i,:));
            %plot Soutirage detente isenthalpique 7-10
            linS710(i,:) = [dat7(4,i) dat10(4,i)];
            linT710(i,:) = [dat7(1,i) dat10(1,i)];
            linH710(i,:) = [dat7(3,i) dat10(3,i)];
        end
       %plot detente entre etat 100 et 110
       linS100110 = [dat100(4) dat110(4)];
       linT100110 = [dat100(1) dat110(1)];
        
        %plot 
        plot(linS12,linT12,'r');
        plot([linS23 linS34 linS45 linS560],[linT23 linT34 linT45 linT560],'r');
        plot(linS6070,linT6070,'r');
        plot(linS701,linT701,'r');
        %plot point
        plot(DAT(4,:),DAT(1,:),'*')
        %plot Soutirage
        linSoutS = [linS67 linS710]';
        linSoutT = [linT67 linT710]';
        plot(linSoutS,linSoutT);
        %plot detente vanne apres subcooler
        plot(linS100110,linT100110,'linewidth',3);
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
        
        %plot Pompe Alim 1-2
        linH12 = [DAT(3,1) DAT(3,2)];
        
        %plot isobare 2-3
        linH23 = arrayfun( @(s) XSteam('h_ps',DAT(2,5),s),linS23);
        
        %linH34 detente 3-4 fait ligne 31
        
        %plot reheat 4-5
        linH45 = arrayfun( @(s) XSteam('h_ps',DAT(2,6),s),linS45);
        
        %linH34 detente 5-60 fait ligne 43
        
        %plot condens 60-70
        linH6070 = [DAT(3,8) DAT(3,9+nsout)];

        %plot sortie condenseur a entree pompe alim 70-1
        linH701 = DAT(3,9+nsout:11+2*nsout);
        
        %plot refroidissement isobare et condens Soutirage 6-7 fait ligne 65 
        
        %plot Soutirage detente isenthalpique 7-10 fait ligne 69
        
        %plot detente isenthalpiuee etat 100 et 110
        linH100110 = [dat100(3) dat110(3)];
        
        %plot
        plot(linS12,linH12,'r');
        plot([linS23 linS34 linS45 linS560],[linH23 linH34 linH45 linH560],'r');
        plot(linS6070,linH6070,'r');
        plot(linS701,linH701,'r');
        %plot point
        plot(DAT(4,:),DAT(3,:),'*')
        %plot Soutirage
        linSoutS = [linS67 linS710]';
        linSoutH = [linH67 linH710]';
        plot(linSoutS,linSoutH);
        %plot detente vanne apres subcooler
        plot(linS100110,linH100110,'linewidth',3);
end