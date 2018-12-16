function [ ETA DATEN DATEX DAT MASSFLOW COMBUSTION ] = mainGT(P_e,options,display )
P_e = 225000; %kW
options.k_mec =0.015;
options.T_0 =15;
options.T_ext = 15;
options.r   = 10;
options.k_cc = 0.95;
options.T_3 = 1050;
option.eta_PiC =0.9;
option.eta_PiT =0.9;

[ ETA DATEN DATEX DAT MASSFLOW COMBUSTION ] = GT(P_e,options,0);

end
