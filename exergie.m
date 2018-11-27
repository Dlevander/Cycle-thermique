function [e] = exergie(h,s);
T0 = 15;
h0 = XSteam('hL_T',T0);
s0 = XSteam('sL_T',T0);
e =(h-h0)-(T0+273.15)*(s-s0);
end