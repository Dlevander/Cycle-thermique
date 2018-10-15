function [e] = exergie(h,s);
T0 = 288.15;
h0 = Xsteam('hL_T',T0);
s0 = Xsteam('sL_T',T0);
e =(h-h0)-T0*(s-s0);
end