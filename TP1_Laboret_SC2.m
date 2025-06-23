close all; clear all; clc;
polos = [-3 , -2];
ceros = -10;
K = 10;
%FdT tiempo continuo sin realimentación
Gcont = zpk(ceros, polos, 10)
Tm = 0.09;
%FdT tiempo discreto sin realimentación
Gdisc = c2d(Gcont,Tm,'zoh')
% Mapa de polos y ceros
pzmap(Gcont)
pzmap(Gdisc)
% 10 veces el período de muestreo.
Gd1 = c2d(Gcont,10*Tm,'zoh')
pzmap(Gd1)
% Respuesta al escalón de ambos sistemas
step(Gcont)
step(Gdisc)
% Sistema con realimentación. Respuesta al escalón
Kp=dcgain(Gdisc)
F=feedback(Gdisc,1)
step(F)
%Respuesta a una rampa
t=0:Tm:100*Tm % genera rampa
lsim(F,t,t)
% Realimentación. Tiempo continuo
pole(F)
zero(F)
Fcont = zpk([0.3917],[0.2653+0.3858i 0.2653-0.3858i],[1.068])
kvec = linspace(-1, 1, 1000)
rlocus(Fcont, kvec)
% Realimentación. Tiempo discreto
Fdisc = c2d(Fcont,Tm,'zoh')
rlocus(Fdisc, kvec)
% 10 veces el período de muestreo
Fd1 = c2d(Fcont,10*Tm,'zoh')
rlocus(Fd1, kvec)
