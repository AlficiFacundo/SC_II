close all; clear all; clc;
%Definición de parámetros%
polos = [-3 , -2];
ceros = -10;
K = 10;
T2pc = 3;
sobrepaso = 15;
Tm = 0.09;
%Cálculo de ξ, frecuencia normalizada w0 y frecuencia amortiguada wd%
Xi = (-log(sobrepaso/100))/(sqrt((pi^2)+(log(sobrepaso/100)^2)))
w0 = 4/((Xi)*(T2pc))
wd = (w0)*sqrt(1-(Xi)^2)
%Muestras por ciclo de la frec amortiguada
m = (1/Tm)*((2*pi)/wd)
%Polos en Z
r = exp((-Xi)*(w0)*(Tm))
omega = (wd)*(Tm)
%Preparación para Controladores. Creo las FdT continuas y discretas con y sin realimentación del sistema%
Gcont = zpk(ceros, polos, 10)
Gdisc = c2d(Gcont,Tm,'zoh')
%Requerimientos : Xi = 0,5169 ; t2pc=3;
F=feedback(C*Gdisc, 1)
