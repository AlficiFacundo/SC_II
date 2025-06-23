close all; clear all; clc;
m = 3;
b = 0.1;
delta = 135;
l = 1;
g = 10;
%%Matrices de estado%%
[A,B,C,D] = linmod('pendulo_mod_tarea',delta*pi/180)
%%Autovalores y rango%%
eig(A)
rank(ctrb(A,B))
%%Matrices ampliadas%%
AA = [[ A ; C ] , (zeros ( 3 , 1 ))]
BA = [ B ; 0 ]
%%Autovalores y rango%%
eig(AA);
rank(ctrb(AA,BA));
%%Diseño del controlador%%
p = -4
K = acker ( AA , BA , [ p p p ] )
k1 = K (1)
k2 = K (2)
k3 = K (3)
eig ( AA - BA*K )
ts=7.5/(-p)
%%Simulación del péndulo con controlador PID%%
sim('pendulo_pid_tarea')
figure(1), plot(tout,yout)
grid on, title('Salida')
xlim([0 5])
figure(2), plot(yout,velocidad) %plano de fase
grid on, title('Plano de fases')
figure(3), plot(tout,torque) % torque total
grid on, title('Torque')
xlim([0 5])
figure(4), plot(tout,-accint) % acción integral
grid on, title('Accion integral')
xlim([0 5])

%%Robustez%%
% correr de nuevo el código de simulación, dibujo y analisis
ymax=max(yout) 
S=(ymax-delta)/delta*100 
erel=(delta-yout)/delta; 
efinal=erel(end) 
ind=find(abs(erel)>.02); 
tss=tout(ind(end))
yte=yout(ind(end))
uf=torque(end)
Intf=-accint(end) 


