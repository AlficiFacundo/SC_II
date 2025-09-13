% Sistemas de Control II -FCEFyN-UNC 
% Profesor: Dr.Ing. Pucheta, Julian
% Alumno: Alfici, Facundo Ezequiel
% Tp N° 2 - Caso de estudio 3 -
% Ítem 3 
% Calcular un sistema controlador que haga evolucionar al péndulo en el equilibrio estable.
%% CÓMO SABER EL TAMAÑO DE PASO Y EL TIEMPO DE SIMULACIÓN
%Usando Nyquist (extendido), puedo decir que para una frecuencia natural wn, su período
%va a ser al menos 10 o 20 muestras por el menor período natural del
%sistema, es decir :
% Tn = (2*pi)/wn
%Con esto, defino Nyquist planteado lo siguiente: 
% tpaso <= (Tn/20) = pi/(10*wn)
%Para el péndulo se tiene que  wn = sqrt(g/l), que son aproximadamente
%2,47 rad/s. Luego se tiene que Tn = 2.54s y finalmente el tamaño de paso Ts = 2.54/20
%Teniendo al final 0.127s. Cualquier valor menor a este es funcional, pero
%ya que computacionalmente me lo permite, puedo usar 0.01 y estoy sobrado.

%Para el tamaño de simulación puedo plantear que, si uso un controlador de
%tipo DLQR(Es el método que uso en el TP), puedo simular las variables de estado hasta que se
%acerquen a cero, con un eror aceptable menor al 1%.
%De todas maneras, voy a usar un tiempo de simulación de 10s, si veo
%aliasing o inestabilidad, bajo hasta encontrar equilibrio entre precisión
%y costo computacional.

%% Inicialización. Tiempo de muestreo, simulación y de integración.
clc; clear; close all;
color = 'r'; %Pone los plots en rojo. Una pavada.
m = 0.1;
Fricc = 0.1;
l = 1.6;
g = 9.8;
M = 1.5;
Ts = 0.01; % Tamaño de paso
T = 10; %Tiempo de simulación
At = 1e-4; % Tiempo de muestreo
Kmax = T/Ts;

%referencias para posición y ángulo.
ref = 10; %Posición de referencia
flag = 0;

%Condiciones Iniciales
alfa(1) = pi; %Ángulo inicial

%% TIEMPO CONTINUO
%Matrices de lazo abierto en el equilibrio estable
Ac=[0 1 0 0;                      %x1=delta - desplazamiento
0 -Fricc/M -m*g/M 0;              %x2=delta_p
0 0 0 1;                          %x3=phi - angulo
0 -Fricc/(l*M) -g*(m+M)/(l*M) 0]; %x4=phi_p
Bc=[ 0; 1/M; 0; 1/(l*M)];
Cc=[1 0 0 0; 0 0 1 0];
Dc=0;

%Pasaje de tiempo continuo a tiempo discreto
sys_c =ss(Ac,Bc,Cc,Dc);
sys_d=c2d(sys_c,Ts,'zoh');
%% Tiempo Discreto
%Matrices de lazo abierto en el equilibrio estable
A = sys_d.a;
B = sys_d.b;
C = sys_d.c;
%Se usará un Integrador. Se mide desplazamiento
Cref = Cc(1,:);
AA=[A,zeros(4,1);-Cref*A,eye(1)];
BB=[B;-Cref*B];
% %Matrices de controlabilidad y alcanzabilidad. Para ver si se pueden
% operar.
% Mc = [BB AA*BB AA^2*BB AA^3*BB AA^4*BB AA^5*BB AA^6]; Controlabili-dad = rank(Mc);
% Ma = [BB AA*BB AA^2*BB AA^3*BB AA^4*BB AA^5*BB]; Alcanzabilidad = rank(Ma);

%Una buena forma de calcular un controlador con el DLQR es usar el
%Hamiltoniano, aplicar euler y verificar si los movimientos de polos y
%ceros se cumplen como deberían.
%% Cálculo del método DLQR para distema Discreto. Caso 0m a 10m y 10m a 0m

%Parametros DLQR de 0m a 10 m
dd = [.1 .01 .01 .1 0.00001]; %Desplazamiento, Velocidad, Angulo, Velocidad angular, Integrador
QQ = diag(dd);
RR = 2.3e-4;
KK = dlqr(AA,BB,QQ,RR);
K = KK(1:4);
KI = -KK(5);
%Controlador de 10 a 0 m
%Matrices para m*10 (Consigna)
m2 = m*10;
Ac_m2=[0 1 0 0;                    %x1=delta - desplazamiento
0 -Fricc/M -m2*g/M 0;              %x2=delta_p
0 0 0 1;                           %x3=phi - angulo
0 -Fricc/(l*M) -g*(m2+M)/(l*M) 0]; %x4=phi_p
sys_c_m2 = ss(Ac_m2, Bc, Cc, Dc);
sys_d_m2=c2d(sys_c_m2,Ts,'zoh');
%Matrices de lazo abierto en el equilibrio estable
A_m2 = sys_d_m2.a;
B_m2 = sys_d_m2.b;
C_m2 = sys_d_m2.c;
AA_m2=[A_m2,zeros(4,1);-Cref*A_m2,eye(1)]; % para el integrador de m2


% Parametros DLQR de 10m a 0m
dd_m2 = [.1 1e-3 1e-3 .1 0.001]; %Desplazamiento, Velocidad, Ángulo, Velocidad angular, Integrador
QQ_m2 = diag(dd);
RR_m2 = .008;
KK_m2 = dlqr(AA_m2,BB,QQ_m2,RR_m2);
K_m2 = KK_m2(1:4);
KI_m2 = -KK_m2(5);
%Observador
Ao = A';
Bo = Cc';
Co = B';

%Parametros DLRQ del Observador calculado
do = [.01 .01 .01 .0001]; %Desplazamiento, Velocidad, Ángulo, Velocidad angular
Qo = diag(do);
Ro = diag([10000 100000]);
Kko = dlqr(Ao,Bo,Qo,Ro);
Ko=Kko';
t = 0;
x = [0;0;alfa(1);0];
p = x(1);
p_p = x(2);
alfa = x(3);
w = x(4);
tita_pp(1) = 0;
h = Ts/20;
u = [];
i = 1;
u_k(1) = 0;
ref = 10;
flag = 0;
v(1) = 0;
x = [0; 0; alfa(1);0];
x_hat = [0;0;pi;0];
xop=[0 0 pi 0]';

%% Simulo el controlador para KMAX iteraciones

for ki=1:Kmax  % Bucle principal de control con Kmax iteraciones (simulación a nivel de controlador)
    
    y_sal = Cc*x;  % Salida del sistema real (mediciones), con dos componentes
    y_sal_o = Cc*(x_hat - xop);  % Estimación de la salida desde el observador (x_hat)

    v(ki+1) = v(ki) + ref - y_sal(1);  
    % Integrador del error (solo con la salida 1), necesario para eliminar el error en estado estacionario

    % --- LEY DE CONTROL ---
    u1(ki) = -K*(x - xop) + KI*v(ki+1);  
    % Control realimentado de estados + integral del error
    % La línea comentada usa observador en lugar de estados reales

    % --- ZONA MUERTA ---
    zona_muerta = 0;  % Umbral para eliminar comandos pequeños (podría representar histéresis o limitaciones del actuador)
    if(abs(u1(ki)) < zona_muerta)
        u1(ki) = 0;  % Si el control está dentro de la zona muerta, se elimina
    else
        u1(ki) = sign(u1(ki)) * (abs(u1(ki)) - zona_muerta);  
        % Se descuenta la zona muerta del control, manteniendo dirección
    end

    % --- SIMULACIÓN DINÁMICA (Loop Interno) ---
    for kii=1:Ts/h  % Loop interno para integrar el sistema entre dos pasos de control (Ts tiempo de muestreo)
        u(i) = u1(ki);  % Aplicación del control en esta subiteración

        % Dinámica del sistema físico (modelo no lineal del carro + péndulo)
        p_pp = (1/(M + m)) * (u(i) - m*l*tita_pp*cos(alfa(i)) + m*l*w(i)^2*sin(alfa(i)) - Fricc*p_p(i));
        tita_pp = (1/l) * (g*sin(alfa(i)) - p_pp*cos(alfa(i)));  % Aceleración angular del péndulo

        % Integración numérica tipo Euler:
        p_p(i+1) = p_p(i) + h*p_pp;      % Velocidad del carro
        p(i+1) = p(i) + h*p_p(i);        % Posición del carro
        w(i+1) = w(i) + h*tita_pp;       % Velocidad angular del péndulo
        alfa(i+1) = alfa(i) + h*w(i);    % Ángulo del péndulo

        % --- CAMBIO DE FASE EN EL CONTROL ---
        if(p(i) >= 9.99)  % Cuando el carro alcanza 10 metros
            if(flag == 0)  % Solo se ejecuta una vez
                ref = 0;       % Cambia la referencia a 0 (volver al origen)
                m = m*10;      % Se incrementa la masa m (simula cambio de carga)
                flag = 1;      % Marca que el cambio ya se hizo
                K = K_m2;      % Se actualizan las ganancias del controlador
                KI = KI_m2;
            end
        end

        i = i + 1;  % Avance del índice del integrador
    end

    % --- ACTUALIZACIÓN DEL ESTADO ---
    x = [p(i-1); p_p(i-1); alfa(i-1); w(i-1)];  % Estado real extraído de simulación
    x_hat = A*x_hat + B*u1(ki) + Ko*(y_sal - y_sal_o) + xop;  
    % Estimación de estado por observador de tipo Luenberger
end
u(i)=u1(ki);
t=0:h:T;
%Qué hace esta última parte del código?
% Se simula el sistema en tiempo discreto del péndulo, usando un control
% realimentado de estados, con acción integral y observador. La idea
% principal del Loop es simular el movimiento del carro desde 0m a 10m.
% Para luego, durante el trayecto, mantener el péndulo en posición
% vertical, es decir, en un equilibrio inestable.

%Cuando se llega a los 10 metros, se quiere volver al origen y se simula
%una carga de 10 veces la masa del péndulo. Con esto se cambian las
%ganancias del controlador para la nueva condición de peso.

%Funcionalidad :
%1- Mide o estima el estado del sistema.
%2- Calcula la salida y el error respecto a referencia (10m al llegar y 0m
%al volver)
%3- Aplica un control basado en la realimentación de estados (con y sin
%observador) y realiza una acción integral para eliminar error en estado
%estacionario (osilaciones).
%4- Simula con un modelo físico no lineal.
%5- Detecta los cambios de las condiciones como la referencia, es decir, de
%cuando pasa de 10m como referencia (subida) a 0m (bajada), también la el
%cambio de masa y de las ganancias del controlador.
%Actuaizo el estado con el observador.

%% Hago los gráficos de las respuestas
figure(1);
subplot(3,2,1); grid on; hold on;
plot(t,w,color,'LineWidth',1.5);grid on; title('Velocidad angular \omega');
subplot(3,2,2); grid on; hold on;
plot(t,alfa,color,'LineWidth',1.5); title('Ángulo \phi');xlabel('Tiempo');
subplot(3,2,3); grid on; hold on;
plot(t,p,color,'LineWidth',1.5);title('Posición grúa \theta');xlabel('Tiempo');
subplot(3,2,4); grid on; hold on;
plot(t,p_p,color,'LineWidth',1.5);title('Velocidad de grúa \theta_p');
subplot(3,1,3); grid on; hold on;
plot(t,u,color,'LineWidth',1.5);title('Acción de control u');xlabel('Tiempo en Seg.');
%
% figure(2);
% subplot(2,1,1);grid on; hold on;
% plot(alfa,w,color,'LineWidth',1.5);
% title('Ángulo vs Velocidad angular');
% xlabel('Ángulo');ylabel('Velocidad angular');
%
% subplot(2,1,2);grid on; hold on;
% plot(p,p_p,color,'LineWidth',1.5);
% title('Distancia vs velocidad');
% xlabel('Distancia');ylabel('Velocidad');
