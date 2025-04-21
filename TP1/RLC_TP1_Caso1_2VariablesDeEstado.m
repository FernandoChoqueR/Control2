%RLC - sistema de dos variables de estado:
%lsim(): (solo sit. linelaes)
% Parámetros
R = 220;            %ohm
L = 500e-3;         %Hy
C = 2.2e-6;         %F

t_final=100e-3;  % tiempo final de simulacion [segundos]
t_puntos=100000;  %cantidad de puntos de simulación, paso de simulacion 10us
t_paso=t_final/t_puntos %paso de la base de tiempo

% Matrices del sistema
A = [-R/L, -1/L;
      1/C,   0 ];

B = [1/L;
     0];

C_out = [1 0;     % Corriente
         0 1;     % Tensión en el capacitor
         R 0];    % Tensión en la resistencia Vo(t)

D = [0; 0; 0];

% Sistema en espacio de estado
sys = ss(A, B, C_out, D);

% Tiempo de simulación
t = linspace(0, t_final, t_puntos);  % Simular 1 ms con buena resolución

% Entrada (escalón)
%u = 12 * ones(size(t));  % Fuente de 12V

% Entrada: onda cuadrada entre 12V y -12V, periodo de T=20ms (50Hz) (cambia cada 10ms)
%   T=1/f --> f=1/T
f_cambio = 1/(2*10e-3);  % Frecuencia de cambio cada 10 ms -->50 Hz
u = 12 * square(2*pi*f_cambio*t);  % Señal cuadrada entre +12 y -12 V


% Simulación
[y, t_out, x] = lsim(sys, u, t);

% Graficar
figure;
subplot(4,1,1);
plot(t_out, y(:,1), 'b', 'LineWidth', 2);
title('Corriente i(t)'); ylabel('i [A]'); grid on; grid minor;

subplot(4,1,2);
plot(t_out, y(:,2), 'r', 'LineWidth', 2);
title('Tensión en el capacitor v_C(t)'); ylabel('v_C [V]'); grid on;

subplot(4,1,3);
plot(t_out, u, 'k', 'LineWidth', 2);
title('Tensión en la entrada v_e(t)'); ylabel('v_e [V]'); grid on;

subplot(4,1,4);
plot(t_out, y(:,3), 'g', 'LineWidth', 2);
title('Tensión en la resistencia v_R(t)'); ylabel('v_R [V]');
xlabel('Tiempo [s]'); grid on;
