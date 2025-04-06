% --- Código principal ---

% Parámetros PID y simulación
Kp = 2; Ki = 1; Kd = 0.5; T = 0.1;
N = 100;
t_total = N * T;
t = 0:T:t_total - T;

% Señal de error
e = zeros(1, N);
e(20:30) = 0.5;
e(50:70) = 1.5;

% Inicialización de términos
u = zeros(1, N);
P = zeros(1, N);
I = zeros(1, N);
D = zeros(1, N);
s = 0;            % Acumulador para integral
e_prev = 0;       % Error anterior

% Cálculo PID discreto
for k = 1:N
    P(k) = Kp * e(k);
    
    s = s + e(k) * T;   % Integral acumulada
    I(k) = Ki * s;
    
    if k > 1
        D(k) = Kd * (e(k) - e_prev) / T;
    else
        D(k) = 0;
    end
    
    e_prev = e(k);
    u(k) = P(k) + I(k) + D(k);
end

% --- Gráficas ---

figure;

subplot(4,1,1);
plot(t, e, 'b', 'LineWidth', 2); grid on;
ylabel('e[k]'); title('Error');

subplot(4,1,2);
plot(t, P, 'm', 'LineWidth', 2); grid on;
ylabel('Proporcional');

subplot(4,1,3);
plot(t, I, 'g', 'LineWidth', 2); grid on;
ylabel('Integral');

subplot(4,1,4);
plot(t, D, 'c', 'LineWidth', 2); hold on;
plot(t, u, 'r--', 'LineWidth', 2); grid on;
ylabel('Derivativa / PID');
xlabel('Tiempo [s]');
legend('Derivativo', 'Salida PID');
title('Salida del PID y su componente derivativa');