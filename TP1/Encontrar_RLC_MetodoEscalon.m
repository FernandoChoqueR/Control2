Estimacion de los 3 puntos equidistantes:
% Cargar datos
datos = readtable("C:\Users\ferok\UNIVERSIDAD\SISTEMAS DE CONTROL 2\PUCHETA\TPs\TP1\Curvas_Medidas_RLC_2025.xls");
datos.Properties.VariableNames = {'Tiempo', 'Corriente', 'V_C','V_e','V_R'};

% Selección de señales
t = datos.Tiempo;
vc = datos.V_C;
I = datos.Corriente;

% Normalizar VC para facilitar el análisis (0 a 1)
vc_norm = (vc - min(vc)) / (max(vc) - min(vc));

% Calcular derivada numérica (aproximación de la pendiente)
dvc = diff(vc_norm) ./ diff(t);
t_diff = t(1:end-1);  % vector de tiempo asociado a la derivada

% Buscar la región donde la derivada supera cierto umbral -> inicio del crecimiento
umbral = 0.01;  % ajustar según la suavidad de la señal
inicio_idx = find(dvc > umbral, 1, 'first');

% Desde el punto de inicio, tomamos tres puntos equiespaciados
h = 45; % separación entre muestras (ajustable o incluso adaptativa)
idx1 = inicio_idx+65;
idx2 = idx1 + h;
idx3 = idx2 + h;

% Validación
if idx3 > length(vc)
    error("Los índices calculados se van fuera del rango. Ajustar 'h'.");
end

% Valores y tiempos
y1 = vc(idx1);
y2 = vc(idx2);
y3 = vc(idx3);

t1 = t(idx1);
t2 = t(idx2);
t3 = t(idx3);

h_tiempo = t2 - t1;

% Método de Chen
ln_ratio = log((y2 - y1) / (y3 - y2));
zeta = 1 / sqrt(1 + (pi / ln_ratio)^2);  %amortiguamiento
omega_d = pi / h_tiempo;
omega_n = omega_d / sqrt(1 - zeta^2);

%Calcular parametros R L C , Hay que suponer un valor? --> C=
C = 5e-6  % en Faradios --> Este valor asignamos nosotros.

L = 1 / (omega_n^2 * C)
R = 2 * zeta * sqrt(L / C)


% Mostrar resultados
fprintf("Tiempos usados: t1 = %.4f s, t2 = %.4f s, t3 = %.4f s\n", t1, t2, t3);
fprintf("Valores VC: y1 = %.4f, y2 = %.4f, y3 = %.4f\n", y1, y2, y3);
fprintf("zeta     = %.4f\n", zeta);
fprintf("omega_d  = %.4f rad/s\n", omega_d);
fprintf("omega_n  = %.4f rad/s\n", omega_n);

% (Opcional) Graficar la selección de puntos en la curva de Vc
figure;
plot(t, vc, 'b'); hold on;
plot([t1 t2 t3], [y1 y2 y3], 'ro', 'MarkerFaceColor','r');
xlabel('Tiempo (s)');
ylabel('V_C (V)');
title('Selección automática de puntos para método de Chen');
legend('V_C(t)', 'Puntos usados');
grid on;

Estimacion de Funcion de transferencia por Chen:

% t1 = 0.0103;  t2 = 0.0120;  t3 = 0.0138;
% y1 = 4.8248; y2 = 11.8104; y3 = 11.9950;

% t1 = 0.0100 ; t2 = 0.0118 ; t3 = 0.0135 ;
% y1 = 0.0000; y2 = 11.6813; y3 = 11.9916 ;

%t1 = 0.0105;  t2 = 0.0115 ; t3 = 0.0125;  %cero cpmplejo
%y1 = 7.7301; y2 = 11.4645; y3 = 11.9328;


%t1 = 0.0105 ; t2 = 0.0112; t3 = 0.0119;     %Grafica--> Fase no minima
%y1 = 7.7301; y2 = 11.0017; y3 = 11.7666;

% t1 = 0.0105 ; t2 = 0.0112 ; t3 = 0.0118 ;
% y1 = 7.7301; y2 = 10.8925; y3 = 11.7127;

% t1 = 0.0105 ; t2 = 0.0110 ; t3 = 0.0114;
% y1 = 7.7301 ; y2 = 10.3224 ; y3 = 11.3409;

u = 12;
y_inf = 12;

sys_est = chen_estimate_tf(t1, y1, t2, y2, t3, y3, u, y_inf);

% Mostrar función de transferencia
disp('Función de transferencia estimada:');
sys_est

% Ceros y polos
z = zero(sys_est);
p = pole(sys_est);

disp('Ceros del sistema estimado:');
disp(z);

disp('Polos del sistema estimado:');
disp(p);

% Graficar respuesta al impulso (útil para ver fase no mínima)
figure;
impulse(sys_est);
title('Respuesta al impulso del sistema estimado');
grid on;

step(u*sys_est);
title('Respuesta al Escalon');
grid on; grid minor;

Metodo de chen para tres puntos equisdistantes:
function sys_est = chen_estimate_tf(t1, y1, t2, y2, t3, y3, u0, y_inf)
% MÉTODO DE CHEN PARA ESTIMACIÓN DE SISTEMA DE 2º ORDEN + CERO
% Entrada:
%   t1, t2, t3  - tiempos de los tres puntos (equiespaciados preferentemente)
%   y1, y2, y3  - valores medidos en esos tiempos
%   u0          - amplitud del escalón aplicado
%   y_inf       - valor final alcanzado por la salida
%
% Salida:
%   sys_est     - función de transferencia estimada (objeto tf)

    % Paso 1: Ganancia estática
    K = y_inf / u0;                 %cuarto punto

    % Paso 2: Calcular k1, k2, k3
    k1 = y1 / (K * u0) - 1;
    k2 = y2 / (K * u0) - 1;
    k3 = y3 / (K * u0) - 1;

    % Paso 3: Discriminante beta_e
    beta_e = 4 * k1^3 * k3 - 3 * k1^2 * k2^2 - 4 * k2^3 + k3^2 + 6 * k1 * k2 * k3;

    % Verificación de validez
    if beta_e < 0
         error('Discriminante negativo. El modelo no es válido con estos puntos.');
    end

    % Paso 4: Calcular alphas
    alpha1 = (k1 * k2 + k3 - sqrt(beta_e)) / (2 * (k1^2 + k2));
    alpha2 = (k1 * k2 + k3 + sqrt(beta_e)) / (2 * (k1^2 + k2));

    % Paso 5: Beta
    beta = (k1 + alpha2) / (alpha1 - alpha2);

    % Paso 6: Constantes de tiempo
    T1 = -t1 / log(alpha1);
    T2 = -t1 / log(alpha2);
    T3 = beta * (T1 - T2) + T1;

    % Paso 7: Función de transferencia
    num = K * [T3 1];
    den = conv([T1 1], [T2 1]);
    sys_est = tf(num, den);
end




