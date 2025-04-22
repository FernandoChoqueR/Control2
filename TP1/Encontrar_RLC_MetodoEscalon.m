datos = readtable("C:\Users\ferok\UNIVERSIDAD\SISTEMAS DE CONTROL 2\PUCHETA\TPs\TP1\Curvas_Medidas_RLC_2025.xls");
datos.Properties.VariableNames = {'Tiempo', 'Corriente', 'V_C','V_e','V_R'}


% cargamos a cada vecotr:
t = datos.Tiempo;
I = datos.Corriente;
vc = datos.V_C;
ve = datos.V_e;
vr = datos.V_R;

% Crear una figura con subplots
figure;

% Subplot 1: Graficar Corriente
subplot(4,1,1);
plot(t, I, 'b', 'LineWidth', 2);
xlabel('Tiempo');
ylabel('Corriente');
title('Corriente vs Tiempo');
grid on;

% Subplot 2: Graficar V_C
subplot(4,1,2);
plot(t, vc , 'r', 'LineWidth', 2);
xlabel('Tiempo');
ylabel('V_C');
title('V_C vs Tiempo');
grid on;

% Subplot 3: Graficar V_e
subplot(4,1,3); 
plot(t, ve,  'k', 'LineWidth', 2);
xlabel('Tiempo');
ylabel('V_e');
title('V_e vs Tiempo');
grid on;

% Subplot 4: Graficar V_R
subplot(4,1,4); 
plot(t, vr,  'g', 'LineWidth', 2);
xlabel('Tiempo');
ylabel('V_R');
title('V_R vs Tiempo');
grid on;

% Normalizar VC para facilitar el análisis (0 a 1)
vc_norm = (vc - min(vc)) / (max(vc) - min(vc));

% Calcular derivada numérica (aproximación de la pendiente)
dvc = diff(vc_norm) ./ diff(t);
t_diff = t(1:end-1);  % vector de tiempo asociado a la derivada

% Buscar la región donde la derivada supera cierto umbral -> inicio del crecimiento
umbral = 0.01;  % ajustar según la suavidad de la señal
inicio_idx = find(dvc > umbral, 1, 'first');

% Desde el punto de inicio, tomamos tres puntos equiespaciados
h = 200; % separación entre muestras (ajustable o incluso adaptativa)
idx1 = inicio_idx;
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
C = 1e-6  % en Faradios --> Este valor asignamos nosotros.

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
