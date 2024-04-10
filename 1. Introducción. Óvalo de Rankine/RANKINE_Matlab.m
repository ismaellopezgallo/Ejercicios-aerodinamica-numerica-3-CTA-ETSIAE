%% 1. Inicialización
close all;
clc;

%% 2. Definición de las funciones de cálculo del campo de velocidades 
   % (La función se encuentra al final del script)
   
%%-----------------------------------------------------------------------%%

%% 3. Cálculo del campo de velocidades, función potencial y función de corriente

% Definimos la magnitud y ángulo de la velocidad incidente
U_inf = 1;
alpha = 0;

% Definimos las coordenadas de los manantiales y sumideros
Q = [20, -20];
x_s = [-4, 4];
z_s = [0, 0];

% Se definen también sus coordenadas en el plano complejo, para el cálculo del potencial.
t_s = x_s + 1i*z_s;
num_sources = length(Q);

% Definimos una malla para las coordenadas x e z
nx = 40;
nz = 40;
x = linspace(-10, 10, nx);
z = linspace(-10, 10, nz);
[Xm, Zm] = meshgrid(x, z);
Tm = Xm + 1i*Zm;

% Calculamos las componentes u y w del campo de velocidades inducido por los manantiales o sumideros
[u, w] = velocity_field(Xm, Zm, Q, x_s, z_s);

% Añadimos el componente de la velocidad incidente
u = u + U_inf * cosd(alpha);
w = w + U_inf * sind(alpha);

% Construimos la función potencial
pot = 0;
for s = 1:num_sources
    pot = pot + Q(s)*log(Tm - t_s(s)) / (2*pi);
end
pot = pot + U_inf * exp(-1i*alpha*pi/180) * Tm;

% De la función potencial extraemos tanto el potencial de velocidades (Phi) 
% como la función de corriente (Psi)
phi = real(pot);
psi = imag(pot);

% ----------------------------------------------------------- %

%% 4. Visualización de resultados
% Visualizamos la función de corriente
figure;

cont = contourf(x, z, psi, 10,'--');
colormap jet                                                      % Pintamos los contornos
hold;
contour(x, z, psi, [0 0], '-k');                                  % Representamos la línea de corriente psi=0    
colorbar('eastoutside')                                           % Pintamos la barra de niveles

N = sqrt(u.^2 + w.^2);                                                
quiver(x, z, u./N, w./N);                                         % Representamos el campo vectorial normalizado

xlabel('x');                                                      % Etiquetamos de forma apropiada, y mostramos en pantalla
ylabel('z');
title('Función de corriente \Psi');
pbaspect([1 1 1])


% Visualizamos el campo de velocidades
figure;

subplot(1,2,1)
cont = contourf(x, z, u, 10);                                 % Pintamos los contornos de U
colormap jet 
hold;
contour(x, z, u, 10, '-k');                                   % Representamos lineas de isocontornos
colorbar('eastoutside')                                       % Pintamos la barra de niveles
xlabel('x');                                                  % Etiquetamos de forma apropiada, y mostramos en pantalla
ylabel('z');
title('Velocidad horizontal \Phi_x');
pbaspect([1 1 1])


subplot(1,2,2)
cont = contourf(x, z, w, 10);                                 % Pintamos los contornos de U
colormap jet 
hold;
contour(x, z, w, 10, '-k');                                   % Representamos lineas de isocontornos
colorbar('eastoutside')                                       % Pintamos la barra de niveles
xlabel('x');                                                  % Etiquetamos de forma apropiada, y mostramos en pantalla
ylabel('z');
title('Velocidad vertical \Phi_y');
pbaspect([1 1 1])



% ----------------------------------------------------------- %

%% Definición de funciones
function [u, w] = velocity_field(Xm, Zm, Q, x_s, z_s)
    % Calcula el campo de velocidades inducido por una serie finita de manantiales y sumideros.
    % Parámetros:
    % Xm, Zm -- coordenadas en el plano xz
    % Q -- magnitud de los manantiales/sumideros
    % x_s, z_s -- coordenadas de los manantiales/sumideros

    % Devuelve:
    % u, w -- componentes del campo de velocidades
    u = zeros(size(Xm));
    w = zeros(size(Zm));
    N = length(Q);
    for i = 1:N
        r = sqrt((Xm - x_s(i)).^2 + (Zm - z_s(i)).^2);
        u = u + Q(i) * (Xm - x_s(i)) ./ (2 * pi * r.^2);
        w = w + Q(i) * (Zm - z_s(i)) ./ (2 * pi * r.^2);
    end

end

