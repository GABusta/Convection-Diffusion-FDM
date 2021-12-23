%------------------------------------
%------- Constantes TP. Nº1  --------
%------------------------------------
h = 10.0;    % [m]   profundidad
u = 0.8;     % [m/s] Velocidad media de la corriente
up = 0;      % Velocidad parabólica, 0 = no   1 = sí  
eL= 30.0;    % cte - Dirección x
eT= 0.3;     % cte - Dirección y
f = 0.04;    % Factor de fricción
lx= 4000.0;  % [m] Largo del rio - x
ly= 350.0;   % [m] Largo del rio - y
DL= eL*h*f*u;% Coef. Difusión "x"
DT= eT*h*f*u;% Coef. Difusión "y"


%----------------------------------------
%------- Tiempo de simulación -----------
%----------------------------------------
ts = 2.0;    % [horas]


%----------------------------------------
%-------- Condición Dirichlet -----------
%----------------------------------------
c0 = 0.0;


%----------------------------------------
%-------- Tolerancia de cálculo ---------
%----------------------------------------
tol=6e-5;    % Tolerancia de cálculo


%----------------------------------------
%---- Descarga contaminante Puntual -----
%----------------------------------------
FC= 60.0;   % [Kg/dia] descarga contaminante
k = 0.1;    % [1/dia]  Tasa decaimiento
posx =100.0; % [m] Ubicación "x", desde nodo extremo izq.
posy = 30.0; % [m] Ubicación "y", desde nodo extremo izq.


%----------------------------------------
%----- Profundidad variable del RIO -----
%----------------------------------------
% Si río no tiene pendiente == 0.0
hv = 0.0;
%hv = 0.0005;    % [m]   profundidad variable


