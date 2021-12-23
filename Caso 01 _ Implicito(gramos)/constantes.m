%------------------------------------
%------- Constantes TP. N�1  --------
%------------------------------------
h = 10.0;    % [m]   profundidad
u = 0.8;     % [m/s] Velocidad media de la corriente
up = 0;      % Velocidad parab�lica, 0 = no   1 = s�  
eL= 30.0;    % cte - Direcci�n x
eT= 0.3;     % cte - Direcci�n y
f = 0.04;    % Factor de fricci�n
lx= 4000.0;  % [m] Largo del rio - x
ly= 350.0;   % [m] Largo del rio - y
DL= eL*h*f*u;% Coef. Difusi�n "x"
DT= eT*h*f*u;% Coef. Difusi�n "y"


%----------------------------------------
%------- Tiempo de simulaci�n -----------
%----------------------------------------
ts = 2.0;    % [horas]


%----------------------------------------
%-------- Condici�n Dirichlet -----------
%----------------------------------------
c0 = 0.0;


%----------------------------------------
%-------- Tolerancia de c�lculo ---------
%----------------------------------------
tol=6e-5;    % Tolerancia de c�lculo


%----------------------------------------
%---- Descarga contaminante Puntual -----
%----------------------------------------
FC= 60.0;   % [Kg/dia] descarga contaminante
k = 0.1;    % [1/dia]  Tasa decaimiento
posx =100.0; % [m] Ubicaci�n "x", desde nodo extremo izq.
posy = 30.0; % [m] Ubicaci�n "y", desde nodo extremo izq.


%----------------------------------------
%----- Profundidad variable del RIO -----
%----------------------------------------
% Si r�o no tiene pendiente == 0.0
hv = 0.0;
%hv = 0.0005;    % [m]   profundidad variable


