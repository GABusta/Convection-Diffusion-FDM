clear all; clc; close all;
constantes;
%------------------------------------------
texto = ['Implicito_gramos_caso01_10.txt';'Implicito_gramos_caso01_07.txt';...
        'Implicito_gramos_caso01_05.txt';'Implicito_gramos_caso01_04.txt';...
        'Implicito_gramos_caso01_03.txt';'Implicito_gramos_caso01_02.txt'];
paso=[10.0, 7.0, 5.0, 4.0, 3.0, 2.0];
dt=[7.0, 5.0, 2.5, 1.5, 0.8, 0.4];

%------- Nodos de testeo -------
ntest = 400; d_test = zeros(1,ntest); d_test(1) = 120.0;
for i=1:ntest
    d_test(i) = d_test(i) + 10*i;
end
Ftest = zeros(length(d_test),length(paso));  % Ubicación nodos testeo
concen = zeros(length(d_test),length(paso)); % concen. en nodo testeo
%------------------------------------------
%----------- Bidimensional ---------------
%------------------------------------------
for j=1:length(paso)
    A=load(texto(j,:));
    dx = paso(j); dy =dx;
    ny=round(ly/dy); dy =ly/ny;  nx=round(lx/dx); dx =lx/nx;
    %-- Determinación del Nodo Fuente ---
    ii = round(posx/dx)+1;% nodo en x
    jj = round(posy/dy)+1;% nodo en y
    %-- Determinación de concentración en nodos testeo ---
    for i =1:length(d_test)
        Ft_x = round(d_test(i)/dx)+1;
        Ftest(i,j) =(nx+1)*(jj -1) +Ft_x;
        concen(i,j) = A(Ftest(i,j));
    end
end
%------------------------------------------
%---------------- Error -------------------
%------------------------------------------
err_max = zeros(1,length(paso)-1); delta = zeros(1,length(paso)-1);
err = zeros(1,length(d_test));  
dx2 = zeros(1,length(paso)-1);% dx^2
for j=1:length(paso)-1
    for i=1:length(d_test)
        err(i) = abs(concen(i,6)-concen(i,j));
    end
    err_max(j) = max(err(:));
    delta(j) = err_max(1)/err_max(j);
    dx2(j) = err_max(1)/((paso(j)-paso(6))^2);
end

plot(paso(1:5),delta(:),'b','LineWidth',2); grid on;
xlabel('dx = dy [metros]');
ylabel('crecimiento error [adimensional]');
title('Paso espacial  vs.  crecimiento del error');


