clear all; clc; close all;
constantes;
%------------------------------------------
% texto = ['Explicito_caso01_paso02.txt';'Explicito_caso01_paso03.txt';
%          'Explicito_caso01_paso04.txt';'Explicito_caso01_paso07.txt';...
%          'Explicito_caso01_paso10.txt'];
% tit = 'Concentraci�nFinal vs. Posici�n (Explicito)';

texto = ['Implicito_eL015.txt';'Implicito_eL020.txt';'Implicito_eL040.txt'];
tit = 'Difusividad vs. Posici�n (Implicito)';
paso=[2.0, 2.0, 2.0];
figure(1); hold on
%------------------------------------------
for j=1:length(paso)
    A=load(texto(j,:));
    dx = paso(j); dy =dx;
    ny=round(ly/dy); dy =ly/ny;  nx=round(lx/dx); dx =lx/nx;
    %-- Determinaci�n del Nodo Fuente ---
    ii = round(posx/dx)+1;% nodo en x
    jj = round(posy/dy)+1;% nodo en y
    ind = (jj-1)*(nx+1) + ii;
    cont=zeros((nx+1)-ii,1);dist=zeros((nx+1)-ii,1);
    for i = 1:((nx+1)-ii)
        cont(i) = A(ind +(i-1));
        dist(i) = ii*dx + dx*i;
    end
    plot(dist,cont,'LineWidth',1.1); grid on; xlim([0 4100]);
end
legend('0.80','0.85','1.00');
xlabel('Distancia [metros]'); ylabel('Concentraci�n');
title(tit);
hold off;