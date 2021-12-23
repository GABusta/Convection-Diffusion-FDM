clear all; clc; close all;
constantes;
%------------------------------------------

texto1 = 'Unidim_Implic_gr_caso01_02.txt';
texto2='Implicito_gramos_caso01_02.txt';
     
tit = 'ConcentraciónFinal  vs.  posición (Implicito)';
paso=2.0;
figure(1); hold on
%------------------------------------------
%----------- Unidimensional ---------------
%------------------------------------------
A1=load(texto1(1,:));
dx1 = paso;
nx1=round(lx/dx1); dx1 =lx/nx1;
%-- Determinación del Nodo Fuente ---
ii1 = round(posx/dx1)+1;% nodo en x
dist1 = 0.0:dx1:nx1*dx1; 
dist1(1:ii1)=[]; A1(1:ii1)=[];
plot(dist1,A1(:),'LineWidth',1.3);
%------------------------------------------
%----------- Bidimensional ---------------
%------------------------------------------
A2=load(texto2(1,:));
dx2 = paso; dy2 =dx2;
ny2=round(ly/dy2); dy2 =ly/ny2;  nx2=round(lx/dx2); dx2 =lx/nx2;
%-- Determinación del Nodo Fuente ---
ii2 = round(posx/dx2)+1;% nodo en x
jj2 = round(posy/dy2)+1;% nodo en y
ind = (jj2-1)*(nx2+1) + ii2;
cont=zeros((nx2+1)-ii2,1);dist2=zeros((nx2+1)-ii2,1);
for i = 1:((nx2+1)-ii2)
    cont(i) = A2(ind +(i-1));
    dist2(i) = ii2*dx2 + dx2*i;
end
plot(dist2,cont,'LineWidth',1.3); grid on; xlim([0 4100]);

%------------------------------------------
%------------------------------------------
legend('1-D (dx=2m)','2-D (dx=2m)');
xlabel('Distancia [metros]'); ylabel('Concentración [gr/m3]');
title(tit);
hold off;