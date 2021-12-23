clear all; clc; close all;
constantes;
%------------------------------------------
% texto = ['Explicito_caso01_paso02.txt';'Explicito_caso01_paso03.txt';
%          'Explicito_caso01_paso04.txt';'Explicito_caso01_paso07.txt';...
%          'Explicito_caso01_paso10.txt'];
% tit = 'ConcentraciónFinal vs. Posición (Explicito)';

% texto = ['Implicito_caso01_paso02.txt';'Implicito_caso01_paso03.txt';
%          'Implicito_caso01_paso04.txt';'Implicito_caso01_paso07.txt';...
%          'Implicito_caso01_paso10.txt'];

% texto = ['Vol_elem_Impli_gramos_caso01_02.txt';'Vol_elem_Impli_gramos_caso01_03.txt';...
%     'Vol_elem_Impli_gramos_caso01_04.txt';'Vol_elem_Impli_gramos_caso01_05.txt';...
%     'Vol_elem_Impli_gramos_caso01_07.txt';'Vol_elem_Impli_gramos_caso01_10.txt'];

texto = ['Implicito_Diff_caso01_eL_10.txt';'Implicito_Diff_caso01_eL_15.txt';...
    'Implicito_Diff_caso01_eL_30.txt';'Implicito_Diff_caso01_eL_50.txt';...
    'Implicito_Diff_caso01_eL_60.txt';'Implicito_Diff_caso01_eL_66.txt'];

% texto = ['Unidim_Implic_gr_caso01_01.txt';'Unidim_Implic_gr_caso01_02.txt';...
%     'Unidim_Implic_gr_caso01_03.txt';'Unidim_Implic_gr_caso01_04.txt';...
%     'Unidim_Implic_gr_caso01_05.txt';'Unidim_Implic_gr_caso01_07.txt';...
%     'Unidim_Implic_gr_caso01_10.txt'];
     
%tit = 'ConcentraciónFinal vs. Posición (Implicito)';
tit = 'ConcentraciónFinal (2-D) vs. Posición (Implicito)';
%paso=[2.0, 3.0, 4.0, 5.0, 7.0, 10.0];
%dt=[0.4, 0.8, 1.5, 2.5, 5.0, 7.0];
paso=[2.0, 2.0, 2.0, 2.0, 2.0, 2.0];
dt=[0.4,0.4,0.4,0.4,0.4,0.4];
figure(1); hold on
%------------------------------------------
%----------- Unidimensional ---------------
%------------------------------------------
% for j=1:length(paso)
%     A=load(texto(j,:));
%     dx = paso(j); dy =dx;
%     ny=round(ly/dy); dy =ly/ny;  nx=round(lx/dx); dx =lx/nx;
%     %-- Determinación del Nodo Fuente ---
%     ii = round(posx/dx)+1;% nodo en x
%     dist = 0.0:dx:nx*dx; %A(1:ii)=0.0;
%     dist(1:ii)=[]; A(1:ii)=[];
%     plot(dist,A(:),'LineWidth',1.1); grid on; xlim([0 4100]);
% end
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
    ind = (jj-1)*(nx+1) + ii;
    cont=zeros((nx+1)-ii,1);dist=zeros((nx+1)-ii,1);
    for i = 1:((nx+1)-ii)
        cont(i) = A(ind +(i-1));
        dist(i) = ii*dx + dx*i;
    end
    %fuente = max(cont);
    fuente = 1.0;
    plot(dist,cont/fuente,'LineWidth',1.1); grid on; 
    xlim([0 4100]); %ylim([0 1.1]);
end
%------------------------------------------
%------------------------------------------
legend('eL = 10','eL = 15','eL = 30','eL = 50','eL = 60','eL = 66');
xlabel('Distancia [metros]'); ylabel('Concentración [gramos / m3]');
title(tit);
hold off;


