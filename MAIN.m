%-----------------------------------------------------
%--------------------   TP. Nº1   --------------------
%-----------------------------------------------------
clear all; clc; close all
constantes; % Carga de constantes comunes

%----------------------------------
%------   Método EXPLICITO   ------
%----------------------------------
% tic
% dx1 = 2.0; dy1 = dx1;
% [dt1,nx1,ny1,dx1,dy1]=pasos_dtdx(DT,DL,lx,ly,u,h,hv,eL,eT,f,up,dx1,dy1);
% tp1 = round(ts*60*60/dt1); % Cant. pasos de tiempo 
% dt1 =ts*60*60/tp1;
% %-- Determinación del Nodo Fuente ---
% FC_nx1 = round(posx/dx1)+1;% nodo en x
% FC_ny1 = round(posy/dy1)+1;% nodo en y
% 
% [C1,u_mod1,cont1]=Explicito(nx1,ny1,lx,ly,dx1,dy1,k,dt1,tp1,FC_nx1,FC_ny1,FC,c0,u,h,hv,ts,eL,eT,f,up);
% titulo1 = 'Método Explícito, ts = '; ta1 = 'Explícito 4[km].png'; tb1 = 'Explícito 2[km].png';
% [aa] = imprimirResultados(C1(:,cont1-1),nx1,ny1,1,FC,dx1,dy1,lx,ly,h,titulo1,ts,ta1,tb1);
% %vidName1='Pluma_contaminante_Explícito.avi';
% %[aa]=GrabarVideo(C1,cont1,nx1,ny1,FC,dx1,dy1,lx,ly,h,titulo1,ts,vidName1);
% close all; clear C1 rx1 ry1
% toc

%----------------------------------
%----   Método IMPLICITO C-N   ----
%----------------------------------
tic
dx2 =4.0; dy2 = dx2; dt2 = 1.5;%---------- incrementos (MODIFICAR)
Adim = 0;%------ 1 = si   /  0 = no
[nx2, ny2, dx2, dy2, Pe, Cr]=CN_pasos(lx, ly, dx2, dy2, dt2, u, hv, h, up, eL, eT, f);
tp2= round(ts*60*60/dt2);                         % Cant. pasos de tiempo
dt2 = ts*60*60/tp2; 
%-- Determinación del Nodo Fuente ---
FC_nx2 = round(posx/dx2)+1;% nodo en x
FC_ny2 = round(posy/dy2)+1;% nodo en y

[C2,u_mod2,cont2,C3,fuente]=Implicito_CN(nx2,ny2,lx,ly,dx2,dy2,k,dt2,tp2,FC_nx2,FC_ny2,FC,c0,tol,u,h,hv,ts,eL,eT,f,up,Adim);
titulo2 = 'Método Implicito C-N, ts = '; ta2= 'Implícito 4[km].png'; tb2= 'Implícito 2[km].png';
[aa]=imprimirResultados(C2(:,cont2-1),nx2,ny2,3,FC,dx2,dy2,lx,ly,h,titulo2,ts,ta2,tb2);
%%
%[a]=Tiempos_analisis(C3, dx2, dy2, dt2, tp2, ly, lx, nx2, ny2,tol,FC_ny2,FC_nx2);
%vidName2='Pluma_contaminante_Implícito.avi';
%[aa]=GrabarVideo(C2,cont2,nx2,ny2,FC,dx2,dy2,lx,ly,h,titulo2,ts,vidName2);
%%
%close all; clear C2 C3 rx2 ry2

toc

