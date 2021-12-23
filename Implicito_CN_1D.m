%----------------------------------------------------
%---------  IMPLICITO, C-N / Unidimensional ---------
%----------------------------------------------------
clear all; clc; close all;
constantes;
dxx =[10.0, 7.0, 5.0, 4.0, 3.0, 2.0, 1.0];
dtt =[7.0, 6.0, 2.5, 1.5, 0.8, 0.4, 0.1];
for m=1:length(dxx(1,:))
    
dx =dxx(m); dt = dtt(m);%---------- incrementos (MODIFICAR)

nx = round(lx/dx); dx = lx/nx;
tp= round(ts*60*60/dt);                         % Cant. pasos de tiempo
dt = ts*60*60/tp; 
%-- Determinación del Nodo Fuente ---
FC_nx = round(posx/dx)+1;% nodo en x

%-------------------------------
%----- Matrices iniciadas ------
%-------------------------------
C = zeros((nx+1),1); C2 = zeros((nx+1),tp);
fuente = FC*1000/(86400);
C(FC_nx)= fuente; %Cond. Inic.
C(1) = c0; C(nx+1) = c0; 
CC = C;
%-------------------------------
%---- solución Gauss seidel ----
%-------------------------------
cr = u*dt/dx; rx = dt*(eL*h*f*u)/(dx^2); ry = dt*(eT*h*f*u)/(dx^2); 
% a = ctes. iteración actual   , b = ctes. iteración anterior
a = [(1 +rx+ry),(cr*0.25 -rx*0.5),(-cr*0.25 -rx*0.5)];
b = [(1 -rx -ry -k*dt/86400),(-cr*0.25 +rx*0.5),(cr*0.25 +rx*0.5)];     
for it =1:tp
    error = 1.0; cont = 1;
    while (error > tol) && (cont < 10000)             
        for i=2:(nx) % Recorre de forma Horizontal - Figura  
            if i == FC_nx
                CC(i) = fuente;               
            else
                Bi = C(i)*b(1) +C(i+1)*b(2) +C(i-1)*b(3);
                a12 = CC(i-1)*a(3);% n+1
                a23 = C(i+1)*a(2) ;% n 
                CC(i) = (Bi -a12 -a23)/a(1);
            end
        end
        error = max(abs(CC-C));
        cont = cont + 1;
        C = CC; C2(:,it) = CC(:);
    end
end

plot(C2(:,tp))

if m==1
    fid=fopen('Unidim_Implic_gr_caso01_10.txt','w');
elseif m==2
    fid=fopen('Unidim_Implic_gr_caso01_07.txt','w');
elseif m==3
    fid=fopen('Unidim_Implic_gr_caso01_05.txt','w');
elseif m==4
    fid=fopen('Unidim_Implic_gr_caso01_04.txt','w');
elseif m==5
    fid=fopen('Unidim_Implic_gr_caso01_03.txt','w');
elseif m==6
    fid=fopen('Unidim_Implic_gr_caso01_02.txt','w');
else
    fid=fopen('Unidim_Implic_gr_caso01_01.txt','w');
end
    
for i=1:(nx+1)
    fprintf(fid,'%f\n',CC(i));
end
fclose(fid);
clear a12 a2 Bi CC C

end
