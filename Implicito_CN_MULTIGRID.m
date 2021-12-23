%-----------------------------------------------
%---------  IMPLICITO, C-N / multigrid ---------
%-----------------------------------------------
clear all; clc; close all;
constantes;
A=load('Implicito_caso01_paso02.txt');
dx1 =2.0; dy1 = dx1; nx1 = round(lx/dx1); ny1 = round(ly/dy1);
dx1 = lx/nx1; dy1 = ly/ny1;

%-------------------------------
%--------    Pasos   -----------
%-------------------------------
dx =1.0; dy = dx; dt = 0.1;%---------- incrementos (MODIFICAR)

Adim = 1;%------ 1 = si   /  0 = no
nx = round(lx/dx); ny = round(ly/dy);dx = lx/nx; dy = ly/ny;
tp= round(ts*60*60/dt); dt = ts*60*60/tp; 
%-- Determinación del Nodo Fuente ---
FC_nx = round(posx/dx)+1;% nodo en x
FC_ny = round(posy/dy)+1;% nodo en y

%-------------------------------
%----- Matrices iniciadas ------
%-------------------------------
C = zeros((nx+1)*(ny+1),1); % Mapeo matriz vieja a nueva
for i=1:(ny1+1)    % filas matriz vieja
   for j=1:(nx1+1) % columnas matriz vieja
       ind1=(nx1+1)*(i-1)+j;
       ind =((nx1+1)*2 -1)*2*(i-1) +(ind1*2-1) ;
       C(ind) = A(ind1);
   end
   for j=1:(nx1)
       ind =((nx1+1)*2-1)*2*(i-1) +(j*2);
       C(ind)= (C(ind-1)+C(ind+1))*0.5;
   end
   if i>1
      for j=1:(nx+1)
          ind =((nx1+1)*2-1)*(i-1) + j;
          C(ind)=(C(ind+(nx+1)) + C(ind-(nx+1)))*0.5;
      end
   end
end

%fuente = FC*1000*dt/(86400*dx*dy*h(1));
fuente = FC/(86400);
C((nx+1)*(FC_ny -1) +FC_nx)= fuente;
CC = C;

cr = u*dt/dx;
rx = dt*(eL*h*f*u)/(dx^2);
ry = dt*(eT*h*f*u)/(dy^2);

%-------------------------------
%---- solución Gauss seidel ----
%-------------------------------
for it =1:tp
    error = 1.0; cont = 1;
    while (error > tol) && (cont < 10000)
        for i =1:(ny+1)     % Recorre de forma Vertical - Figura

            for j =1:(nx+1) % Recorre de forma Horizontal - Figura 
                ind =(i-1)*(nx+1) +j; % Indice recorrido matriz
                % a = ctes. iteración actual   , b = ctes. iteración anterior
                a = [(1 +rx +ry),(cr*0.25 -rx*0.5),(-cr*0.25 -rx*0.5),(-0.5*ry),(-0.5*ry)];
                b = [(1 -rx -ry -k*dt/86400),(-cr*0.25 +rx*0.5),(cr*0.25 +rx*0.5),(0.5*ry),(0.5*ry)];
           
                if j==1 || j==(nx+1) %- Primer y ultima Col. - Figura
                    CC(ind) = c0;% Cond. Dirichlet
                else
                    if i ==1        %--- PRIMER FILA - Figura          
                        Bi = C(ind)*b(1) +C(ind+1)*b(2) +C(ind-1)*b(3) +...
                        C(ind +(nx+1))*(b(4) +b(5));    % n
                        a12 = CC(ind-1)*a(3);                             % n+1
                        a23 = C(ind+1)*a(2) +C(ind +(nx+1))*(a(4) + a(5));% n 
                        CC(ind) = (Bi -a12 -a23)/a(1);               

                    elseif i==(ny+1)%--- ULTIMA FILA - Figura
                        Bi = C(ind)*b(1) +C(ind+1)*b(2) +C(ind-1)*b(3) +...
                        C(ind -(nx+1))*(b(4) +b(5));    % n
                        a12 = CC(ind-1)*a(3) +CC(ind -(nx+1))*(a(4) + a(5));% n+1
                        a23 = C(ind+1)*a(2);% n 
                        CC(ind) = (Bi -a12 -a23)/a(1);
            
                    elseif i==FC_ny%--- FILA Contaminante - Figura
                        if j==FC_nx     %- Contaminante - Figura
                            CC(ind) = fuente;
                        else                %- Col. inter.  - Figura
                            Bi = C(ind)*b(1) +C(ind+1)*b(2) +C(ind-1)*b(3) +...
                            C(ind +(nx+1))*b(4) +C(ind -(nx+1))*b(5);    % n
                            a12 = CC(ind-1)*a(3) +CC(ind -(nx+1))*a(5);% n+1
                            a23 = C(ind+1)*a(2)  +C(ind +(nx+1))*a(4);% n 
                            CC(ind) = (Bi -a12 -a23)/a(1);
                        end

                    else            %--- FILA INTER. - Figura
                        Bi = C(ind)*b(1) +C(ind+1)*b(2) +C(ind-1)*b(3) +...
                        C(ind +(nx+1))*b(4) +C(ind -(nx+1))*b(5);    % n
                        a12 = CC(ind-1)*a(3) +CC(ind -(nx+1))*a(5);% n+1
                        a23 = C(ind+1)*a(2) +C(ind +(nx+1))*a(4);% n 
                        CC(ind) = (Bi -a12 -a23)/a(1);
                    end
                end
            end
        end
        error = max(abs(CC-C));
        cont = cont + 1;
        C = CC;
    end
end


if adim==0
    fuente = 1.0;
    fid=fopen('Implicito.txt','w');
else
    fid=fopen('Implicito_Adimen.txt','w');
end

for i=1:(nx+1)*(ny+1)
    fprintf(fid,'%f\n',CC(i)/fuente);
end
fclose(fid);

