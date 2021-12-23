%-----------------------------------
%---------  IMPLICITO, C-N ---------
%-----------------------------------
clear all; clc; close all;
constantes;
%dxx =[10.0, 7.0, 5.0, 4.0, 3.0, 2.0];
%dtt =[7.0, 6.0, 2.5, 1.5, 0.8, 0.4];
dxx =[2.0] ; dtt=[0.4];
%k11=[k/3, k/2, k, k*2, k*3]; 
eL = 60.0;
for m=1:length(dxx(1,:))
tic
%-------------------------------
%--------    Pasos   -----------
%-------------------------------
dx = dxx(m); dy = dx; dt = dtt(m); tol=1e-4;%---- incrementos (MODIFICAR)
%dx = 2.0; dy = dx; dt = 0.4; tol=1e-4;
%k=k11(m);
Adim = 0;%------ 1 = si   /  0 = no
nx = round(lx/dx); ny = round(ly/dy);dx = lx/nx; dy = ly/ny;
tp= round(ts*60*60/dt); dt = ts*60*60/tp; 
%-- Determinación del Nodo Fuente ---
FC_nx = round(posx/dx)+1;% nodo en x
FC_ny = round(posy/dy)+1;% nodo en y
%-- Determinación de Nodos Testeo ---
d_test = [120,150.0,200.0,300.0,500.0,600.0,800.0,1000.0,1200.0,1500.0,1900.0];
for i=1:length(d_test)
    Ft_x = round(d_test(i)/dx)+1;
    Ftest(i) =(nx+1)*(FC_ny -1) +Ft_x;
end
%-------------------------------
%----- Matrices iniciadas ------
%-------------------------------
C = zeros((nx+1)*(ny+1),1); % Mapeo matriz vieja a nueva

%fuente = FC*1000*dt/(86400*dx*dy*h);
fuente = FC*1000/(86400);
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
    C_ant = C;
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
%     for ff=1:length(d_test)
%         anterior=C_ant(Ftest(ff)); actual=CC(Ftest(ff));
%     	if (anterior > actual)||(actual > anterior)
%             tf(m) = it*dt; 
%         end
%     end

end


if m==1
    %fid=fopen('Vol_elem_Impli_gramos_caso01_10.txt','w');
    fid=fopen('Implicito_Diff_caso01_eL_50.txt','w');
elseif m==2
    fid=fopen('Vol_elem_Impli_gramos_caso01_07.txt','w');
    %fid=fopen('Implicito_decai_caso01_2.txt','w');
elseif m==3
    fid=fopen('Vol_elem_Impli_gramos_caso01_05.txt','w');
    %fid=fopen('Implicito_decai_caso01_3.txt','w');
elseif m==4
    fid=fopen('Vol_elem_Impli_gramos_caso01_04.txt','w');
    %fid=fopen('Implicito_decai_caso01_4.txt','w');
elseif m==5
    fid=fopen('Vol_elem_Impli_gramos_caso01_03.txt','w');
    %fid=fopen('Implicito_decai_caso01_5.txt','w');
elseif m==6
    fid=fopen('Vol_elem_Impli_gramos_caso01_02.txt','w');
else
    fid=fopen('Vol_elem_Impli_gramos_caso01_01.txt','w');
end
    
for i=1:(nx+1)*(ny+1)
    fprintf(fid,'%f\n',CC(i));
end
fclose(fid);
clear a12 a2 Bi CC C

tf (m)= toc;
end