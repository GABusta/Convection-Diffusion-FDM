%------------------------------------
%---------  IMPLICITO, C-N  ---------
%------------------------------------
function [CC2, u_mod, cont_ts, CC3,fuente]=Implicito_CN(nx,ny,lx,ly,dx,dy,k,dt,tp,FC_nx,FC_ny,FC,c0,tol,u,h0,hv,ts,eL,eT,f,up,adim)
%-------------------------------
%--- Altura del rio variable ---
%-------------------------------
h=zeros(nx+1,1);% rx = h; ry = h;
if hv == 0.0
    h(:) = h0;
else
    for i=1:(nx+1)
        h(i) = h0 + hv*dx*(i-1);
    end
end

%-------------------------------
%----- Matrices iniciadas ------
%-------------------------------
cont_ts = 1; 
CC2 = zeros((nx+1)*(ny+1), round(ts));
%CC3 = zeros((nx+1)*(ny+1), tp);
CC3=0.0;

C = zeros((nx+1)*(ny+1),1);
if adim ==0
    fuente = FC*1000/(86400);
    C((nx+1)*(FC_ny -1) +FC_nx)= fuente;
else
    fuente = FC*dt*1000/(86400*dx*dy*h(1));
    C((nx+1)*(FC_ny -1) +FC_nx)= fuente;
end
CC = C;
u_mod = zeros(ny+1,1); % perfil de velocidad del río superficial

%-------------------------------
%---- solución Gauss seidel ----
%-------------------------------
for it =1:tp
    error = 1.0; cont = 1;
    while (error > tol) && (cont < 10000)
        for i =1:(ny+1)     % Recorre de forma Vertical - Figura
            
            %----------------------------------------------------------
            %-- Distribución Parabólica o constante de la Velocidad ---
            %----------------------------------------------------------
            dist = (i-1)*dy;
            if up == 1 
                u_mod(i) = -dist*(dist -ly)*(u*2)/((ly/2)^2);
            else
                u_mod(i) = u;  
            end        
            cr = u_mod(i)*dt/dx;  % "cr" Modificado
            %----------------------------------------------------------
            
            for j =1:(nx+1) % Recorre de forma Horizontal - Figura 
                ind =(i-1)*(nx+1) +j; % Indice recorrido matriz
                %%%%%%%%%%%%%%%%%%%
                rx = dt*(eL*h(j)*f*u_mod(i))/(dx^2);
                ry = dt*(eT*h(j)*f*u_mod(i))/(dy^2);
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
    %----------------------------------------
    %------ Muestreo para visualización -----
    %----------------------------------------
    %CC3(:,it) = CC(:);
    if cont_ts*60 == round(it*dt)
        CC2(:,cont_ts) = CC;
        cont_ts = cont_ts + 1;
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


return
