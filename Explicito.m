%------------------------------------------
%----------- METODO EXPLICITO  ------------
%------------------------------------------

function [CC2,u_mod,cont_ts]=Explicito(nx,ny,lx,ly,dx,dy,k,dt,tp,FC_nx,FC_ny,FC,c0,u,h0,hv,ts,eL,eT,f,up)
%-------------------------------
%--- Altura del rio variable ---
%-------------------------------
h=zeros(nx+1,1); %rx = h; ry = h;
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
cont_ts = 1; CC2 = zeros((nx+1)*(ny+1), round(ts));

C=zeros((nx+1)*(ny+1),1); % Matriz inicializada
C((nx+1)*(FC_ny -1) +FC_nx)= FC/(86400);%FC*1000/(86400*dx*dy*h(1));
CC = C; 
u_mod = zeros(ny+1,1); % perfil de velocidad del río superficial

%-------------------------------
%------ solución Explícita -----
%-------------------------------
for it =1:tp
    for i =1:(ny+1)     % Recorre de forma Vertical - Figura
        dist = (i-1)*dy;
        
        if up == 1 %-- Distribución Parabólica o constante de la Velocidad
            u_mod(i) = -dist*(dist -ly)*(u*2)/((ly/2)^2);
        else
            u_mod(i) = u;  
        end
        cr = u_mod(i)*dt/dx;  % "cr" Modificado
    
        for j =1:(nx+1) % Recorre de forma Horizontal - Figura 
            ind = (i-1)*(nx+1) + j; % Indice recorrido matriz
            rx = dt*(eL*h(j)*f*u_mod(i))/(dx^2);
            ry = dt*(eT*h(j)*f*u_mod(i))/(dy^2);
            
            xi = 0.5*cr +rx; xd = -0.5*cr +rx;
            yi = ry;         yd = ry; 
            xy = 1 -2*rx -2*ry -k*dt/86400;
        
            if j==1 || j==(nx+1) %- Primer ò Ultima Col. - Figura
                C(ind) = c0;% Cond. Dirichlet
            else 
                if i ==1 %----------- PRIMER FILA - Figura
                    CC(ind) = C(ind-1)*(xi) + C(ind+1)*(xd) +...
                          C(ind +nx+1)*(yi+yd) + C(ind)*(xy);
                     
                elseif i==(ny+1) %--- ULTIMA FILA - Figura
                    CC(ind) = C(ind-1)*(xi) + C(ind+1)*(xd) +...
                         C(ind -(nx+1))*(yi+yd) + C(ind)*(xy);
                     
                elseif i==FC_ny %---- FILA Contaminante - Figura
                    if j==FC_nx 
                        %CC(ind) = FC*1000/(86400*dx*dy*h(j));
                        CC(ind) = FC/(86400);
                    else            %- Col. inter.  - Figura
                        CC(ind) =C(ind-1)*(xi) +C(ind+1)*(xd) +C(ind)*(xy) +...
                         C(ind -(nx+1))*(yi) +C(ind +(nx+1))*(yd);
                    end
                
                else %--------------- FILA INTER. - Figura    
                    CC(ind) =C(ind-1)*(xi) + C(ind+1)*(xd) + C(ind)*(xy) +...
                         C(ind -(nx+1))*(yi) +C(ind +(nx+1))*(yd);
                end
            end
        end
    end
    C = CC;
    %----------------------------------------
    %------ Muestreo para visualización -----
    %----------------------------------------
    if cont_ts*60 == round(it*dt)
    	CC2(:,cont_ts) = CC;
        cont_ts = cont_ts + 1;
    end
end
fid=fopen('Explicito.txt','w');
for i=1:(nx+1)*(ny+1)
    fprintf(fid,'%f\n',CC(i));
end
fclose(fid);
return

