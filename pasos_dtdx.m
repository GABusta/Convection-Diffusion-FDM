%----------------------------------------
%------ Pasos espacial y temporal -------
%----------------------------------------
function [dt,nx,ny,dx,dy]=pasos_dtdx(DT,DL,lx,ly,u,h0,hv,eL,eT,f,up,dx,dy)

ny=round(ly/dy); % Intervalos 
dy =ly/ny;        % paso en -x

nx=round(lx/dx); % Intervalos 
dx =lx/nx;        % paso en -x

%--- paso tiempo "estimado"  rx es el limitante
dt2= 0.5*dx^2/DL;
dt = dt2*0.4; 

%--- rx, ry, cr
rx = dt*DL/(dx^2); 
ry = dt*DT/(dy^2);
cr = u*2*dt/dx;
sqr_cr = sqrt(2*rx);

%----------------------------------------------------------------
%--- GRAFICA - Justificación de la elección del paso temporal ---
%----------------------------------------------------------------
dxx=0.0:0.5:10.0;
dt_DL=(0.5*dxx.^2)/DL;
figure(1);
plot(dxx,dt_DL,'LineWidth',2);
legend('dt-DL'); grid on;
xlabel('Tamaño elemento dx = dx = [m]');
ylabel('dt = [segundos]');
title('dt  vs.  dx=dy  -- valores límites');
ylim([0 6.0]);
saveas(gcf,'PasosExplicito.png'); close;

%---------------------------------------
%--- Difusividades, R y u variables ----
%---------------------------------------
h=zeros(nx+1,1); 
if hv == 0.0
    rx = zeros(ny+1,1); ry = rx; lim = rx ;
    h = h0; yy = 0.0:dy:ly;
    for i =1:ny+1
        if up == 1 % Velocidad = Parabólica o constante
            dist = (i-1)*dy;
            u_mod(i) = -dist*(dist -ly)*(u*2)/((ly/2)^2);% Distribución Parabólica
        else
            u_mod(i) = u;
        end
        rx(i) = dt*(eL*h*f*u_mod(i))/(dx^2);
        ry(i) = dt*(eT*h*f*u_mod(i))/(dy^2);
        cr(i) = u_mod(i)*dt/dx; 
        cr_lim(i) = sqrt((1-2*0.5)^2 +1);
    end
    figure(2); 
    subplot(1,2,1);hold on; grid on;
    plot(yy,rx,'-b','LineWidth',2);
    plot(yy,ry,'-r','LineWidth',2);
    plot(yy,lim+0.5,'--k','LineWidth',2);
    ylim([0.0 0.55]);legend('rx','ry','límite');
    xlabel('Distancia transversal [metros]'); ylabel('valores de " r "');
    hold off;
    if up ==1 %--- Velocidad constante o parabólica
        title('Difusividades/Pasos, h = cte , u = parabólica');
    else
        title('Difusividades/Pasos, h = cte , u = cte.');
    end
    
    subplot(1,2,2);hold on; grid on;
    plot(yy,cr,'b','LineWidth',2);
    plot(yy,cr_lim,'--k','LineWidth',2);
    ylim([0.0 (max(cr_lim)+0.3)]);legend('Cr','Cr - Limite','Cr - Limite');
    ylabel('valores de " Cr "');

    hold off; saveas(gcf,'R_Limites_h(cte).png'); close;
else
    for i =1:ny+1
        if up == 1 % Velocidad = Parabólica o constante
            dist = (i-1)*dy;
            u_mod(i) = -dist*(dist -ly)*(u*2)/((ly/2)^2);
        else
            u_mod(i) = u;
        end
        for j=1:(nx+1)
            h(j) = h0 + hv*dx*(j-1);
            rx(i,j) = dt*(eL*h(j)*f*u_mod(i))/(dx^2);
            ry(i,j) = dt*(eT*h(j)*f*u_mod(i))/(dy^2);
        end  
    end
    xx = 0.0:dx:lx; yy = 0.0:dy:ly;
    [X,Y] = meshgrid(xx,yy);
    figure(2); 
    contourf(X,Y,rx);colorbar; caxis([0.0 0.3]);
    if up ==1 %--- Velocidad constante o parabólica
        title('Rx -> h = variable, u = parabolica');
    else
        title('Rx -> h = variable, u = cte.');
    end
    xlabel('Largo [metros]'); ylabel('Ancho [metros]');
    saveas(gcf,'Rx_Limite_h(var).png'); close;
    
    figure(3); 
    contourf(X,Y,ry);colorbar; caxis([0.0 0.5]);
    if up ==1 %--- Velocidad constante o parabólica
        title('Ry -> h = variable, u = parabolica');
    else
        title('Ry -> h = variable, u = cte.');
    end
    xlabel('Largo [metros]'); ylabel('Ancho [metros]');
    saveas(gcf,'Ry_Limite_h(var).png'); close;
    
    figure(4); hold on; grid on;
    plot(xx,-h,'k','LineWidth',2);
    linea = zeros(length(xx(1,:)),1);
    plot(xx,linea,'--b','LineWidth',2);
    ylim([-h(length(xx(1,:))) 0.5]);
    hold off;
    title('profundidad, h = variable');
    legend('Lecho del río','Línea de Agua');
    xlabel('Largo [metros]'); ylabel('Profundidad [metros]');
    saveas(gcf,'h(var).png'); close;
end

return
