function [nx, ny, dx, dy, Pe, Cr]=CN_pasos(lx, ly, dx, dy, dt, u, hv, h0, up, eL, eT, f)

nx = round(lx/dx); ny = round(ly/dy);
dx = lx/nx; dy = ly/ny;

l1x = [0.0, 3.1]; l1y = [9.0, 0.0];
l2x = [0.0, 6.1]; l2y = [16.5, 0.0];

%---------------------------------------------
%--- Calidad de malla vs. Pe-Cr calculados ---
%---------------------------------------------
if hv == 0.0
    h = h0;
    if up == 1 % Velocidad = Parabólica o constante
        u_mod = -dy*(dy -ly)*(u*2)/((ly/2)^2);
        Pe = u*2*dx/(eL*h*f*u_mod);
        Cr= u_mod*dt/dx;
        rx(1:5) = dt*(eL*h0*f*u)/(dx^2);
    else
        Pe = dx/(eL*h*f);
        Cr = u*dt/dx;
        rx(1:5) = dt*(eL*h0*f*u)/(dx^2);
    end
else
    figure(3);hold on; grid on;
    if up == 1 % Velocidad = Parabólica o constante
        u_mod = -dy*(dy -ly)*(u*2)/((ly/2)^2);
        Cr = u*2*dt/dx;
        Pe = dx*2*u/(eL*h0*f*u_mod);
        rx(1:5) = dt*(eL*h0*f*u)/(dx^2);
    else
        Cr = u*dt/dx;
        Pe = dx/(eL*h0*f);
        rx(1:5) = dt*(eL*h0*f*u)/(dx^2);
    end
end
xx(1:5) = 1.0;

subplot(1,2,1);
plot(Cr,Pe,'xb','LineWidth',3); 
line(l1x, l1y,'Color','red','LineWidth',2); 
line(l2x, l2y,'Color','red','LineWidth',2); 
xlabel('Cr'); ylabel('Pe'); title('calidad mallado, relación Pe vs. Cr');
legend('Pe - Cr (Calculados)');hold off;

subplot(1,2,2);
hold on;plot(xx,'-k','LineWidth',3);plot(rx,'r','LineWidth',3);
ylim([0.0 (rx(1)+0.5)]); 
title('rx límite');legend('rx límite','rx -actual');hold off;

saveas(gcf,'Implicito_calidadMalla.bmp');  close;
    
return

