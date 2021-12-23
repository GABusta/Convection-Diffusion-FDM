clear all; clc; close all;
tsim=[18*60*60, 4346.0, 1003, 357, 187, 42, 17.8];
dx=[1.0, 2.0, 3.0, 4.0, 5.0, 7.0, 10.0];
plot(dx, tsim,'LineWidth',1.5); grid on;
ylim([0 4700]); xlabel('Pasos dx = dy [m]');ylabel('Tiempo [seg.]');
title('Tiempo de Cálculo vs. Paso espacial');

dt=[0.1, 0.4, 0.8, 1.5, 2.5, 5.0, 7.0]; lx=4000.0; ly=350.0;
for i=1:length(dx(1,:))
    nx = round(lx/dx(i)); ny = round(ly/dx(i));
    tp(i)= round(2.0*60*60/dt(i));
    nodos(i)=(nx+1)*(ny+1);
end




