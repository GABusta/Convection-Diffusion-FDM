function [aa]=imprimirResultados(C,nx,ny,it,cc,dx,dy,lx,ly,h,titulo1,ts,ta,tb)
Rc = zeros(ny+1,nx+1);

%-----------------------------------------
%--- Reordenamiento para visualización ---
%-----------------------------------------
for i=1:ny+1
    inic = (nx+1)*(i-1)+1;  
    Rc(i,:) = C(inic:(inic+nx));          
end
%--- máxima concentración ---
%concen = cc*1000/(86400*dx*dy*h); 
concen = cc/(86400); 

%---------------------------------
%--- gráfica de Dominio grande ---
%---------------------------------
figure(it);    
xx = 0.0:dx:lx; yy = 0.0:dy:ly;
[X,Y] = meshgrid(xx,yy);
contourf(X,Y,Rc,40,'LineStyle','none'); colorbar; caxis([0.0 concen]);
xlabel('Largo [metros]'); ylabel('Ancho [metros]');
title([titulo1 num2str(ts) '  horas']);
%set(gca,'DataAspectRatio',[1 1 1]);
saveas(gcf,ta);

%---------------------------------
%--- gráfica de Dominio chico --- 
%---------------------------------
figure(it+1);
contourf(X,Y,Rc,40,'LineStyle','none'); colorbar; caxis([0.0 concen]);
xlabel('Largo [metros]'); ylabel('Ancho [metros]');
title([titulo1 num2str(ts) '  horas']); 
xlim([0 xx(round(length(xx(1,:))/2))]);
%set(gca,'DataAspectRatio',[1 1 1]);
saveas(gcf,tb);


aa = 1;
return

