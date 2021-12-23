function  [aa]=GrabarVideo(C,cont,nx,ny,FC,dx,dy,lx,ly,h,titulo,ts,vidName)

v = VideoWriter(vidName,'Uncompressed AVI');
v.FrameRate = 2;
open(v);
%concen = FC*1000/(86400*dx*dy*h); % máxima concentración
concen = FC/(86400); % máxima concentración  

for i=1:cont-1
    Rc = zeros(ny+1,nx+1);
    %-----------------------------------------
    %--- Reordenamiento para visualización ---
    %-----------------------------------------
    for j=1:ny+1
        inic = (nx+1)*(j-1)+1;  
        Rc(j,:) = C(inic:(inic+nx),i);          
    end
    %-----------------------------------------

    xx = 0.0:dx:lx; yy = 0.0:dy:ly;
    [X,Y] = meshgrid(xx,yy);
    contourf(X,Y,Rc,40,'LineStyle','none'); colorbar; caxis([0.0 concen]);
    xlabel('Largo [metros]'); ylabel('Ancho [metros]');
    title([titulo num2str(ts) '  horas']); 
    xlim([0 xx(round(length(xx(1,:))/2))]);
    %set(gca,'DataAspectRatio',[2 1 1]);
    frame = getframe(gcf);
    writeVideo(v,frame);
    
end

close(v);
aa = 1;
return