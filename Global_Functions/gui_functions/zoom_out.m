function zoom_out(src,eventdata,h,fac)


q=get(h,'CurrentAxes');

xlim=get(q,'Xlim');
ylim=get(q,'Ylim');
zlim=get(q,'Zlim');

Lx=abs(xlim(2)-xlim(1));
Ly=abs(ylim(2)-ylim(1));
Lz=abs(zlim(2)-zlim(1));

set(q,'Xlim',[xlim(1)-Lx*fac xlim(2)+Lx*fac])
set(q,'Ylim',[ylim(1)-Ly*fac ylim(2)+Ly*fac])
set(q,'Zlim',[zlim(1)-Lz*fac zlim(2)+Lz*fac])





end

