function [] = contourfplots()
% variables
global x y u v T eps k f mfu mox mpr time

subplot(4,3,1),quiver(x,y,u',v',0.5),xlabel('Width [m]'),...
    ylabel('Height [m]'),title('Velocity vectors')

subplot(4,3,2)
contourf(x,y,u','EdgeColor','none');
colormap(jet)
colorbar
title('u [m/s]'),xlabel('Width [m]'),ylabel('Height [m]');

subplot(4,3,3)
contourf(x,y,v','EdgeColor','none');
colormap(jet)
colorbar
title('v [m/s]'),xlabel('Width [m]'),ylabel('Height [m]')

subplot(4,3,4)
contourf(x,y,f','EdgeColor','none');
colormap(jet)
colorbar
title('f [-]'),xlabel('Width [m]'),ylabel('Height [m]')

subplot(4,3,5)
contourf(x,y,mfu','EdgeColor','none');
colormap(jet)
colorbar
title('m_{fu} [-]'),xlabel('Width [m]'),ylabel('Height [m]')

subplot(4,3,6)
contourf(x,y,mox','EdgeColor','none');
colormap(jet)
colorbar
title('m_{ox} [-]'),xlabel('Width [m]'),ylabel('Height [m]')

subplot(4,3,7)
contourf(x,y,mpr','EdgeColor','none');
colormap(jet)
colorbar
title('m_{pr} [-]'),xlabel('Width [m]'),ylabel('Height [m]')

subplot(4,3,8)
contourf(x,y,T','EdgeColor','none');
colormap(jet)
colorbar
title('T [K]'),xlabel('Width [m]'),ylabel('Height [m]')

subplot(4,3,9)
contourf(x,y,eps','EdgeColor','none');
colormap(jet)
colorbar
title('epsilon [m^2/s^3]'),xlabel('Width [m]'),ylabel('Height [m]')

subplot(4,3,11)
contourf(x,y,k','EdgeColor','none');
colormap(jet)
colorbar
title('k [m^2/s^2]'),xlabel('Width [m]'),ylabel('Height [m]')

str = sprintf('Time = %0.4f s',time);
title(str);

end