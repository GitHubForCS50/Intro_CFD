function [] = createGIF(variable)

% variables
global x y Dt time 

contourf(x,y,variable','EdgeColor','none')
str = sprintf('Time = %0.6f s',time);
title(str)
colorbar
F = getframe(gcf);
[A,map] = rgb2ind(F.cdata,256);
if time == Dt
    imwrite(A,map,'result.gif','gif','LoopCount',Inf,'DelayTime',0);
else
    imwrite(A,map,'result.gif','gif','WriteMode','append','DelayTime',0)
end

end

