% x=-0.5:0.01:0.5;
% y=-1:0.01:1;
% for i=1:length(x)
%     for j=1:length(y)
%         r=sqrt(x(i)^2+(y(j)-0.25)^2);
%         V1(j,i)=100*((r-0.25)^2);
%         r=sqrt(x(i)^2+(y(j)+0.25)^2);
%         V2(j,i)=10*((r-0.25)^2);
%         V=V1.*V2;
%     end
% end
% pcolor(log(V+0.01))
% shading interp
% daspect([1 1 1])
% colorbar
clear all
x=-0.5:0.01:0.5;
y=-0.5:0.01:0.5;
for i=1:length(x)
    for j=1:length(y)
        r=sqrt(0.5*x(i)^2+(y(j))^2);
        V1(j,i)=(min(3*(r-0.25)^2,0.01));
        V2(j,i)=0.05*exp(-1000*x(i)^2);
        if abs(y(j))<0.29
          V3(j,i)=0.01*exp(-1000*(x(i)-0.075)^2);
        else
          V3(j,i)=0;
        end
        if abs(y(j))<0.29
          V4(j,i)=0.01*exp(-1000*(x(i)+0.075)^2);
        else
          V4(j,i)=0;
        end
    end
end
V=max(V1-V3-V4,0);
pcolor(x,y,100*V)
shading interp
daspect([1 1 1])
colorbar