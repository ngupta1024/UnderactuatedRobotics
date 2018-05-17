function [x,y]=CreateGrid
close all;
figure1=figure;
ax=gca;
ax.YTick=-7:8;
ax.XTick=-7:8;
ax.XLim=[-7 8];
ax.YLim=[-7 8];
ax.GridColor='k';
ax.GridAlpha=1;
grid on
hold on
x=linspace(-6.5,7.5,15);
y=linspace(-6.5,7.5,15);
[x,y]=ndgrid(x,y);
x=reshape(x,[1,225]);
y=reshape(y,[1,225]);

%let's denote position of robot by O and goal position by X
scatter(x(112),y(112),200,'filled');
% scatter(x(100),y(100),500,'X');

end
