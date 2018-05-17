function [x,y]=CreateGrid
close all;
figure1=figure;
ax=gca;
ax.YTick=-5:5;
ax.XTick=-5:5;
ax.XLim=[-5 5];
ax.YLim=[-5 5];
ax.GridColor='k';
ax.GridAlpha=1;
grid on
hold on
x=linspace(-4.5,4.5,10);
y=linspace(-4.5,4.5,10);
[x,y]=ndgrid(x,y);
x=reshape(x,[1,100]);
y=reshape(y,[1,100]);

%let's denote position of robot by O and goal position by X
scatter(x(100),y(100),500,'filled');
% scatter(x(100),y(100),500,'X');

end
