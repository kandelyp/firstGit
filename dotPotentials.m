%% 4 Quantum Dots
x = linspace(-2,12, 1000);
dotLoc=[2 4 6  8 ]
barLoc=[1 3 5 7 9];
ydot1= -normpdf(x,dotLoc',0.5*ones(1,4)');
ybar1 = normpdf(x,barLoc',0.5*ones(1,5)');
%ydot1([1 end],:)=ydot1([1 end],:)+1
figure(1);clf;%ax1=subplot(1,3,1);
plot(x,0.2*mean(ydot1)+mean(ybar1),'LineWidth',2, 'color','black')
%xticks(dotLoc);
%grid on;
yticks([]);xticks([]);
box 'off'

width=7.2;
height=width*0.25;
set(h,'Units','inches','Position',[1 1 width height])
set([ax1, ax2 , ax3], ...
    'FontSize', 8,...
    'FontName'   , 'Arial',...
    'Box'         , 'on'     , ...
    'TickDir'     , 'in'     , ...
    'TickLength'  , [.03 .03] , ...
    'XMinorTick'  , 'off'      , ...
    'YMinorTick'  , 'off'      , ...
    'LineWidth'   , 1 ,...
    'XTickLabel' , [],...
    'XTick',[2 4 6 8],...
    'ylim',[-0.3 0.3],...    
    'YTick',[] );
%saveas(h,'multiJCartoon','jpg');

%% Load and readout potentials
%% load
 x= linspace(-2,12, 1000);
dotLoc=[2 4 6  8 ]
barLoc=[1 3 5 7 9];
modLoc=[2 8];
ydot1= -normpdf(x,dotLoc',0.5*ones(1,4)');
ybar1 = normpdf(x,barLoc',0.5*ones(1,5)');
ym1   = normpdf(x,modLoc',0.5*ones(1,2)');
figure(1);clf;
plot(x,0.3*mean(ydot1)+mean(ybar1) - 0.5*mean(ym1) ,'LineWidth',2, 'color','black')
%xticks(dotLoc);
%grid on;
yticks([]);xticks([]);
box 'off'
%% separate:
 x= linspace(-2,12, 1000);
dotLoc=[2 4 6  8 ]
barLoc=[1 3 5 7 9];
modLoc=[2 4 6 8];
ydot1= -normpdf(x,dotLoc',0.5*ones(1,4)');
ybar1 = normpdf(x,barLoc',0.5*ones(1,5)');
ym1   = normpdf(x,modLoc',0.5*ones(1,4)');
figure(2);clf;
plot(x,0.4*mean(ydot1)+mean(ybar1) - 0*mean(ym1) ,'LineWidth',2, 'color','black')
yticks([]);xticks([]);
box 'off'

%% manipulate: 12 34
x= linspace(-2,12, 1000);
dotLoc=[2 4 6  8 ]
barLoc=[1 3 5 7 9];
modLoc=[3 7];
ydot1= -normpdf(x,dotLoc',0.5*ones(1,4)');
ybar1 = normpdf(x,barLoc',0.5*ones(1,5)');
ym1   = normpdf(x,modLoc',0.5*ones(1,2)');
figure(3);clf;
plot(x,0.5*mean(ydot1)+mean(ybar1) - 0.2*mean(ym1) ,'LineWidth',2, 'color','black')
yticks([]);xticks([]);
box 'off'
%% manipulate: 12
x= linspace(-2,12, 1000);
dotLoc=[2 4 6  8 ]
barLoc=[1 3 5 7 9];
modLoc=[3 ];
ydot1= -normpdf(x,dotLoc',0.5*ones(1,4)');
ybar1 = normpdf(x,barLoc',0.5*ones(1,5)');
ym1   = normpdf(x,modLoc',0.5*ones(1,2)');
figure(4);clf;
plot(x,0.5*mean(ydot1)+mean(ybar1) - 0.2*mean(ym1) ,'LineWidth',2, 'color','black')
yticks([]);xticks([]);
box 'off'

%% TS load
 x= linspace(-2,12, 1000);
dotLoc=[2 4 6  8 ]
barLoc=[1 3 5 7 9];
modLoc=[2 4 8];
ydot1= -normpdf(x,dotLoc',0.5*ones(1,4)');
ybar1 = normpdf(x,barLoc',0.5*ones(1,5)');
ym1   = normpdf(x,modLoc',0.5*ones(1,3)');
figure(5);clf;
plot(x,0.3*mean(ydot1)+mean(ybar1) - 0.6*mean(ym1) ,'LineWidth',2, 'color','black')
%xticks(dotLoc);
%grid on;
yticks([]);xticks([]);
box 'off'
%% 2 quantum dots
%% 2 quantum dots
x = linspace(-1,7, 1000);
dotLoc=[2 4 ]
barLoc=[1 3 5 ];
ydot1= -normpdf(x,dotLoc',0.5*ones(1,2)');
ybar1 = normpdf(x,barLoc',0.5*ones(1,3)');
%ydot1([1 end],:)=ydot1([1 end],:)+1
figure(100);clf;%ax1=subplot(1,3,1);
plot(x,0.2*mean(ydot1)+mean(ybar1),'LineWidth',2, 'color','black')
%xticks(dotLoc);
%grid on;
yticks([]);xticks([]);
box 'off'
%% 2 dot load
x= linspace(-1,7, 1000);
dotLoc=[2 4];
barLoc=[1 3 5];
modLoc=[2];
ydot1= -normpdf(x,dotLoc',0.5*ones(1,2)');
ybar1 = normpdf(x,barLoc',0.5*ones(1,3)');
ym1   = normpdf(x,modLoc',0.5*ones(1,1)');
figure(101);clf;
plot(x,0.1*mean(ydot1)+mean(ybar1) -0.5* (ym1) ,'LineWidth',2, 'color','black')
%xticks(dotLoc);
%grid on;
yticks([]);xticks([]);
box 'off'

%% 2 dot evolve
x= linspace(-1,7, 1000);
dotLoc=[2 4];
barLoc=[1 3 5];
modLoc=[3];
ydot1= -normpdf(x,dotLoc',0.5*ones(1,2)');
ybar1 = normpdf(x,barLoc',0.5*ones(1,3)');
ym1   = normpdf(x,modLoc',0.5*ones(1,1)');
figure(103);clf;
plot(x,0.3*mean(ydot1)+mean(ybar1) - 0.3* (ym1) ,'LineWidth',2, 'color','black')
%xticks(dotLoc);
%grid on;
yticks([]);xticks([]);
box 'off'

%%
%%

dotLoc=[2.1 3 4 6 7 8]
barLoc=[1 3 5 7 9];
ydot1= -normpdf(x,dotLoc',0.35*ones(1,6)');
ybar1 = normpdf(x,barLoc',[0.35 0.35 0.35 0.35 0.35]');

figure(2);ax3=subplot(1,4,3);
plot(x,mean(ydot1)+mean(ybar1),'LineWidth',2,'Color','b')




dotLoc=[2.0 2.85 7.15  8]
barLoc=[1 3 5 7 9];
ydot1= -normpdf(x,dotLoc',0.35*ones(1,4)');
ybar1 = normpdf(x,barLoc',[0.35 0.35 0.3 0.35 0.35]');
ybar1(1,:)=ybar(1,:)*0.5;
ybar1(3,:)=ybar(3,:)*0.5;
ybar1(5,:)=ybar1(5,:)*0.5;
figure(2);ax2=subplot(1,4,4);
plot(x,mean(ydot1)+mean(ybar1),'LineWidth',2,'Color','b');hold on;
plot(x,0.6+mean(ydot1)+mean(ybar1),'LineWidth',2,'Color','b');h
width=7.2;
height=width*0.25;
set(h,'Units','inches','Position',[1 1 width height])
set([ax1, ax2,ax3,ax4 ], ...
    'FontSize', 8,...
    'FontName'   , 'Arial',...
    'Box'         , 'on'     , ...
    'TickDir'     , 'in'     , ...
    'TickLength'  , [.03 .03] , ...
    'XMinorTick'  , 'off'      , ...
    'YMinorTick'  , 'off'      , ...
    'LineWidth'   , 1 ,...
    'XTickLabel' , [],...
    'XTick',[ ],...
    'ylim',[-0.4 0.3],...    
    'YTick',[],...
    'XLim',[-3 13]);
figure(2);subplot(1,4,4);ylim([-0.4 0.9]);
%saveas(h,'loadReadCartoon','jpg');
%% Tilt and Barrier pulse
dotLoc=[2 4 6 8];
barLoc = [ 1 3 5 7 9];
x=linspace(-1, 11, 500);
ydot= -normpdf(x,dotLoc',0.35*ones(1,4)');
ybar = normpdf(x,barLoc',0.35*ones(1,5)');
h=figure(1);clf;ax1=subplot(1,2,1);
plot(x,mean(ydot)+mean(ybar),'LineWidth',3,'Color','b')
grid on;
title('Symmetric Potential');
%text(-0.1,1.1,'(a)','Units', 'Normalized', 'VerticalAlignment', 'Top','fontweight','bold')


dotLoc=[2.1 2.6 6.3 7 7.7]
barLoc=[1 2.8 4.9 7 9];
ydot1= -normpdf(x,dotLoc',0.35*ones(1,5)');
ybar1 = normpdf(x,barLoc',[0.35 0.35 0.5 0.35 0.35]');
ybar1(3,:)=ybar1(3,:)*1.2;
figure(1);ax2=subplot(1,2,2);
plot(x,mean(ydot1)+mean(ybar1),'LineWidth',3,'Color','b')
%xticks(dotLoc);
grid on;
title('Tilt and Symmetric Exchange');

width=7.2;
height=width*0.25;
set(h,'Units','inches','Position',[1 1 width height])
set([ax1, ax2 ,ax3], ...
    'FontSize', 8,...
    'FontName'   , 'Arial',...
    'Box'         , 'on'     , ...
    'TickDir'     , 'in'     , ...
    'TickLength'  , [.03 .03] , ...
    'XMinorTick'  , 'off'      , ...
    'YMinorTick'  , 'off'      , ...
    'LineWidth'   , 1 ,...
    'XTickLabel' , [],...
    'XTick',[ 2 4 6 8],...
    'ylim',[-0.4 0.3],...    
    'YTick',[] );
%saveas(h,'potential','jpg');

%% Load and readout potentials
dotLoc=[2 4 6 8];
barLoc = [ 1 3 5 7 9];
x=linspace(-3, 13, 500);
ydot= -normpdf(x,dotLoc',0.35*ones(1,4)');
ybar = normpdf(x,barLoc',0.35*ones(1,5)');
h=figure(2);clf;ax1=subplot(1,4,1);
plot(x,mean(ydot)+mean(ybar),'LineWidth',3,'Color','b')

dotLoc=[2.0 2.85 7.15  8]
barLoc=[1 3 5 7 9];
ydot1= -normpdf(x,dotLoc',0.35*ones(1,4)');
ybar1 = normpdf(x,barLoc',[0.35 0.35 0.3 0.35 0.35]');
ybar1(1,:)=ybar(1,:)*0.5;
ybar1(3,:)=ybar(3,:)*0.5;
ybar1(5,:)=ybar1(5,:)*0.5;
figure(2);ax2=subplot(1,4,2);
plot(x,mean(ydot1)+mean(ybar1),'LineWidth',2,'Color','b')




dotLoc=[2.1 3 4 6 7 8]
barLoc=[1 3 5 7 9];
ydot1= -normpdf(x,dotLoc',0.35*ones(1,6)');
ybar1 = normpdf(x,barLoc',[0.35 0.35 0.35 0.35 0.35]');

figure(1);
plot(x,mean(ydot1)+mean(ybar1),'LineWidth',2,'Color','b')




dotLoc=[2.0 2.85 7.15  8]
barLoc=[1 3 5 7 9];
ydot1= -normpdf(x,dotLoc',0.35*ones(1,4)');
ybar1 = normpdf(x,barLoc',[0.35 0.35 0.3 0.35 0.35]');
ybar1(1,:)=ybar(1,:)*0.5;
ybar1(3,:)=ybar(3,:)*0.5;
ybar1(5,:)=ybar1(5,:)*0.5;
figure(3);clf;
plot(x,mean(ydot1)+mean(ybar1),'LineWidth',2,'Color','b');


figure(4);clf;
plot(x,0.6+mean(ydot1)+mean(ybar1),'LineWidth',2,'Color','b');
%%
h=figure(3);
width=7.2;
height=width*0.5;
set(h,'Units','inches','Position',[1 1 width height])
ax=get(gcf,'children');
ax(contains(get(ax,'type'),'colorbar'))=[];

set(ax, ...
    'FontSize', 8,...
    'FontName'   , 'Arial',...
    'Box'         , 'on'     , ...
    'TickDir'     , 'in'     , ...
    'TickLength'  , [.03 .03] , ...
    'XMinorTick'  , 'off'      , ...
    'YMinorTick'  , 'off'      , ...
    'LineWidth'   , 1 ,...
    'XTickLabel' , [],...
    'XTick',[ ],...
    'ylim',[-0.4 0.3],...    
    'YTick',[],...
    'XLim',[-3 13]);
%figure(2);subplot(1,4,4);ylim([-0.4 0.9]);
%saveas(h,'loadReadCartoon','jpg');

%edited in visual studio