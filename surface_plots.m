%%
clear all;
clc;
%Fig 4B - Minimal regulation network for slope antagonism


hold off;
figure(1); clf;
[s1,s2] = meshgrid(0:.5:10);
mreg=(1+s1+s2)./(2+s1+s2);
h1=surf(s1,s2,mreg);
xlabel('S_1/S_{1,max}');
ylabel('S_2/S_{2,max}');
zlabel('m/m_{max}')
yticks([0 5 10]);
yticklabels({'0','0.5','1'})
xticks([0 5 10]);
xticklabels({'0','0.5','1'})
zlim([0.5,1])
zticks([0.5 0.75 1]);
colormap(parula)
set(gca,'fontsize',25);
set(gca,'fontweight','bold','fontname','arial')
box on;
set(gca, 'LineWidth', 3.5)
saveas(h1,'Fig4B.png')
%%
clear all;
clc;

hold off;
prp=[.5 0 .5];
%dm=2./(2+s1+s2).^2;
s11=10;
s21=0;
dm1=2/(2+s11+s21)^2;
s12=0;
s22=10;
dm2=2/(2+s12+s22)^2;
s13=10;
s23=10;
dm3=2/(2+s13+s23)^2;

lw=5; %linewidth
figure(2); clf;
h1 =bar(0,dm1,'FaceColor','r','LineWidth',lw); hold on; 
bar(1,dm2,'FaceColor','b','LineWidth',lw); hold on; 
bar(2,dm3,'FaceColor',prp,'LineWidth',lw);  
ylim([0 0.018]);
xlim([-1 3]);
xticks([0 1 2 2.5])
xticklabels({'I_1','I_2','I_{1,2}'})
ylabel('slope (\times 10^{-2})')
yticks([0 0.005 0.01 0.015])
ytickangle(0)
yticklabels({'0','0.5','1','1.5'})
xtickangle(0)
set(gca, 'LineWidth', 3.5)
set(gca,'fontsize',30);
x0=10;
y0=10;
width=550;
height=550;
set(gcf,'position',[x0,y0,width,height])
set(gca,'fontweight','bold','fontname','arial')  % Set it to times
saveas(h1,'Fig4B_slope_bar.png')


%%
clear all;
clc;
%Fig 4D - Minimal conversion network for slope antagonism


hold off;
figure(1); clf;
[s1,s2] = meshgrid(0:.5:10);
mreg=(3+2*s1+2*s2)./(2+s1+s2);
h1=surf(s1,s2,mreg/2);
xlabel('S_1/S_{1,max}');
ylabel('S_2/S_{2,max}');
zlabel('m/m_{max}')
yticks([0 5 10]);
yticklabels({'0','0.5','1'})
xticks([0 5 10]);
xticklabels({'0','0.5','1'})
zlim([0.5,1])
zticks([0.5 0.75 1]);
set(gca,'fontsize',25);
set(gca,'fontweight','bold','fontname','arial')
box on;
set(gca, 'LineWidth', 3.5)
saveas(h1,'Fig4D.png')

%%
clear all;
clc;

hold off;
prp=[.5 0 .5];
%dm=2./(2+s1+s2).^2;
s11=10;
s21=0;
dm1=1/(2+s11+s21)^2;
s12=0;
s22=10;
dm2=1/(2+s12+s22)^2;
s13=10;
s23=10;
dm3=1/(2+s13+s23)^2;

lw=5; %linewidth
figure(2); clf;
h1 =bar(0,dm1,'FaceColor','r','LineWidth',lw); hold on; 
bar(1,dm2,'FaceColor','b','LineWidth',lw); hold on; 
bar(2,dm3,'FaceColor',prp,'LineWidth',lw);  
ylim([0 0.01]);
xlim([-1 3]);
xticks([0 1 2 2.5])
xticklabels({'I_1','I_2','I_{1,2}'})
ylabel('slope')
ylabel('slope (\times 10^{-2})')
yticks([0 0.005 0.01])
ytickangle(0)
yticklabels({'0','0.5','1'})
xtickangle(0)
ytickangle(0)
set(gca, 'LineWidth', 3.5)
set(gca,'fontsize',30);
x0=10;
y0=10;
width=550;
height=550;
set(gcf,'position',[x0,y0,width,height])
set(gca,'fontweight','bold','fontname','arial')  % Set it to times
saveas(h1,'Fig4D_slope_bar.png')

%%
clear all;
clc;
%Fig 4F - Minimal binding network for slope antagonism


hold off;
figure(1); clf;
[s1,s2] = meshgrid(0:.5:10);
mreg=0.5*(5+s1+s2-(4+(1+s1+s2).^2).^(0.5));
h1=surf(s1,s2,mreg/2);
xlabel('S_1/S_{1,max}');
ylabel('S_2/S_{2,max}');
zlabel('m/m_{max}')
yticks([0 5 10]);
yticklabels({'0','0.5','1'})
xticks([0 5 10]);
xticklabels({'0','0.5','1'})
zlim([0.5,1])
zticks([0.5 0.75 1]);
set(gca,'fontsize',25);
set(gca,'fontweight','bold','fontname','arial')
box on;
set(gca, 'LineWidth', 3.5)
saveas(h1,'Fig4F.png')


%%
clear all;
clc;

hold off;
prp=[.5 0 .5];
%dm=1-(1+s1+s2)./(4+(1+s1+s2).^2).^0.5;
s11=10;
s21=0;
dm1=1-(1+s11+s21)/(4+(1+s11+s21)^2)^0.5;
s12=0;
s22=10;
dm2=1-(1+s12+s22)/(4+(1+s12+s22)^2)^0.5;
s13=10;
s23=10;
dm3=1-(1+s13+s23)/(4+(1+s13+s23)^2)^0.5;

lw=5; %linewidth
figure(2); clf;
h1 =bar(0,dm1,'FaceColor','r','LineWidth',lw); hold on; 
bar(1,dm2,'FaceColor','b','LineWidth',lw); hold on; 
bar(2,dm3,'FaceColor',prp,'LineWidth',lw);  
ylim([0 0.018]);
xlim([-1 3]);
xticks([0 1 2 2.5])
xticklabels({'I_1','I_2','I_{1,2}'})
ylabel('slope')
ylabel('slope (\times 10^{-2})')
yticks([0 0.005 0.01 0.015])
ytickangle(0)
yticklabels({'0','0.5','1','1.5'})
xtickangle(0)
ytickangle(0)
set(gca, 'LineWidth', 3.5)
set(gca,'fontsize',30);
x0=10;
y0=10;
width=550;
height=550;
set(gcf,'position',[x0,y0,width,height])
set(gca,'fontweight','bold','fontname','arial')  % Set it to times
saveas(h1,'Fig4F_slope_bar.png')
%%
clear all;
clc;

%Fig 5B - Minimal regulation network for value antagonism


hold off;
figure(1); clf;
[s1,s2] = meshgrid(0:.5:10);
mreg=(s1.^2-2*s1.*(s2-1)+(s2+1).^2).^0.5;
h1=surf(s1,s2,mreg/12);
xlabel('S_1/S_{1,max}');
ylabel('S_2/S_{2,max}');
zlabel('m/m_{max}')
yticks([0 5 10]);
yticklabels({'0','0.5','1'})
xticks([0 5 10]);
xticklabels({'0','0.5','1'})
zlim([0.,1])
zticks([0. 0.5 1]);
set(gca,'fontsize',25);
set(gca,'fontweight','bold','fontname','arial')
box on;
set(gca, 'LineWidth', 3.5)
saveas(h1,'Fig5B.png')


%%
clear all;
clc;

hold off;
prp=[.5 0 .5];
%dm=(3+3*s1+s1.^2+3*s2+s1.*s2+s2.^2)./((1+s1).*(1+s2));
s11=10;
s21=0;
dm1=(1/12)*(s11^2-2*s11*(s21-1)+(s21+1)^2)^0.5;
s12=0;
s22=10;
dm2=(1/12)*(s12^2-2*s12*(s22-1)+(s22+1)^2)^0.5;
s13=10;
s23=10;
dm3=(1/12)*(s13^2-2*s13*(s23-1)+(s23+1)^2)^0.5;

lw=5; %linewidth
figure(2); clf;
h1 =bar(0,dm1,'FaceColor','r','LineWidth',lw); hold on; 
bar(1,dm2,'FaceColor','b','LineWidth',lw); hold on; 
bar(2,dm3,'FaceColor',prp,'LineWidth',lw);  
ylim([0 1.2]);
xlim([-1 3]);
xticks([0 1 2 2.5])
xticklabels({'I_1','I_2','I_{1,2}'})
ylabel('slope')
ylabel('Value')
yticks([0 0.5 1])
xtickangle(0)
ytickangle(0)
set(gca, 'LineWidth', 3.5)
set(gca,'fontsize',30);
x0=10;
y0=10;
width=550;
height=550;
set(gcf,'position',[x0,y0,width,height])
set(gca,'fontweight','bold','fontname','arial')  % Set it to times
saveas(h1,'Fig5B_value_bar.png')

%%
clear all;
clc;

%Fig 5D - Minimal conversion network for value antagonism


hold off;
figure(1); clf;
[s1,s2] = meshgrid(0:.5:10);
mreg=(5+2*s1.^2+s1.*(5+s2)+s2.*(5+2*s2))./(3+s1.^2+s1.*(3+s2)+s2.*(3+s2));
h1=surf(s1,s2,mreg/2);
xlabel('S_1/S_{1,max}');
ylabel('S_2/S_{2,max}');
zlabel('m/m_{max}')
yticks([0 5 10]);
yticklabels({'0','0.5','1'})
xticks([0 5 10]);
xticklabels({'0','0.5','1'})
set(gca,'fontsize',20);
zlim([0.7,1])
zticks([0.7 0.85 1]);
set(gca,'fontsize',25);
set(gca,'fontweight','bold','fontname','arial')
box on;
set(gca, 'LineWidth', 3.5)
saveas(h1,'Fig5D.png')

%%
clear all;
clc;

hold off;
prp=[.5 0 .5];
%dm=(5+2*s1.^2+s1.*(5+s2)+s2.*(5+2*s2))./(3+s1.^2+s1.*(3+s2)+s2.*(3+s2));
s11=10;
s21=0;
dm1=0.5*(5+2*s11^2+s11*(5+s21)+s21*(5+2*s21))/(3+s11^2+s11*(3+s21)+s21*(3+s21));
s12=0;
s22=10;
dm2=0.5*(5+2*s12^2+s12*(5+s22)+s22*(5+2*s22))/(3+s12^2+s12*(3+s22)+s22*(3+s22));
s13=10;
s23=10;
dm3=0.5*(5+2*s13^2+s13*(5+s23)+s23*(5+2*s23))/(3+s13^2+s13*(3+s23)+s23*(3+s23));

lw=5; %linewidth
figure(2); clf;
h1 =bar(0,dm1,'FaceColor','r','LineWidth',lw); hold on; 
bar(1,dm2,'FaceColor','b','LineWidth',lw); hold on; 
bar(2,dm3,'FaceColor',prp,'LineWidth',lw);  
ylim([0 1.2]);
xlim([-1 3]);
xticks([0 1 2 2.5])
xticklabels({'I_1','I_2','I_{1,2}'})
ylabel('slope')
ylabel('Value')
yticks([0 0.5 1])
xtickangle(0)
ytickangle(0)
set(gca, 'LineWidth', 3.5)
set(gca,'fontsize',30);
x0=10;
y0=10;
width=550;
height=550;
set(gcf,'position',[x0,y0,width,height])
set(gca,'fontweight','bold','fontname','arial')  % Set it to times
saveas(h1,'Fig5D_value_bar.png')
%%
clear all;
clc;


%Fig 5F - Minimal binding network for value antagonism


hold off;
figure(1); clf;
[s1,s2] = meshgrid(0:.5:10);
mreg=(4*(1+s1)+(1-s1+s2).^2).^(0.5);
h1=surf(s1,s2,mreg/11);
xlabel('S_1/S_{1,max}');
ylabel('S_2/S_{2,max}');
zlabel('m/m_{max}')
yticks([0 5 10]);
yticklabels({'0','0.5','1'})
xticks([0 5 10]);
xticklabels({'0','0.5','1'})
set(gca,'fontsize',20);
zlim([0,1])
zticks([0. 0.5 1]);
set(gca,'fontsize',25);
set(gca,'fontweight','bold','fontname','arial')
box on;
set(gca, 'LineWidth', 3.5)
saveas(h1,'Fig5F.png')

%%
clear all;
clc;

hold off;
prp=[.5 0 .5];
%dm=(4*(1+s1)+(1-s1+s2).^2).^(0.5);
s11=10;
s21=0;
dm1=(4*(1+s11)+(1-s11+s21)^2)^(0.5)/11;
s12=0;
s22=10;
dm2=(4*(1+s12)+(1-s12+s22)^2)^(0.5)/11;
s13=10;
s23=10;
dm3=(4*(1+s13)+(1-s13+s23)^2)^(0.5)/11;
lw=5; %linewidth
figure(2); clf;
h1 =bar(0,dm1,'FaceColor','r','LineWidth',lw); hold on; 
bar(1,dm2,'FaceColor','b','LineWidth',lw); hold on; 
bar(2,dm3,'FaceColor',prp,'LineWidth',lw);  
ylim([0 1.2]);
xlim([-1 3]);
xticks([0 1 2 2.5])
xticklabels({'I_1','I_2','I_{1,2}'})
ylabel('slope')
ylabel('Value')
yticks([0 0.5 1])
xtickangle(0)
ytickangle(0)
set(gca, 'LineWidth', 3.5)
set(gca,'fontsize',30);
x0=10;
y0=10;
width=550;
height=550;
set(gcf,'position',[x0,y0,width,height])
set(gca,'fontweight','bold','fontname','arial')  % Set it to times
saveas(h1,'Fig5F_value_bar.png')
