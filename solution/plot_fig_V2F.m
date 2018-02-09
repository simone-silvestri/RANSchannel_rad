
%% Plotting figures from solutions

close all
clear all

%% IMPORT DATA

path = 'Figures/';

b_VDN = importdata('V2F-DWX/b');
b_VN  = importdata('V2F-NO/b');
b_SN  = importdata('SA-NO/b');
b_DNS = importdata('DNS/bm');
f_DNS = importdata('DNS/bf');
c_DNS = importdata('DNS/bc');

t01_VDR = importdata('V2F-DWX/t01_rad');
t01_VDN = importdata('V2F-DWX/t01');
t01_VN  = importdata('V2F-NO/t01');
t01_SN  = importdata('SA-NO/t01');
t01_DNS = importdata('DNS/t01m');
f01_DNS = importdata('DNS/t01f');
c01_DNS = importdata('DNS/t01c');

t1_VDR = importdata('V2F-DWX/t1_rad');
t1_VDN = importdata('V2F-DWX/t1');
t1_VN  = importdata('V2F-NO/t1');
t1_SN  = importdata('SA-NO/t1');
t1_DNS = importdata('DNS/t1m');
f1_DNS = importdata('DNS/t1f');
c1_DNS = importdata('DNS/t1c');

t5_VDR = importdata('V2F-DWX/t5_rad');
t5_VDN = importdata('V2F-DWX/t5');
t5_VN  = importdata('V2F-NO/t5');
t5_SN  = importdata('SA-NO/t5');
t5_DNS = importdata('DNS/t5m');
f5_DNS = importdata('DNS/t5f');
c5_DNS = importdata('DNS/t5c');

t10_VDR = importdata('V2F-DWX/t10_rad');
t10_VDN = importdata('V2F-DWX/t10');
t10_VN  = importdata('V2F-NO/t10');
t10_SN  = importdata('SA-NO/t10');
t10_DNS = importdata('DNS/t10m');
f10_DNS = importdata('DNS/t10f');
c10_DNS = importdata('DNS/t10c');

t20_VDR = importdata('V2F-DWX/t20_rad');
t20_VDN = importdata('V2F-DWX/t20');
t20_VN  = importdata('V2F-NO/t20');
t20_SN  = importdata('SA-NO/t20');
t20_DNS = importdata('DNS/t20m');
f20_DNS = importdata('DNS/t20f');
c20_DNS = importdata('DNS/t20c');

y = t01_SN(:,1);
yd = t01_DNS(1:2:end,1);

%% PLOTTING SECTION

fig=figure();
plot(yd,b_DNS(1:2:end,5),'ko','MarkerSize',4);
hold on
plot(y,b_SN(:,3),'k:');
plot(y,b_VN(:,3),'k-.');
plot(y,b_VDN(:,3),'k--');
set(gca,'FontSize',14);
xlabel('$$y$$','Interpreter','latex','FontSize',18);
ylabel('$$T$$','Interpreter','latex','FontSize',18);
text(1,0.9,'$$\tau=0.0$$','Interpreter','latex','FontSize',18,'HorizontalAlignment','center');
str1 = 'DNS';
str2 = 'SA';
str3 = 'V2F-NO';
str4 = 'V2F-DW';
legend({str1,str2,str3,str4},'Interpreter','latex','FontSize',12,'location','south');
legend boxoff
print(fig,strcat(path,'Tb.eps'),'-depsc');

fig=figure();
plot(yd,t01_DNS(1:2:end,5),'ko','MarkerSize',4);
hold on
plot(y,t01_SN(:,3),'k:');
plot(y,t01_VN(:,3),'k-.');
plot(y,t01_VDN(:,3),'k--');
plot(y,t01_VDR(:,3),'k-');
set(gca,'FontSize',14);
xlabel('$$y$$','Interpreter','latex','FontSize',18);
ylabel('$$T$$','Interpreter','latex','FontSize',18);
text(1,0.9,'$$\tau=0.1$$','Interpreter','latex','FontSize',18,'HorizontalAlignment','center');
str1 = 'DNS';
str2 = 'SA';
str3 = 'V2F-NO';
str4 = 'V2F-DW';
str5 = 'V2F-DWR';
legend({str1,str2,str3,str4,str5},'Interpreter','latex','FontSize',12,'location','south');
legend boxoff
print(fig,strcat(path,'Tt01.eps'),'-depsc');

fig=figure();
plot(yd,t1_DNS(1:2:end,5),'ko','MarkerSize',4);
hold on
plot(y,t1_SN(:,3),'k:');
plot(y,t1_VN(:,3),'k-.');
plot(y,t1_VDN(:,3),'k--');
plot(y,t1_VDR(:,3),'k-');
xlabel('$$y$$','Interpreter','latex','FontSize',18);
ylabel('$$T$$','Interpreter','latex','FontSize',18);
text(1,0.9,'$$\tau=1$$','Interpreter','latex','FontSize',18,'horizontalalignment','center');
str1 = 'DNS';
str2 = 'SA';
str3 = 'V2F-NO';
str4 = 'V2F-DW';
str5 = 'V2F-DWR';
legend({str1,str2,str3,str4,str5},'Interpreter','latex','FontSize',12,'location','south');
legend boxoff

print(fig,strcat(path,'Tt1.eps'),'-depsc');

fig=figure();
plot(yd,t5_DNS(1:2:end,5),'ko','MarkerSize',4);
hold on
plot(y,t5_SN(:,3),'k:');
plot(y,t5_VN(:,3),'k-.');
plot(y,t5_VDN(:,3),'k--');
plot(y,t5_VDR(:,3),'k-');
set(gca,'FontSize',14);
xlabel('$$y$$','Interpreter','latex','FontSize',18);
ylabel('$$T$$','Interpreter','latex','FontSize',18);
text(1,0.9,'$$\tau=5$$','Interpreter','latex','FontSize',18,'horizontalalignment','center');
str1 = 'DNS';
str2 = 'SA';
str3 = 'V2F-NO';
str4 = 'V2F-DW';
str5 = 'V2F-DWR';
legend({str1,str2,str3,str4,str5},'Interpreter','latex','FontSize',12,'location','south');
legend boxoff
print(fig,strcat(path,'Tt5.eps'),'-depsc');

fig=figure();
plot(yd,t10_DNS(1:2:end,5),'ko','MarkerSize',4);
hold on
plot(y,t10_SN(:,3),'k:');
plot(y,t10_VN(:,3),'k-.');
plot(y,t10_VDN(:,3),'k--');
plot(y,t10_VDR(:,3),'k-');
set(gca,'FontSize',14);
xlabel('$$y$$','Interpreter','latex','FontSize',18);
ylabel('$$T$$','Interpreter','latex','FontSize',18);
text(1,0.9,'$$\tau=10$$','Interpreter','latex','FontSize',18,'horizontalalignment','center');
str1 = 'DNS';
str2 = 'SA';
str3 = 'V2F-NO';
str4 = 'V2F-DW';
str5 = 'V2F-DWR';
legend({str1,str2,str3,str4,str5},'Interpreter','latex','FontSize',12,'location','south');
legend boxoff
print(fig,strcat(path,'Tt10.eps'),'-depsc');

fig=figure();
plot(yd,t20_DNS(1:2:end,5),'ko','MarkerSize',4);
hold on
plot(y,t20_SN(:,3),'k:');
plot(y,t20_VN(:,3),'k-.');
plot(y,t20_VDN(:,3),'k--');
plot(y,t20_VDR(:,3),'k-');
set(gca,'FontSize',14);
xlabel('$$y$$','Interpreter','latex','FontSize',18);
ylabel('$$T$$','Interpreter','latex','FontSize',18);
text(1,0.9,'$$\tau=20$$','Interpreter','latex','FontSize',18,'horizontalalignment','center');
str1 = 'DNS';
str2 = 'SA';
str3 = 'V2F-NO';
str4 = 'V2F-DW';
str5 = 'V2F-DWR';
legend({str1,str2,str3,str4,str5},'Interpreter','latex','FontSize',12,'location','south');
legend boxoff
print(fig,strcat(path,'Tt20.eps'),'-depsc');

fig=figure();
plot(yd,f_DNS(1:2:end,2),'ko','MarkerSize',4);
hold on
plot(y,b_SN(:,9),'k:');
plot(y,b_VN(:,9),'k-.');
plot(y,b_VDN(:,9),'k--');
set(gca,'FontSize',14);
xlabel('$$y$$','Interpreter','latex','FontSize',18);
ylabel('$$\overline{v^\prime T^\prime}$$','Interpreter','latex','FontSize',18);
text(1,f01_DNS(84,2)/2,'$$\tau=0.0$$','Interpreter','latex','FontSize',18,'horizontalalignment','center');
str1 = 'DNS';
str2 = 'SA';
str3 = 'V2F-NO';
str4 = 'V2F-DW';
legend({str1,str2,str3,str4},'Interpreter','latex','FontSize',12,'location','northwest');
legend boxoff
print(fig,strcat(path,'uTb.eps'),'-depsc');

fig=figure();
plot(yd,f01_DNS(1:2:end,2),'ko','MarkerSize',4);
hold on
plot(y,t01_SN(:,9),'k:');
plot(y,t01_VN(:,9),'k-.');
plot(y,t01_VDN(:,9),'k--');
plot(y,t01_VDR(:,9),'k-');
set(gca,'FontSize',14);
xlabel('$$y$$','Interpreter','latex','FontSize',18);
ylabel('$$\overline{v^\prime T^\prime}$$','Interpreter','latex','FontSize',18);
text(1,f01_DNS(84,2)/2,'$$\tau=0.1$$','Interpreter','latex','FontSize',18,'horizontalalignment','center');
str1 = 'DNS';
str2 = 'SA';
str3 = 'V2F-NO';
str4 = 'V2F-DW';
str5 = 'V2F-DWR';
legend({str1,str2,str3,str4,str5},'Interpreter','latex','FontSize',12,'location','northwest');
legend boxoff
print(fig,strcat(path,'uTt01.eps'),'-depsc');

fig=figure();
plot(yd,f1_DNS(1:2:end,2),'ko','MarkerSize',4);
hold on
plot(y,t1_SN(:,9),'k:');
plot(y,t1_VN(:,9),'k-.');
plot(y,t1_VDN(:,9),'k--');
plot(y,t1_VDR(:,9),'k-');
set(gca,'FontSize',14);
xlabel('$$y$$','Interpreter','latex','FontSize',18);
ylabel('$$\overline{v^\prime T^\prime}$$','Interpreter','latex','FontSize',18);
text(1,f1_DNS(84,2)/2,'$$\tau=1$$','Interpreter','latex','FontSize',18,'horizontalalignment','center');
str1 = 'DNS';
str2 = 'SA';
str3 = 'V2F-NO';
str4 = 'V2F-DW';
str5 = 'V2F-DWR';
legend({str1,str2,str3,str4,str5},'Interpreter','latex','FontSize',12,'location','northwest');
legend boxoff
print(fig,strcat(path,'uTt1.eps'),'-depsc');

fig=figure();
plot(yd,f5_DNS(1:2:end,2),'ko','MarkerSize',4);
hold on
plot(y,t5_SN(:,9),'k:');
plot(y,t5_VN(:,9),'k-.');
plot(y,t5_VDN(:,9),'k--');
plot(y,t5_VDR(:,9),'k-');
set(gca,'FontSize',14);
xlabel('$$y$$','Interpreter','latex','FontSize',18);
ylabel('$$\overline{v^\prime T^\prime}$$','Interpreter','latex','FontSize',18);
text(1,f5_DNS(84,2)/2,'$$\tau=5$$','Interpreter','latex','FontSize',18,'horizontalalignment','center');
str1 = 'DNS';
str2 = 'SA';
str3 = 'V2F-NO';
str4 = 'V2F-DW';
str5 = 'V2F-DWR';
legend({str1,str2,str3,str4,str5},'Interpreter','latex','FontSize',12,'location','northwest');
legend boxoff
print(fig,strcat(path,'uTt5.eps'),'-depsc');

fig=figure();
plot(yd,f10_DNS(1:2:end,2),'ko','MarkerSize',4);
hold on
plot(y,t10_SN(:,9),'k:');
plot(y,t10_VN(:,9),'k-.');
plot(y,t10_VDN(:,9),'k--');
plot(y,t10_VDR(:,9),'k-');
set(gca,'FontSize',14);
xlabel('$$y$$','Interpreter','latex','FontSize',18);
ylabel('$$\overline{v^\prime T^\prime}$$','Interpreter','latex','FontSize',18);
text(1,f10_DNS(84,2)/2,'$$\tau=10$$','Interpreter','latex','FontSize',18,'horizontalalignment','center');
str1 = 'DNS';
str2 = 'SA';
str3 = 'V2F-NO';
str4 = 'V2F-DW';
str5 = 'V2F-DWR';
legend({str1,str2,str3,str4,str5},'Interpreter','latex','FontSize',12,'location','northwest');
legend boxoff
print(fig,strcat(path,'uTt10.eps'),'-depsc');

fig=figure();
plot(yd,f20_DNS(1:2:end,2),'ko','MarkerSize',4);
hold on
plot(y,t20_SN(:,9),'k:');
plot(y,t20_VN(:,9),'k-.');
plot(y,t20_VDN(:,9),'k--');
plot(y,t20_VDR(:,9),'k-');
set(gca,'FontSize',14);
xlabel('$$y$$','Interpreter','latex','FontSize',18);
ylabel('$$\overline{v^\prime T^\prime}$$','Interpreter','latex','FontSize',18);
text(1,f20_DNS(84,2)/2,'$$\tau=20$$','Interpreter','latex','FontSize',18,'horizontalalignment','center');
str1 = 'DNS';
str2 = 'SA';
str3 = 'V2F-NO';
str4 = 'V2F-DW';
str5 = 'V2F-DWR';
legend({str1,str2,str3,str4,str5},'Interpreter','latex','FontSize',12,'location','northwest');
legend boxoff
print(fig,strcat(path,'uTt20.eps'),'-depsc');



%% PLOTTING ONE FIGURE Temperature

fig=figure('Position',[0 0 1200 800]);
axp = axes('Position',[0 0 1 1],'Visible','off');
[ha1,pos1] = tight_subplot(2,3,[0.03 0.015],[0.15 0.1],[0.07 0.07]);

axes(ha1(1))
pos1 = get(gca,'Position');
plot(yd,b_DNS(1:2:end,5),'ko','MarkerSize',4);
hold on
plot(y,b_SN(:,3),'k:');
plot(y,b_VN(:,3),'k-.');
plot(y,b_VDN(:,3),'k--');
set(gca,'FontSize',14);
xlabel('');
set(gca,'XTickLabel','');
ylabel('$$T$$','Interpreter','latex','FontSize',18);
text(1,0.9,'$$\tau=0.0$$','Interpreter','latex','FontSize',18,'HorizontalAlignment','center');
str1 = 'DNS';
str2 = 'SA';
str3 = 'V2F-NO';
str4 = 'V2F-DW';
legend({str1,str2,str3,str4},'Interpreter','latex','FontSize',12,'location','south');
legend boxoff

axes(ha1(2))
pos1 = get(gca,'Position');
plot(yd,t01_DNS(1:2:end,5),'ko','MarkerSize',4);
hold on
plot(y,t01_SN(:,3),'k:');
plot(y,t01_VN(:,3),'k-.');
plot(y,t01_VDN(:,3),'k--');
plot(y,t01_VDR(:,3),'k-');
set(gca,'FontSize',14);
xlabel('');
set(gca,'XTickLabel','')
ylabel('');
set(gca,'YTickLabel','')
text(1,0.9,'$$\tau=0.1$$','Interpreter','latex','FontSize',18,'HorizontalAlignment','center');
str1 = 'DNS';
str2 = 'SA';
str3 = 'V2F-NO';
str4 = 'V2F-DW';
str5 = 'V2F-DWR';
legend({str1,str2,str3,str4,str5},'Interpreter','latex','FontSize',12,'location','south');
legend boxoff


axes(ha1(3))
pos1 = get(gca,'Position');
plot(yd,t1_DNS(1:2:end,5),'ko','MarkerSize',4);
hold on
plot(y,t1_SN(:,3),'k:');
plot(y,t1_VN(:,3),'k-.');
plot(y,t1_VDN(:,3),'k--');
plot(y,t1_VDR(:,3),'k-');
set(gca,'FontSize',14);
xlabel('');
ylabel('');
set(gca,'XTickLabel','')
set(gca,'YTickLabel','')
text(1,0.9,'$$\tau=1$$','Interpreter','latex','FontSize',18,'HorizontalAlignment','center');
str1 = 'DNS';
str2 = 'SA';
str3 = 'V2F-NO';
str4 = 'V2F-DW';
str5 = 'V2F-DWR';
legend({str1,str2,str3,str4,str5},'Interpreter','latex','FontSize',12,'location','south');
legend boxoff

axes(ha1(4))
pos1 = get(gca,'Position');
plot(yd,t5_DNS(1:2:end,5),'ko','MarkerSize',4);
hold on
plot(y,t5_SN(:,3),'k:');
plot(y,t5_VN(:,3),'k-.');
plot(y,t5_VDN(:,3),'k--');
plot(y,t5_VDR(:,3),'k-');
set(gca,'FontSize',14);
xlabel('$$y$$','Interpreter','latex','FontSize',18);
ylabel('$$T$$','Interpreter','latex','FontSize',18);
text(1,0.9,'$$\tau=5$$','Interpreter','latex','FontSize',18,'HorizontalAlignment','center');
str1 = 'DNS';
str2 = 'SA';
str3 = 'V2F-NO';
str4 = 'V2F-DW';
str5 = 'V2F-DWR';
legend({str1,str2,str3,str4,str5},'Interpreter','latex','FontSize',12,'location','south');
legend boxoff

axes(ha1(5))
pos1 = get(gca,'Position');
plot(yd,t10_DNS(1:2:end,5),'ko','MarkerSize',4);
hold on
plot(y,t10_SN(:,3),'k:');
plot(y,t10_VN(:,3),'k-.');
plot(y,t10_VDN(:,3),'k--');
plot(y,t10_VDR(:,3),'k-');
set(gca,'FontSize',14);
xlabel('$$y$$','Interpreter','latex','FontSize',18);
ylabel('');
set(gca,'YTickLabel','')
text(1,0.9,'$$\tau=10$$','Interpreter','latex','FontSize',18,'HorizontalAlignment','center');
str1 = 'DNS';
str2 = 'SA';
str3 = 'V2F-NO';
str4 = 'V2F-DW';
str5 = 'V2F-DWR';
legend({str1,str2,str3,str4,str5},'Interpreter','latex','FontSize',12,'location','south');
legend boxoff

axes(ha1(6))
pos1 = get(gca,'Position');
plot(yd,t20_DNS(1:2:end,5),'ko','MarkerSize',4);
hold on
plot(y,t20_SN(:,3),'k:');
plot(y,t20_VN(:,3),'k-.');
plot(y,t20_VDN(:,3),'k--');
plot(y,t20_VDR(:,3),'k-');
set(gca,'FontSize',14);
xlabel('$$y$$','Interpreter','latex','FontSize',18);
ylabel('');
set(gca,'YTickLabel','')
text(1,0.9,'$$\tau=20$$','Interpreter','latex','FontSize',18,'HorizontalAlignment','center');
str1 = 'DNS';
str2 = 'SA';
str3 = 'V2F-NO';
str4 = 'V2F-DW';
str5 = 'V2F-DWR';
legend({str1,str2,str3,str4,str5},'Interpreter','latex','FontSize',12,'location','south');
legend boxoff


print(fig,strcat(path,'Tempgrey.eps'),'-depsc')

%% PLOTTING ONE FIGURE

fig=figure('Position',[0 0 1200 800]);
axp = axes('Position',[0 0 1 1],'Visible','off');
[ha1,pos1] = tight_subplot(2,3,[0.03 0.015],[0.15 0.1],[0.07 0.07]);

axes(ha1(1))
pos1 = get(gca,'Position');
plot(yd,f_DNS(1:2:end,2),'ko','MarkerSize',4);
hold on
plot(y,b_SN(:,9),'k:');
plot(y,b_VN(:,9),'k-.');
plot(y,b_VDN(:,9),'k--');
set(gca,'FontSize',14);
xlabel('');
set(gca,'XTickLabel','');
ylabel('$$\overline{v^\prime \theta^\prime}$$','Interpreter','latex','FontSize',18);
text(1,0.9,'$$\tau=0.0$$','Interpreter','latex','FontSize',18,'HorizontalAlignment','center');
str1 = 'DNS';
str2 = 'SA';
str3 = 'V2F-NO';
str4 = 'V2F-DW';
legend({str1,str2,str3,str4},'Interpreter','latex','FontSize',12,'location','south');
legend boxoff

axes(ha1(2))
pos1 = get(gca,'Position');
plot(yd,f01_DNS(1:2:end,2),'ko','MarkerSize',4);
hold on
plot(y,t01_SN(:,9),'k:');
plot(y,t01_VN(:,9),'k-.');
plot(y,t01_VDN(:,9),'k--');
plot(y,t01_VDR(:,9),'k-');
set(gca,'FontSize',14);
xlabel('');
set(gca,'XTickLabel','')
ylabel('');
set(gca,'YTickLabel','')
text(1,0.9,'$$\tau=0.1$$','Interpreter','latex','FontSize',18,'HorizontalAlignment','center');
str1 = 'DNS';
str2 = 'SA';
str3 = 'V2F-NO';
str4 = 'V2F-DW';
str5 = 'V2F-DWR';
legend({str1,str2,str3,str4,str5},'Interpreter','latex','FontSize',12,'location','south');
legend boxoff


axes(ha1(3))
pos1 = get(gca,'Position');
plot(yd,f1_DNS(1:2:end,2),'ko','MarkerSize',4);
hold on
plot(y,t1_SN(:,9),'k:');
plot(y,t1_VN(:,9),'k-.');
plot(y,t1_VDN(:,9),'k--');
plot(y,t1_VDR(:,9),'k-');
set(gca,'FontSize',14);
xlabel('');
ylabel('');
set(gca,'XTickLabel','')
set(gca,'YTickLabel','')
text(1,0.9,'$$\tau=1$$','Interpreter','latex','FontSize',18,'HorizontalAlignment','center');
str1 = 'DNS';
str2 = 'SA';
str3 = 'V2F-NO';
str4 = 'V2F-DW';
str5 = 'V2F-DWR';
legend({str1,str2,str3,str4,str5},'Interpreter','latex','FontSize',12,'location','south');
legend boxoff

axes(ha1(4))
pos1 = get(gca,'Position');
plot(yd,f5_DNS(1:2:end,2),'ko','MarkerSize',4);
hold on
plot(y,t5_SN(:,9),'k:');
plot(y,t5_VN(:,9),'k-.');
plot(y,t5_VDN(:,9),'k--');
plot(y,t5_VDR(:,9),'k-');
set(gca,'FontSize',14);
xlabel('$$y$$','Interpreter','latex','FontSize',18);
ylabel('$$\overline{v^\prime \theta^\prime}$$','Interpreter','latex','FontSize',18);
text(1,0.9,'$$\tau=5$$','Interpreter','latex','FontSize',18,'HorizontalAlignment','center');
str1 = 'DNS';
str2 = 'SA';
str3 = 'V2F-NO';
str4 = 'V2F-DW';
str5 = 'V2F-DWR';
legend({str1,str2,str3,str4,str5},'Interpreter','latex','FontSize',12,'location','south');
legend boxoff

axes(ha1(5))
pos1 = get(gca,'Position');
plot(yd,f10_DNS(1:2:end,2),'ko','MarkerSize',4);
hold on
plot(y,t10_SN(:,9),'k:');
plot(y,t10_VN(:,9),'k-.');
plot(y,t10_VDN(:,9),'k--');
plot(y,t10_VDR(:,9),'k-');
set(gca,'FontSize',14);
xlabel('$$y$$','Interpreter','latex','FontSize',18);
ylabel('');
set(gca,'YTickLabel','')
text(1,0.9,'$$\tau=10$$','Interpreter','latex','FontSize',18,'HorizontalAlignment','center');
str1 = 'DNS';
str2 = 'SA';
str3 = 'V2F-NO';
str4 = 'V2F-DW';
str5 = 'V2F-DWR';
legend({str1,str2,str3,str4,str5},'Interpreter','latex','FontSize',12,'location','south');
legend boxoff

axes(ha1(6))
pos1 = get(gca,'Position');
plot(yd,f20_DNS(1:2:end,2),'ko','MarkerSize',4);
hold on
plot(y,t20_SN(:,9),'k:');
plot(y,t20_VN(:,9),'k-.');
plot(y,t20_VDN(:,9),'k--');
plot(y,t20_VDR(:,9),'k-');
set(gca,'FontSize',14);
xlabel('$$y$$','Interpreter','latex','FontSize',18);
ylabel('');
set(gca,'YTickLabel','')
text(1,0.9,'$$\tau=20$$','Interpreter','latex','FontSize',18,'HorizontalAlignment','center');
str1 = 'DNS';
str2 = 'SA';
str3 = 'V2F-NO';
str4 = 'V2F-DW';
str5 = 'V2F-DWR';
legend({str1,str2,str3,str4,str5},'Interpreter','latex','FontSize',12,'location','south');
legend boxoff


print(fig,strcat(path,'fluxgrey.eps'),'-depsc')

