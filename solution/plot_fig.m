
%% Plotting figures from solutions

close all
clear all

%% IMPORT DATA

t01_VDR = importdata('V2F-DWX/t01_rad');
t01_VDN = importdata('V2F-DWX/t01');
t01_VN  = importdata('V2F-NO/t01');
t01_SN  = importdata('SA-NO/t01');
t01_DNS = importdata('DNS/t01m');
f01_DNS = importdata('DNS/t01f');

t1_VDR = importdata('V2F-DWX/t1_rad');
t1_VDN = importdata('V2F-DWX/t1');
t1_VN  = importdata('V2F-NO/t1');
t1_SN  = importdata('SA-NO/t1');
t1_DNS = importdata('DNS/t1m');
f1_DNS = importdata('DNS/t1f');

t5_VDR = importdata('V2F-DWX/t5_rad');
t5_VDN = importdata('V2F-DWX/t5');
t5_VN  = importdata('V2F-NO/t5');
t5_SN  = importdata('SA-NO/t5');
t5_DNS = importdata('DNS/t5m');
f5_DNS = importdata('DNS/t5f');

t10_VDR = importdata('V2F-DWX/t10_rad');
t10_VDN = importdata('V2F-DWX/t10');
t10_VN  = importdata('V2F-NO/t10');
t10_SN  = importdata('SA-NO/t10');
t10_DNS = importdata('DNS/t10m');
f10_DNS = importdata('DNS/t10f');

t20_VDR = importdata('V2F-DWX/t20_rad');
t20_VDN = importdata('V2F-DWX/t20');
t20_VN  = importdata('V2F-NO/t20');
t20_SN  = importdata('SA-NO/t20');
t20_DNS = importdata('DNS/t20m');
f20_DNS = importdata('DNS/t20f');

y = t01_SN(:,1);
yd = t01_DNS(1:2:end,1);

%% PLOTTING SECTION

fig=figure();
plot(yd,t01_DNS(1:2:end,5),'ko','MarkerSize',4);
hold on
plot(y,t01_SN(:,3),'k:');
plot(y,t01_VN(:,3),'k-.');
plot(y,t01_VDN(:,3),'k--');
plot(y,t01_VDR(:,3),'k-');
set(gca,'FontSize',14);
xlabel('$$y$$','Interpreter','latex');
ylabel('$$T$$','Interpreter','latex');
text(1,0.9,'$$\tau=0.1$$','Interpreter','latex','FontSize',14);
print(fig,'Figures/Tt01.eps','-depsc');

fig=figure();
plot(yd,t1_DNS(1:2:end,5),'ko','MarkerSize',4);
hold on
plot(y,t1_SN(:,3),'k:');
plot(y,t1_VN(:,3),'k-.');
plot(y,t1_VDN(:,3),'k--');
plot(y,t1_VDR(:,3),'k-');
xlabel('$$y$$','Interpreter','latex');
ylabel('$$T$$','Interpreter','latex');
text(1,0.9,'$$\tau=1$$','Interpreter','latex','FontSize',14);
print(fig,'Figures/Tt1.eps','-depsc');

fig=figure();
plot(yd,t5_DNS(1:2:end,5),'ko','MarkerSize',4);
hold on
plot(y,t5_SN(:,3),'k:');
plot(y,t5_VN(:,3),'k-.');
plot(y,t5_VDN(:,3),'k--');
plot(y,t5_VDR(:,3),'k-');
set(gca,'FontSize',14);
xlabel('$$y$$','Interpreter','latex');
ylabel('$$T$$','Interpreter','latex');
text(1,0.9,'$$\tau=5$$','Interpreter','latex','FontSize',14);
print(fig,'Figures/Tt5.eps','-depsc');

fig=figure();
plot(yd,t10_DNS(1:2:end,5),'ko','MarkerSize',4);
hold on
plot(y,t10_SN(:,3),'k:');
plot(y,t10_VN(:,3),'k-.');
plot(y,t10_VDN(:,3),'k--');
plot(y,t10_VDR(:,3),'k-');
set(gca,'FontSize',14);
xlabel('$$y$$','Interpreter','latex');
ylabel('$$T$$','Interpreter','latex');
text(1,0.9,'$$\tau=10$$','Interpreter','latex','FontSize',14);
print(fig,'Figures/Tt10.eps','-depsc');

fig=figure();
plot(yd,t20_DNS(1:2:end,5),'ko','MarkerSize',4);
hold on
plot(y,t20_SN(:,3),'k:');
plot(y,t20_VN(:,3),'k-.');
plot(y,t20_VDN(:,3),'k--');
plot(y,t20_VDR(:,3),'k-');
set(gca,'FontSize',14);
xlabel('$$y$$','Interpreter','latex');
ylabel('$$T$$','Interpreter','latex');
text(1,0.9,'$$\tau=20$$','Interpreter','latex','FontSize',14);
print(fig,'Figures/Tt20.eps','-depsc');

fig=figure();
plot(yd,f01_DNS(1:2:end,2),'ko','MarkerSize',4);
hold on
plot(y,t01_SN(:,9),'k:');
plot(y,t01_VN(:,9),'k-.');
plot(y,t01_VDN(:,9),'k--');
plot(y,t01_VDR(:,9),'k-');
set(gca,'FontSize',14);
xlabel('$$y$$','Interpreter','latex');
ylabel('$$\overline{v^\prime T^\prime}$$','Interpreter','latex');
text(1,f01_DNS(84,2)/2,'$$\tau=0.1$$','Interpreter','latex','FontSize',14);
print(fig,'Figures/uTt01.eps','-depsc');

fig=figure();
plot(yd,f1_DNS(1:2:end,2),'ko','MarkerSize',4);
hold on
plot(y,t1_SN(:,9),'k:');
plot(y,t1_VN(:,9),'k-.');
plot(y,t1_VDN(:,9),'k--');
plot(y,t1_VDR(:,9),'k-');
set(gca,'FontSize',14);
xlabel('$$y$$','Interpreter','latex');
ylabel('$$\overline{v^\prime T^\prime}$$','Interpreter','latex');
text(1,f1_DNS(84,2)/2,'$$\tau=1$$','Interpreter','latex','FontSize',14);
print(fig,'Figures/uTt1.eps','-depsc');

fig=figure();
plot(yd,f5_DNS(1:2:end,2),'ko','MarkerSize',4);
hold on
plot(y,t5_SN(:,9),'k:');
plot(y,t5_VN(:,9),'k-.');
plot(y,t5_VDN(:,9),'k--');
plot(y,t5_VDR(:,9),'k-');
set(gca,'FontSize',14);
xlabel('$$y$$','Interpreter','latex');
ylabel('$$\overline{v^\prime T^\prime}$$','Interpreter','latex');
text(1,f5_DNS(84,2)/2,'$$\tau=5$$','Interpreter','latex','FontSize',14);
print(fig,'Figures/uTt5.eps','-depsc');

fig=figure();
plot(yd,f10_DNS(1:2:end,2),'ko','MarkerSize',4);
hold on
plot(y,t10_SN(:,9),'k:');
plot(y,t10_VN(:,9),'k-.');
plot(y,t10_VDN(:,9),'k--');
plot(y,t10_VDR(:,9),'k-');
set(gca,'FontSize',14);
xlabel('$$y$$','Interpreter','latex');
ylabel('$$\overline{v^\prime T^\prime}$$','Interpreter','latex');
text(1,f10_DNS(84,2)/2,'$$\tau=10$$','Interpreter','latex','FontSize',14);
print(fig,'Figures/uTt10.eps','-depsc');

fig=figure();
plot(yd,f20_DNS(1:2:end,2),'ko','MarkerSize',4);
hold on
plot(y,t20_SN(:,9),'k:');
plot(y,t20_VN(:,9),'k-.');
plot(y,t20_VDN(:,9),'k--');
plot(y,t20_VDR(:,9),'k-');
set(gca,'FontSize',14);
xlabel('$$y$$','Interpreter','latex');
ylabel('$$\overline{v^\prime T^\prime}$$','Interpreter','latex');
text(1,f20_DNS(84,2)/2,'$$\tau=20$$','Interpreter','latex','FontSize',14);
print(fig,'Figures/uTt20.eps','-depsc');













