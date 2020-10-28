%% Test the velocity transformation
close all; 
% clear all;

%--------------------------------------------------------------------------
%       Include folders
%--------------------------------------------------------------------------

addpath('mesh');                % functions for the mesh
addpath('function');            % functions for the solver
addpath('turbmodels');          % functions of the turbulence models
addpath('fluxmodels');          % functions of the turbulent flux models
addpath('radiation');           % functions for the radiative calculations

% -----  channel height  -----
height = 2;

% -----  number of mesh points  -----
n = 400;

% -----  discretization  -----
% discr = 'finitediff' ... finite difference discretization; 
%                          requires additional parameters: 
%                          fact ... stretching factor for hyperbolic tan 
%                          ns   ... number of stencil points to the left 
%                                   and right of central differencing
%                                   scheme. So, ns = 1 is second order.
% discr = 'chebyshev' ...  Chebyshev discretization
%                          Note, this discretization is very unstable.
%                          Therefore it is best to first start with second
%                          order finite difference scheme and then switch
discr = 'finitediff';

% -----  streching factor and stencil for finite difference discretization
fact = 6;
ns = 1;

%% ------------------------------------------------------------------------
%
%      Generate mesh and angular mesh
%
[MESH] = mesh(height, n, fact, ns, discr);

%% START IMPORTING
% 
ctb = importdata('solution/DNS/bm');
cth = importdata('solution/DNS/H2Om');
ct01= importdata('solution/DNS/t01rm');
ct1 = importdata('solution/DNS/t1rm');
ct10= importdata('solution/DNS/t10rm');
ctRb= importdata('solution/DNS/Rbm');
cthr= importdata('solution/DNS/H2ORm');



restb  = importdata('~/Dropbox/Allresults/JFM/bench/resid');
rest01 = importdata('~/Dropbox/Allresults/JQSRT/tau01/resid');
rest1  = importdata('~/Dropbox/Allresults/JQSRT/tau1/resid');
rest10 = importdata('~/Dropbox/Allresults/JQSRT/tau10/resid');
resth  = importdata('~/Dropbox/Allresults/JQSRT/H2O/resid');
restH  = importdata('~/Dropbox/Allresults/new_higRe/H2O/resid');
restR  = importdata('~/Dropbox/Allresults/new_higRe/bench/resid');

set(groot,'DefaultLineLineWidth',0.8);

%%

for i=2:19
    for j=1:166
        cb(j,i) = ctb(166-j+1,i);
    end
end
for i=2:18
    for j=1:166
        resb(j,i) = restb(166-j+1,i);
    end
end
for i=2:18
    for j=1:190
        ch(j,i) = cth(190-j+1,i);
        resh(j,i) = resth(190-j+1,i);        
    end
end
for i=2:18
    for j=1:414
        cR(j,i) = ctRb(414-j+1,i);
        cH(j,i) = cthr(414-j+1,i);
        resH(j,i) = restH(414-j+1,i);        
        resR(j,i) = restR(414-j+1,i);        
    end
end
for i=2:18
    for j=1:166
        c01(j,i) = ct01(166-j+1,i);
        c1(j,i)  = ct1(166-j+1,i);
        res01(j,i) = rest01(166-j+1,i);
        res1(j,i) = rest1(166-j+1,i);
    end
    for j=1:190
        c10(j,i) = ct10(190-j+1,i);
        res10(j,i) = rest10(190-j+1,i);
    end
end
cb(:,1) = ctb(:,1);
ch(:,1) = cth(:,1);
c01(:,1) = ct01(:,1);
c1(:,1) = ct1(:,1);
c10(:,1) = ct10(:,1);
cH(:,1) = cthr(:,1);
cR(:,1) = ctRb(:,1);
% 
% cb  = ctb;
% ch  = cth;
% c01 = ct01;
% c1  = ct1;
% c10 = ct10;
% cH  = cthr;
% cR  = ctRb;

y = MESH.y;

ub = interp1(cb(:,1) ,cb(:,4) ,y,'pchip');
uh = interp1(ch(:,1) ,ch(:,4) ,y,'pchip');
u01= interp1(c01(:,1),c01(:,4),y,'pchip');
u1 = interp1(c1(:,1) ,c1(:,4) ,y,'pchip');
u10= interp1(c10(:,1),c10(:,4),y,'pchip');
uR = interp1(cR(:,1) ,cR(:,4) ,y,'pchip');
uH = interp1(cH(:,1) ,cH(:,4) ,y,'pchip');

rb = ones(n,1);
r01= interp1(c01(:,1),c01(:,14),y,'pchip');
r1 = interp1(c1(:,1) ,c1(:,14) ,y,'pchip');
r10= interp1(c10(:,1),c10(:,14),y,'pchip');
rh = interp1(ch(:,1) ,ch(:,14) ,y,'pchip');
rR = interp1(cR(:,1) ,cR(:,14) ,y,'pchip');
rH = interp1(cH(:,1) ,cH(:,14) ,y,'pchip');
mR = (rR).^(-1.15);
mH = (rH).^(-1.15);

tvb  = (MESH.ddy*ub) /2800;
tv01 = (MESH.ddy*u01)/3750;
tv1  = (MESH.ddy*u1) /3750;
tv10 = (MESH.ddy*u10)/3750;
tvh  = (MESH.ddy*uh) /3750;
tvR  = mR.*(MESH.ddy*uR) /16700;
tvH  = mH.*(MESH.ddy*uH) /16700;


utb  = sqrt(tvb(1)/rb(1));
uth  = sqrt(tvh(1)/rh(1));
ut01 = sqrt(tv01(1)/r01(1));
ut1  = sqrt(tv1(1)/r1(1));
ut10 = sqrt(tv10(1)/r10(1));
utR  = sqrt(tvR(1)/rR(1));
utH  = sqrt(tvH(1)/rH(1));

ubVD  = velTransVD(ub./utb,rb);
uhVD  = velTransVD(uh./uth,rh);
u01VD = velTransVD(u01./ut01,r01);
u1VD  = velTransVD(u1./ut1,r1);
u10VD = velTransVD(u10./ut10,r10);
uRVD = velTransVD(uR./utR,rR);
uHVD = velTransVD(uH./utH,rH);

Reb = 2900*utb;
Reh = 3750*uth;
Re01 = 3750*ut01;
Re1 = 3750*ut1;
Re10 = 3750*ut10;
ReR = 16700*utR;
ReH = 16700*utH;

utsb  = sqrt(tvb(1)./rb);
utsh  = sqrt(tvh(1)./rh);
uts01 = sqrt(tv01(1)./r01);
uts1  = sqrt(tv1(1)./r1);
uts10  = sqrt(tv10(1)./r10);
utsR  = sqrt(tvR(1)/rR);
utsH  = sqrt(tvH(1)/rH);

Resb = 2900*utb*sqrt(rb/rb(1));
Resh = 3750*uth*sqrt(rh/rh(1));
Res01 = 3750*ut01*sqrt(r01/r01(1));
Res1 = 3750*ut1*sqrt(r1/r1(1));
Res10 = 3750*ut10*sqrt(r10/r10(1));
ResR = 16700*utR*sqrt(rR/rR(1))./mR;
ResH = 16700*utH*sqrt(rH/rH(1))./mH;

usb  = velTransSLS(ubVD,Resb,MESH);
ush  = velTransSLS(uhVD,Resh,MESH);
us01 = velTransSLS(u01VD,Res01,MESH);
us1  = velTransSLS(u1VD,Res1,MESH);
us10 = velTransSLS(u10VD,Res10,MESH);
usR = velTransSLS(uRVD,ResR,MESH);
usH = velTransSLS(uHVD,ResH,MESH);

ysb = y.*Resb;
ysh = y.*Resh;
ys01= y.*Res01;
ys1 = y.*Res1;
ys10= y.*Res10;
ysR = y.*ResR;
ysH = y.*ResH;

ypb = y.*Reb;
yph = y.*Reh;
yp01= y.*Re01;
yp1 = y.*Re1;
yp10= y.*Re10;
ypR = y.*ReR;
ypH = y.*ReH;

figure(1)
plot(ysb(1:end/2),usb(1:end/2),'k-','LineWidth',1.1);
hold on
plot(ysh(1:end/2),ush(1:end/2),'ro','Markersize',1);
plot(ys01(1:end/2),us01(1:end/2),'bd','Markersize',1);
plot(ys1(1:end/2),us1(1:end/2),'gs','Markersize',1);
plot(ys10(1:end/2),us10(1:end/2),'k^','Markersize',1);
plot(ysR(1:end/2),usR(1:end/2),'--','Markersize',1);
plot(ysH(1:end/2),usH(1:end/2),'-.','Markersize',1);
set(gca,'xscale','log')

figure(2)
hold on
plot(yph(1:end/2),uhVD(1:end/2),'ro','Markersize',1);
plot(yp01(1:end/2),u01VD(1:end/2),'bd','Markersize',1);
plot(yp1(1:end/2),u1VD(1:end/2),'gs','Markersize',1);
plot(yp10(1:end/2),u10VD(1:end/2),'k^','Markersize',1);
plot(ypR(1:end/2),uRVD(1:end/2),'--','Markersize',1);
plot(ypH(1:end/2),uHVD(1:end/2),'-.','Markersize',1);
set(gca,'xscale','log')

figure(3)
hold on
plot(yph(1:end/2),uh(1:end/2)/uth);
plot(yp01(1:end/2),u01(1:end/2)/ut01);
plot(yp1(1:end/2),u1(1:end/2)/ut1);
plot(yp10(1:end/2),u10(1:end/2)/ut10);
plot(ypR(1:end/2),uR(1:end/2)/utR);
plot(ypH(1:end/2),uH(1:end/2)/utH);
set(gca,'xscale','log')

%DIAGNOSTIC FUNCTION
Dh = (MESH.ddy*uhVD)./Resh;
D01 = (MESH.ddy*u01VD)./Res01;
D1 = (MESH.ddy*u1VD)./Res1;
D10 = (MESH.ddy*u10VD)./Res10;
DR = (MESH.ddy*uRVD)./ResR;
DH = (MESH.ddy*uHVD)./ResH;

figure(4)
hold on
plot(ysh(1:end/2) ,Dh(1:end/2),'ro','Markersize',1);
plot(ys01(1:end/2),D01(1:end/2),'bd','Markersize',1);
plot(ys1(1:end/2) ,D1(1:end/2),'gs','Markersize',1);
plot(ysR(1:end/2) ,DR(1:end/2),'--','Markersize',1);
plot(ysH(1:end/2) ,DH(1:end/2),'-.','Markersize',1);
set(gca,'xscale','log','yscale','log')
box on
%set(gca,'yscale','log')

%% ANISOTROPIES

% clear ut01 ut1 ut10 Res01 Res1 Res10;
% 
% ut01 = zeros(length(r01),1);
% ut1  = zeros(length(r1),1);
% ut10 = zeros(length(r10),1);
% utR = zeros(length(rR),1);
% utH = zeros(length(rH),1);
% Res01 = zeros(length(r01),1);
% Res1 = zeros(length(r1),1);
% Res10 = zeros(length(r10),1);
% ResR = zeros(length(rR),1);
% ResH = zeros(length(rH),1);

% ut01(1:end/2) = sqrt(tv01(1)./r01(1:end/2));
% ut1(1:end/2)  = sqrt(tv1(1)./r1(1:end/2));
% ut10(1:end/2)  = sqrt(tv10(1)./r10(1:end/2));
% utR(1:end/2)  = sqrt(tvR(1)./rR(1:end/2));
% utH(1:end/2)  = sqrt(tvH(1)./rH(1:end/2));
% ut01(end/2+1:end) = sqrt(-tv01(end)./r01(end/2+1:end));
% ut1(end/2+1:end)  = sqrt(-tv1(end)./r1(end/2+1:end));
% ut10(end/2+1:end)  = sqrt(-tv10(end)./r10(end/2+1:end));
% utR(end/2+1:end)  = sqrt(-tvR(end)./rR(end/2+1:end));
% utH(end/2+1:end)  = sqrt(-tvH(end)./rH(end/2+1:end));
% Res01(1:end/2) = 3750*ut01(1:end/2).*r01(1:end/2); %.*sqrt(r01(1:end/2)/r01(1));
% Res1(1:end/2) = 3750*ut1(1:end/2).*r1(1:end/2); %.*sqrt(r1(1:end/2)/r1(1));
% Res10(1:end/2) = 3750*ut10(1:end/2).*r10(1:end/2); %.*sqrt(r10(1:end/2)/r10(1));
% ResR(1:end/2) = 16700*utR(1:end/2).*rR(1:end/2); %.*sqrt(r10(1:end/2)/r10(1));
% ResH(1:end/2) = 16700*utH(1:end/2).*rH(1:end/2); %.*sqrt(r10(1:end/2)/r10(1));
% Res01(end/2+1:end) = 3750*ut01(end/2+1:end).*r01(end/2+1:end); %.*sqrt(r01(end/2:end)/r01(end));
% Res1(end/2+1:end) = 3750*ut1(end/2+1:end).*r1(end/2+1:end); %.*sqrt(r1(end/2:end)/r1(end));
% Res10(end/2+1:end) = 3750*ut10(end/2+1:end).*r10(end/2+1:end); %.*sqrt(r10(end/2:end)/r10(end));
% ResR(end/2+1:end) = 16700*utR(end/2+1:end).*rR(end/2+1:end); %.*sqrt(r10(end/2:end)/r10(end));
% ResH(end/2+1:end) = 16700*utH(end/2+1:end).*rH(end/2+1:end); %.*sqrt(r10(end/2:end)/r10(end));
% 
% ysb = y.*Resb;
% ysh = y.*Resh;
% ys01= y.*Res01;
% ys1 = y.*Res1;
% ys10= y.*Res10;
% ysR = y.*ResR;
% ysH = y.*ResH;

figure(5)
hold on
plot(y(1:end/2),rb(1:end/2),'k-','LineWidth',1.1);
plot(y(1:end/2),rh(1:end/2),'ro','Markersize',4);
plot(y(1:end/2),r01(1:end/2),'bd','Markersize',4);
plot(y(1:end/2),r1(1:end/2),'gs','Markersize',4);
plot(y(1:end/2),r10(1:end/2),'k^','Markersize',4);
plot(y(1:end/2),rR(1:end/2),'--','Markersize',4);
plot(y(1:end/2),rH(1:end/2),'-.','Markersize',4);


%%

wnb  = interp1(cb(:,1),resb(:,2).^2./(resb(:,2).^2+resb(:,3).^2+resb(:,4).^2),y,'pchip');
wn01 = interp1(c01(:,1),res01(:,2).^2./(res01(:,2).^2+res01(:,3).^2+res01(:,4).^2),y,'pchip');
wn1  = interp1(c1(:,1) ,res1(:,2).^2./(res1(:,2).^2+res1(:,3).^2+res1(:,4).^2) ,y,'pchip');
wn10 = interp1(c10(:,1),res10(:,2).^2./(res10(:,2).^2+res10(:,3).^2+res10(:,4).^2),y,'pchip');
wnh  = interp1(ch(:,1) ,resh(:,2).^2./(resh(:,2).^2+resh(:,3).^2+resh(:,4).^2) ,y,'pchip');
wnH  = interp1(cH(:,1),resH(:,2).^2./(resH(:,2).^2+resH(:,3).^2+resH(:,4).^2),y,'pchip');
wnR  = interp1(cR(:,1) ,resR(:,2).^2./(resR(:,2).^2+resR(:,3).^2+resR(:,4).^2) ,y,'pchip');

spb  = interp1(cb(:,1),resb(:,3).^2./(resb(:,2).^2+resb(:,3).^2+resb(:,4).^2),y,'pchip');
sp01= interp1(c01(:,1),res01(:,3).^2./(res01(:,2).^2+res01(:,3).^2+res01(:,4).^2),y,'pchip');
sp1 = interp1(c1(:,1) ,res1(:,3).^2./(res1(:,2).^2+res1(:,3).^2+res1(:,4).^2) ,y,'pchip');
sp10= interp1(c10(:,1),res10(:,3).^2./(res10(:,2).^2+res10(:,3).^2+res10(:,4).^2),y,'pchip');
sph  = interp1(ch(:,1) ,resh(:,3).^2./(resh(:,2).^2+resh(:,3).^2+resh(:,4).^2) ,y,'pchip');
spH  = interp1(cH(:,1),resH(:,3).^2./(resH(:,2).^2+resH(:,3).^2+resH(:,4).^2),y,'pchip');
spR  = interp1(cR(:,1) ,resR(:,3).^2./(resR(:,2).^2+resR(:,3).^2+resR(:,4).^2) ,y,'pchip');

stb  = interp1(cb(:,1),resb(:,4).^2./(resb(:,2).^2+resb(:,3).^2+resb(:,4).^2),y,'pchip');
st01= interp1(c01(:,1),res01(:,4).^2./(res01(:,2).^2+res01(:,3).^2+res01(:,4).^2),y,'pchip');
st1 = interp1(c1(:,1) ,res1(:,4).^2./(res1(:,2).^2+res1(:,3).^2+res1(:,4).^2) ,y,'pchip');
st10= interp1(c10(:,1),res10(:,4).^2./(res10(:,2).^2+res10(:,3).^2+res10(:,4).^2),y,'pchip');
sth  = interp1(ch(:,1) ,resh(:,4).^2./(resh(:,2).^2+resh(:,3).^2+resh(:,4).^2) ,y,'pchip');
stH  = interp1(cH(:,1),resH(:,4).^2./(resH(:,2).^2+resH(:,3).^2+resH(:,4).^2),y,'pchip');
stR  = interp1(cR(:,1) ,resR(:,4).^2./(resR(:,2).^2+resR(:,3).^2+resR(:,4).^2) ,y,'pchip');


figure(9)
plot(ysb(1:end/2),wnb(1:end/2)-1./3,'k-','Markersize',1);
hold on
plot(ysh(1:end/2),wnh(1:end/2)-1./3,'ro','Markersize',1);
plot(ys01(1:end/2),wn01(1:end/2)-1./3,'bd','Markersize',1);
plot(ys1(1:end/2),wn1(1:end/2)-1./3,'gs','Markersize',1);
plot(ys10(1:end/2),wn10(1:end/2)-1./3,'k^','Markersize',1);
plot(ysR(1:end/2),wnR(1:end/2)-1./3,'--','Markersize',1);
plot(ysH(1:end/2),wnH(1:end/2)-1./3,'-.','Markersize',1);
plot(ysb(1:end/2),spb(1:end/2)-1./3,'k-','Markersize',1);
plot(ysh(1:end/2),sph(1:end/2)-1./3,'ro','Markersize',1);
plot(ys01(1:end/2),sp01(1:end/2)-1./3,'bd','Markersize',1);
plot(ys1(1:end/2),sp1(1:end/2)-1./3,'gs','Markersize',1);
plot(ys10(1:end/2),sp10(1:end/2)-1./3,'k^','Markersize',1);
plot(ysR(1:end/2),spR(1:end/2)-1./3,'--','Markersize',1);
plot(ysH(1:end/2),spH(1:end/2)-1./3,'-.','Markersize',1);
plot(ysb(1:end/2),stb(1:end/2)-1./3,'k-','Markersize',1);
plot(ysh(1:end/2),sth(1:end/2)-1./3,'ro','Markersize',1);
plot(ys01(1:end/2),st01(1:end/2)-1./3,'bd','Markersize',1);
plot(ys1(1:end/2),st1(1:end/2)-1./3,'gs','Markersize',1);
plot(ys10(1:end/2),st10(1:end/2)-1./3,'k^','Markersize',1);
plot(ysR(1:end/2),stR(1:end/2)-1./3,'--','Markersize',1);
plot(ysH(1:end/2),stH(1:end/2)-1./3,'-.','Markersize',1);
set(gca,'xscale','log')
xlim([0.1 1e2])

%%

strb  = derivativeX2(usb,ysb);
str01 = derivativeX2(us01,ys01);
str1 = derivativeX2(us1,ys1);
str10 = derivativeX2(us10,ys10);
strh = derivativeX2(ush,ysh);
strH = derivativeX2(usH,ysH);
strR = derivativeX2(usR,ysR);

strb = strb./interp1(y,strb,0.0,'pchip');
str01 = str01./interp1(y,str01,0.0,'pchip');
str1 = str1./interp1(y,str1,0.0,'pchip');
str10 = str10./interp1(y,str10,0.0,'pchip');
strh = strh./interp1(y,strh,0.0,'pchip');
strH = strH./interp1(y,strH,0.0,'pchip');
strR = strR./interp1(y,strR,0.0,'pchip');

dst01 = str01-strb;
dst1  = str1-strb;
dst10 = str10-strb;
dsth  = strh-strb;
dstH  = strH-strb;
dstR  = strR-strb;
% 
% dst01 = integralX4(dst01,y);
% dst1  = integralX4(dst1,y);
% dst10 = integralX4(dst10,y);
% dsth  = integralX4(dsth,y);
% dstH  = integralX4(dstH,y);
% dstR  = integralX4(dstR,y);


s01 = min(dst01(1:end/2));
s1  = min(dst1(1:end/2));
s10 = min(dst10(1:end/2));
sh  = min(dsth(1:end/2));
sH  = min(dstH(1:end/2));
sR  = min(dstR(1:end/2));


    
ubv2  = velTransVS(ubVD,1./Resb);
u01v2 = velTransVS(ubVD,1./Res01);
u1v2  = velTransVS(ubVD,1./Res1);
u10v2 = velTransVS(ubVD,1./Res10);
uhv2  = velTransVS(ubVD,1./Resh);
uHv2  = velTransVS(ubVD,1./ResH);
uRv2  = velTransVS(ubVD,1./ResR);

Ud01 = interp1(y,u01v2 - ubv2,1.0,'pchip');
Ud1  = interp1(y,u1v2  - ubv2,1.0,'pchip');
Ud10 = interp1(y,u10v2 - ubv2,1.0,'pchip');
Udh  = interp1(y,uhv2  - ubv2,1.0,'pchip');
UdH  = interp1(y,uHv2  - ubv2,1.0,'pchip');
UdR  = interp1(y,uRv2  - ubv2,1.0,'pchip');
% 
% m01 = s01.*exp(1./4).*pi^(1./2)./Ud01;
% m1  = s1 .*exp(1./4).*pi^(1./2)./Ud1;
% m10 = s10.*exp(1./4).*pi^(1./2)./Ud10;
% mh  = sh .*exp(1./4).*pi^(1./2)./Udh;
% mH  = sH .*exp(1./4).*pi^(1./2)./UdH;
% mR  = sR .*exp(1./4).*pi^(1./2)./UdR;



m01 = 1./y(s01 == dst01);
m1  = 1./y(s1  == dst1);
m10 = 1./y(s10 == dst10);
mh  = 1./y(sh  == dsth);
mH  = 1./y(sH  == dstH);
mR  = 1./y(sR  == dstR);

F01 = s01.*exp(-log(y.*m01).^2);
F1  = s1.*exp(-log(y.*m1).^2);
F10 = s10.*exp(-log(y.*m10).^2);
Fh  = sh.*exp(-log(y.*mh).^2);
FH  = sH.*exp(-log(y.*mH).^2);
FR  = sR.*exp(-log(y.*mR).^2);

u01vd = integralX4(Res01.*(str01 + F01),y);
u1vd  = integralX4(Res1 .*(str1  + F1),y);
u10vd = integralX4(Res10.*(str10 + F10),y);
uhvd  = integralX4(Resh .*(strh  + Fh),y);
uHvd  = integralX4(ResH .*(strH  + FH),y);
uRvd  = integralX4(ResR .*(strR  + FR),y);


figure(5)
hold on
plot(y(1:end/2),uhVD(1:end/2) ,'ro','Markersize',4);
plot(y(1:end/2),u01VD(1:end/2),'bd','Markersize',4);
plot(y(1:end/2),u1VD(1:end/2) ,'gs','Markersize',4);
plot(y(1:end/2),u10VD(1:end/2),'k^','Markersize',4);
plot(y(1:end/2),uRVD(1:end/2) ,'--','Markersize',4);
plot(y(1:end/2),uHVD(1:end/2) ,'-.','Markersize',4);

plot(y(1:end/2),uhvd(1:end/2) ,'ro','Markersize',4);
plot(y(1:end/2),u01vd(1:end/2),'bd','Markersize',4);
plot(y(1:end/2),u1vd(1:end/2) ,'gs','Markersize',4);
plot(y(1:end/2),u10vd(1:end/2),'k^','Markersize',4);
plot(y(1:end/2),uRvd(1:end/2) ,'--','Markersize',4);
plot(y(1:end/2),uHvd(1:end/2) ,'-.','Markersize',4);
set(gca,'xscale','log')




figure(6)
hold on
plot(y(1:end/2),dsth(1:end/2) ,'ro','Markersize',4);
plot(y(1:end/2),dst01(1:end/2),'bd','Markersize',4);
plot(y(1:end/2),dst1(1:end/2) ,'gs','Markersize',4);
plot(y(1:end/2),dst10(1:end/2),'k^','Markersize',4);
plot(y(1:end/2),dstR(1:end/2) ,'--','Markersize',4);
plot(y(1:end/2),dstH(1:end/2) ,'-.','Markersize',4);

plot(y(1:end/2),Fh(1:end/2) ,'ro','Markersize',4);
plot(y(1:end/2),F01(1:end/2),'bd','Markersize',4);
plot(y(1:end/2),F1(1:end/2) ,'gs','Markersize',4);
plot(y(1:end/2),F10(1:end/2),'k^','Markersize',4);
plot(y(1:end/2),FR(1:end/2) ,'--','Markersize',4);
plot(y(1:end/2),FH(1:end/2) ,'-.','Markersize',4);
set(gca,'xscale','log')





