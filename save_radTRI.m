clear all
close all

radCase  = 't10p';
RadMod   = 1;
Tc = 573;
Th = 955;
T0 = 1.5;
n=100;
fact = 5;
ns = 1;
discr = 'finitediff';
[MESH] = mesh(2.0, n, fact, ns, discr);
y = MESH.y;
% radiative quantities
pathm    = strcat('solution/DNS/',radCase,'m');
radm     = strcat('solution/DNS/radiation/',radCase,'m');
radg     = strcat('solution/DNS/radiation/',radCase,'kg');
pathf    = strcat('solution/DNS/',radCase,'f');
pathc    = strcat('solution/DNS/',radCase,'c');
pathk    = strcat('solution/DNS/',radCase,'k');
pathr    = strcat('solution/DNS/',radCase,'r');
Rm = importdata(radm);
Dm = importdata(pathm);
Df = importdata(pathf);
Dc = importdata(pathc);
Dk = importdata(pathk);
Dr = importdata(pathr);
cP = zeros(6,1);
Tdns = (Th-Tc)*interp1(Dm(:,1),Dm(:,5),MESH.y,'pchip')+Tc;

kP    = interp1(Rm(:,1),Rm(:,2),MESH.y,'pchip');
QR    = interp1(Rm(:,1),Rm(:,3),MESH.y,'pchip');
Q_old = QR;
Em    = interp1(Rm(:,1),Rm(:,4),MESH.y,'pchip');
E_old = Em;
G     = interp1(Rm(:,1),Rm(:,5),MESH.y,'pchip');
G_old = G;
qy    = integralX(MESH.y,QR);
kG    = kP;


Tdns = (Tdns -Tc)/(Th-Tc);
t2dns = interp1(Dc(:,1),Dc(:,2),MESH.y,'pchip');
if RadMod==1
    Em = 4.*(Tdns./T0+1).^4 + t2dns.*(24.*Tdns.^2./T0.^4 + 48./T0.^3.*Tdns + 24./T0.^2);
else
    Em = 4.*(Tdns./T0+1).^4;
end

G = rad_analytical(Em,MESH.y,kP,Tc,Th);

save(strcat('radiation/noTRI/G',radCase,'_',num2str(RadMod),'.mat'),'G')

