function G = rad_analytical(E,r,k,Tc,Th)
%UNTITLED14 Summary of this function goes here
%   Detailed explanation goes here


m = length(E);

Q = zeros(m,1);
G = zeros(m,1);

Iw1=(Th./Tc).^4;
Iw2=(Tc./Tc).^4;

abs = k;

tau = r.*abs;

%tau(1) = 1e-20;

tauL = abs(end)*2.0;
%tau(end) = tauL - 1e-20;

% Tt    = @(x) 955 + (573-955)*x./10.*0.5;

n=500;
fact = 15;
% y - coordinate: tanh clustering at the walls
i = ((0:n-1)')/(n-1) - 0.5;
y = 2.0*(1.0 + tanh(fact*i)/tanh(fact/2))/2.0;
y = (y + (2-y(end:-1:1)))/4.0;

E = E./4;
%Tt = @(x) -(955-573).*x.*0.5./abs(end)+955;

%% CALCULATING ABSORPTION


E2_tau  = zeros(m,1);
E2_tau2 = zeros(m,1);
Int1    = zeros(m,1);
Int2    = zeros(m,1);
% Int1t   = zeros(m,1);
% Int2t   = zeros(m,1);
% Gt      = zeros(m,1);

for i=1:m

   E2_tau(i)=exp(-tau(i))-tau(i)*igamma(0,tau(i));
   E2_tau2(i)=exp(-(tauL-tau(i)))-(tauL-tau(i))*igamma(0,tauL-tau(i));    

   Int1(i) = 0;
   Int2(i) = 0;

   xtemp1  = tau(i).*y;
   xtemp2  = (tauL-tau(i)).*y+tau(i);
   
   T1     = interp1(tau,E,xtemp1,'pchip');
   fun1   = T1.*igamma(0,tau(i) - xtemp1);
   fun1(fun1==inf) = T1(fun1==inf).*igamma(0,1e-100);
   Int1(i)= trapz(xtemp1,fun1);
   
 
   T2     = interp1(tau,E,xtemp2,'pchip');
   fun2   = T2.*igamma(0,xtemp2 - tau(i));
   fun2(fun2==inf) = T2(fun2==inf).*igamma(0,1e-100);
   Int2(i)= trapz(xtemp2,fun2);

%    fun1t=@(x) Tt(x).^4.*igamma(0,tau(i)-x);
%    Int1t(i) = integral(fun1t,0,tau(i),'RelTol',1e-08);
%    fun2t=@(x) Tt(x).^4.*igamma(0,x-tau(i));
%    Int2t(i) = integral(fun2t,tau(i),tauL,'RelTol',1e-08);
%    
%    Int1t(i) = Int1t(i)*2*stefan;
%    Int2t(i) = Int2t(i)*2*stefan;
%    grey isothermal walls
    G(i) = 2*(Iw1*E2_tau(i)+Iw2*E2_tau2(i)+Int1(i)+Int2(i));
%    Gt(i)= 2*(stefan*Tw1^4*E2_tau(i)+stefan*Tw2^4*E2_tau2(i))+(Int1t(i)+Int2t(i));
    if mod(i,10)==0 
        disp(i)
    end
   
end

end

