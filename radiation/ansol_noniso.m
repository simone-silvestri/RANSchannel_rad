clear all
close all
clc

stefan=5.670373e-08;
Tw1=1000;
Tw2=500;
r = 0.005:0.005:1.995;
abs=5.;
tau =r*abs;

T = @(t) 500+500*(1-0.5*t/abs);

tauL = abs*2.0;

Lx=2.;


eps=1.;
divQ=zeros(1,length(tau));
dtau = tau(end)-tau(end-1);

%% CALCULATING WALL RADIOSITIES

Int3=0;
fun = @(t) T(t).^4.*(exp(-t)-t.*igamma(0,t));
Int3 = integral(fun,0,tauL,'RelTol',1e-08);

Int3 = Int3 * 2 * (1-eps) * stefan;
Int4 = stefan * eps * Tw1^4;

   fun3=@(x) (exp(-tauL./x)).*x;
   E3=integral(fun3,0,1,'RelTol',1e-08);

deno = (1 + (1-eps) * 2*E3);

Jw = (Int3+Int4)/deno;

%% CALCULATING ABSORPTION

for i=1:length(tau)

   E2_tau(i)=exp(-tau(i))-tau(i)*igamma(0,tau(i));
   E2_tau2(i)=exp(-(tauL-tau(i)))-(tauL-tau(i))*igamma(0,tauL-tau(i));    

   Int1(i) = 0;
   Int2(i) = 0;
   
   fun1=@(x) T(x).^4.*igamma(0,tau(i)-x);
   Int1(i) = integral(fun1,0,tau(i),'RelTol',1e-08);
   fun1=@(x) T(x).^4.*igamma(0,x-tau(i));
   Int2(i) = integral(fun1,tau(i),tauL,'RelTol',1e-08);
   
   Int1(i) = Int1(i)*2*stefan;
   Int2(i) = Int2(i)*2*stefan;
   
   % grey isothermal walls
    G(i) = 2*(stefan*Tw1^4*E2_tau(i)+stefan*Tw2^4*E2_tau2(i))+(Int1(i)+Int2(i));
    
    disp(i)
    
end
divQ = 4*stefan*abs*T(tau).^4 - abs*G;
rad = [r divQ];
plot(rad(:,1),rad(:,2))

