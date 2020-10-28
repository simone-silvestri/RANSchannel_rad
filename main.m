%**************************************************************************
%
% RANS solver for fully developed turbulent channel with varying properties
% 
%       Created on: Nov 22, 2017
%          Authors: Rene Pecnik         (R.Pecnik@tudelft.nl)
%                   Gustavo J. Otero R. (G.J.OteroRodriguez@tudelft.nl)
%                   Process & Energy Department, Faculty of 3mE
%                   Delft University of Technology, the Netherlands.
%       Literature: Otero et al., 2017. Heat and fluid flow
% Last modified on: Jan 07, 2018
%               By: Rene Pecnik 
%**************************************************************************
% 
close all
clear
clc

%--------------------------------------------------------------------------
%       Include folders
%--------------------------------------------------------------------------

addpath('mesh');                % functions for the mesh
addpath('function');            % functions for the solver
addpath('turbmodels');          % functions of the turbulence models
addpath('fluxmodels');          % functions of the turbulent flux models
addpath('radiation');           % functions for the radiative calculations




%% define the models
% % 
% model(1).tmod = 'SA';
% model(1).fmod = {'NO'};
% model(1).cases= {'b','t01','t1','t5','t10','t20','br','t01r','t1r','t10r','H2O','CO2','H2OP','Rb','H2OR'};
% model(1).rmod = zeros(length(model(1).cases),1);
% model(1).kmod = zeros(length(model(1).cases),1);
% model(1).dens = [zeros(6,1); ones(7,1); 2*ones(2,1)];
% 
% 
% model(1).cases= {'b','br','Rb','t01','t1','t5','t10','t20','t01r','t1r','t10r',...
%     't01','t1','t5','t10','t20','t01r','t1r','t10r'};
% model(1).rmod = [zeros(11,1);...
%     ones(8,1)];
% model(1).kmod = [zeros(11,1);...
%     zeros(8,1)];
% model(1).dens = [0 1 2 0 0 0 0 0 1 1 1 ...
%     0 0 0 0 0 1 1 1];
% model.smod = 2*ones(19,1);
% 
% model(1).tmod = 'V2F';
% model(1).fmod = {'DWX'};

model(1).tmod = 'V2F';
model(1).fmod = {'NO'};
model(1).cases= {'Rb'};
model(1).rmod = 1;
model(1).kmod = 0;
model(1).dens = 2;
model.smod    = 2;


visual = 1;
                               
% STARTING LOOP!
for mdel = 1:length(model)
    for flu = 1:length(model(mdel).fmod)
        for cas = 1:length(model(mdel).cases)


%% ------------------------------------------------------------------------
%
%    User defined inputs
%  ------------------------------------------------------------------------
%


% -----  choose turbulence model  -----
% 'Cess'... Cess, R.D., "A survery of the literature on heat transfer in 
%           turbulent tube flow", Tech. Rep. 8-0529-R24, Westinghouse, 1958.
% 'SA'  ... Spalart, A. and Allmaras, S., "One equation turbulence model for 
%           aerodynamic flows", Recherche Aerospatiale-French edition, 1994.
% 'MK'  ... Myong, H.K. and Kasagi, N., "A new approach to the improvement of
%           k-epsilon turbulence models for wall bounded shear flow", JSME 
%           Internationla Journal, 1990.
% 'SST' ... Menter, F.R., "Zonal Two equation k-omega turbulence models for 
%           aerodynamic flows", AIAA 93-2906, 1993.
% 'V2F' ... Medic, G. and Durbin, P.A., "Towards improved prediction of heat 
%           transfer on turbine blades", ASME, J. Turbomach. 2012.
% 'DNS' ... without turbulence model; k,e,u taken from DNS solution
% 'NO'  ... without turbulence model; laminar
turbMod = model(mdel).tmod;


% -----  choose flux model  -----
% '1EQ'... Cess, R.D., "A survery of the literature on heat transfer in 
%           turbulent tube flow", Tech. Rep. 8-0529-R24, Westinghouse, 1958.
% '2EQ'  ... Spalart, A. and Allmaras, S., "One equation turbulence model for 
%           aerodynamic flows", Recherche Aerospatiale-French edition, 1994.
% 'PRT'  ... Myong, H.K. and Kasagi, N., "A new approach to the improvement of
%           k-epsilon turbulence models for wall bounded shear flow", JSME 
%           Internationla Journal, 1990.
turbPrT = model(mdel).fmod{flu};

% 0 ...  constant rho
% 1 ...  variable rho
% 2 ...  rho from DNS
varDens= model(mdel).dens(cas);
underrelaxU=0.7;

% -----  compressible modification  -----
% 0 ... Conventional models without compressible modifications
% 1 ... Otero et al.
% 2 ... Catris, S. and Aupoix, B., "Density corrections for turbulence
%       models", Aerosp. Sci. Techn., 2000.  
compMod = 0;

% -----  solve energy equation  ----- %%
% 0 ... energy eq not solved, density and viscosity taken from DNS
% 1 ... energy eq solved with under relaxation underrelaxT 
solveEnergy = 1;
underrelaxT = 0.9;

% -----  solve radiation equation  ----- 
% 0 ... radiation eq not solved (non radiative channel)
% 1 ... radiation eq solved every stepRad energy equation iterations
% 2 ... radiative heat source taken from DNS calculations (radCase)
solveRad = model(mdel).smod(cas);
radCase  = model(mdel).cases{cas};
underrelaxQ = 0.1;

% -----  choose Radiation model modification -----
% 0 ...  Conventional t2 - et equations
% 1 ...  Radiative source term in t2 and et equations
RadMod  = model(mdel).rmod(cas);
stepRad = 1;
startRad = 1000;
kPMod   = model(mdel).kmod(cas); % has to be 0 for using kP and 1 for using kG and 2 kG simplified 

% -----  channel height  -----
height = 2;

% -----  number of mesh points  -----
n = 100;
% -----  number of angular points  -----
Tc     = 573;
Th     = 955;
T0     = 1.5;
Pr     = 1.0;

if strcmp(radCase,'t10p') || strcmp(radCase,'t1p')
    Pr = 0.7;
end

if varDens==1
    ReT    = 3750;
else
    ReT    = 2900;
end

if strcmp(radCase,'H2OR') || strcmp(radCase,'Rb')
    ReT = 16700;
    T0  = 0.5;
    Tc  = 600;
    Th  = 1800;
    Pr  = 0.93;
end

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
fact = 5;
ns = 1;

% -----  Parameter definition
Pl  = 0.03;
Prt = ones(n,1)*1.0; 
b_old = -ones(n-2,1)*0.005;
b = b_old;
switch varDens
    case 0; casename = 'constant';
    case 1; casename = 'vardens';          
    case 2; casename = 'varvisc';
    case 3; casename = 'nodensity';
end

%% ------------------------------------------------------------------------
%
%      Generate mesh and angular mesh
%
[MESH] = mesh(height, n, fact, ns, discr);

%% ------------------------------------------------------------------------
%
%       Solve RANS 
%

u      = zeros(n,1);
uT     = zeros(n,1);
T      = 1 - MESH.y/height;
T(1)   = 1;
T(end) = 0;

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
Qdns  = QR;
Em    = interp1(Rm(:,1),Rm(:,4),MESH.y,'pchip');
Edns  = Em;
G     = interp1(Rm(:,1),Rm(:,5),MESH.y,'pchip');
Gdns = G;
qy    = integralX(MESH.y,QR);
kG    = kP;

%if strcmp(radCase,'H2OR')
%    cr33 = 22;
%    cr22 = 22;
%     WV = interp1(Dm(:,1),WV,MESH.y,'pchip');
%     WV   = ((cr33-cr22).*MESH.y.^2 - 2*(cr33-cr22).*MESH.y +cr33);
%else
cr33 = 7.*(ReT/2900).^(3./4); %.*Pr.^(1./2);
cr22 = 7.*(ReT/2900).^(3./4); %.*Pr.^(1./2);
WV   = ((cr33-cr22).*MESH.y.^2 - 2*(cr33-cr22).*MESH.y +cr33);

%end

% 
% if kPMod==1
%     temp = importdata(radg);
%     kG   = interp1(temp(:,1),temp(:,2),MESH.y,'pchip');
% end
% if kPMod==2
%     F    = interp1(Rm(:,1),Rm(:,6),MESH.y,'pchip');
%     for i=1:length(MESH.y)
%         fun = @(x) x./WV(i).*atan(WV(i)./x)-F(i);
%         kG(i) = fzero(fun,1);
%     end
% end


cP1  = polyfit(1./Tdns,kP,5);
cP3  = polyfit(1./Tdns,kG,5);
for i = 1:6
    cP(6-i+1) = cP1(i);
    cP2(6-i+1)= cP3(i);
end

Tdns = (Tdns -Tc)/(Th-Tc);

if solveRad == 3
    if RadMod==1
        t2dns = interp1(Dc(:,1),Dc(:,2),MESH.y,'pchip');
        Em    = 4.*(Tdns./T0+1).^4 + t2dns.*(24.*Tdns.^2./T0.^4 + 48./T0.^3.*Tdns + 24./T0.^2);
        G     = rad_analytical(Em,MESH.y,kP,Tc,Th);
        QR    = kP.*(Em-G);
    else
        Em = 4.*(Tdns./T0+1).^4;
        G = rad_analytical(Em,MESH.y,kP,Tc,Th);
        QR = kP.*(Em-G);
    end
end


if (solveEnergy==0)  % transport properties
    r=ones(n,1); mu=ones(n,1)/ReT; T=zeros(n,1);
else
    [r, mu, alpha] = calcProp(T, Dm, ReT, Pr, casename, T0, MESH.y);
end

% turbulent scalars
t2   = 0.1*ones(n,1); t2(1)= 0.1; t2(n)= 0.1;
k    = 0.1*ones(n,1); k(1) = 0.0; k(n) = 0.0;
e    = 0.001*ones(n,1);
et   = 0.001*ones(n,1);
er   = 0.001*ones(n,1);
v2   = 2/3*k;
om   = ones(n,1);
mut  = zeros(n,1);
alphat  = zeros(n,1);
nuSA = mu./r; nuSA(1) = 0.0; nuSA(n) = 0.0;

%% SOLVE THE EQUATIONS
%--------------------------------------------------------------------------
%
%       Iterate RANS equations
%
nmax   = 1000000;   tol  = 1.e-8;  % iteration limits
nResid = 100;                      % interval to print residuals

residual = 1e20; residualT = 1e20; residualQ = 1e20; iter = 1;

calckG  = 1200;

figure('Position',[1000 1000 800 600])
figure('Position',[50 1000 800 600])

if solveRad == 1
    QR = 0.0*QR;
end
while (residual > tol || residualT > tol || residualQ > tol*1e2) && (iter<nmax)
 
    if solveRad == 1 && iter >= startRad
        nResid=1;
    end
    
    if kPMod == 2 && mod(iter,calckG)==0
        ls =real(ls);
        T  =real(T);
        kG = fixkG(T,radCase,ls,kP);
    end

    if(visual==1)
        if(mod(iter,nResid)==0)
            figure(1)
            clf;
            plot(MESH.y,T,MESH.y,Tdns);
            figure(2)
            semilogy(iter/20,residualT,'bo','MarkerSize',6)
            ylim([1e-09 1])
            hold on
            pause(0.00001);
        end
    end

    % Solve turbulence model to calculate eddy viscosity
    switch turbMod
        case 'V2F';   [k,e,v2,mut] = V2F(u,k,e,v2,r,mu,MESH,compMod);
        case 'MK';    [k,e,mut]    = MK(u,k,e,r,mu,ReT,MESH,compMod,iter);
        case 'SST';   [k,om,mut]   = KOmSST(u,k,om,r,mu,MESH,compMod);
        case 'SA';    [nuSA,mut]   = SA(u,nuSA,r,mu,MESH,compMod);
        case 'Cess';  mut          = Cess(r,mu,ReT,MESH,compMod);
        case 'DNS';   [k,e,u,mut,r]= DNS(Dk,Dm,mu,ReT,MESH,kPMod);
        otherwise;	  mut          = zeros(n,0.9);
    end
    
    % Solve energy equation
    if (solveEnergy == 1)

        % Solve turbulent flux model to calculate eddy diffusivity 
        switch turbPrT
            case 'V2T'; [uT,lam] = V2T_new(uT,k,e,v2,mu,alpha,ReT,Pr,Pl,T,kP,kG,r,MESH,RadMod,Em,G,Th,Tc,T0,cP,cP2,WV);
            case 'PRT'; [lam,alphat] = PRT( mu,mut,alpha,T,r,qy,ReT,Pr,Pl,MESH,RadMod);
            case 'DWX'; [lam,t2,et,er,alphat,ls,fk,fe,fg] = DWX_rad(T,r,u,t2,et,er,k,e,alpha,mu,ReT,Pr,Pl,MESH,Em,G,kP,kG,RadMod,Th,Tc,T0,cP,cP2,alphat,WV);
            case 'DWV'; [lam,t2,et,alphat] = DWV(T,r,t2,et,k,e,alpha,MESH,v2);
            otherwise;  lam = alpha + (mut./0.9);
        end
%         T  =real(T);
%         et =real(et);
%         t2 =real(t2);
        
        T_old = T;
        
        % diffusion matrix: lam*d2phi/dy2 + dlam/dy dphi/dy
        A =    bsxfun(@times, lam, MESH.d2dy2) ... 
             + bsxfun(@times, MESH.ddy*lam, MESH.ddy);
         
         turbQ = MESH.ddy*(r.*uT);
         % source term: -Qt/RePr + duT/dy
         switch turbPrT
             case 'V2T'; b = QR(2:n-1)./(ReT*Pr*Pl) + turbQ(2:n-1);
             otherwise;  b = QR(2:n-1)./(ReT*Pr*Pl);
         end
         
         % Solve
         T = solveEq(T,A,b,underrelaxT);
         
         residualT = norm(T - T_old);
    else
        residualT = 0;
    end
    
    if solveRad == 1 && iter >= startRad
        if RadMod==1
            Em = 4.*(T./T0+1).^4 + t2.*(24.*T.^2./T0.^4 + 48./T0.^3.*T + 24./T0.^2);
            if mod(iter,stepRad)==0
                G = rad_analytical(Em,MESH.y,kP,Tc,Th);
                Q_old = QR;
            end
            QR = kP.*(Em-G) + t2.*(fk.*(fe-fg));
            residualQ = norm(QR(2:n-1)-Q_old(2:n-1));
        else
            Em = 4.*(T./T0+1).^4;
            if mod(iter,stepRad)==0
                G = rad_analytical(Em,MESH.y,kP,Tc,Th);
                Q_old = QR;
            end
            QR = kP.*(Em-G) + t2.*(fk.*(fe-fg));
            residualQ = norm(QR(2:n-1)-Q_old(2:n-1)); 
        end
    else
        residualQ = 0.0;
    end
    
    if strcmp(turbMod,'DNS')==0
        
        % calculate density
        [r,mu,alpha] = calcProp(abs(T), Dm, ReT, Pr, casename, T0, MESH.y);
         
        % Solve momentum equation:  0 = d/dy[(mu+mut)dudy] - rho fx
        mueff = mu + mut;
        
        % diffusion matrix: mueff*d2phi/dy2 + dmueff/dy dphi/dy
        A =   bsxfun(@times, mueff, MESH.d2dy2) ...
            + bsxfun(@times, MESH.ddy*mueff, MESH.ddy);
        
        % Right hand side
        %b = -ones(n-2,1);
        bulk  = trapz(MESH.y,u)/2;
        b     = (0.99*b_old + (bulk-1)*0.005);
        
        % Solve
        u_old = u;
        u = underrelaxU*u+(1-underrelaxU)*solveEq(u,A,b,1);
        residual = norm(u-u_old);
        
        b_old = b;
    else
        residual = 0;
    end
    
    if(visual==1) || (visual==0)
        % Printing residuals
        if (mod(iter,nResid) == 0)
            fprintf('%d\t%12.6e\t%12.6e\t%12.6e\n', iter, residual, residualT, residualQ);
        end
    end
    iter = iter + 1;
end
%fprintf('%d\t%12.6e\n\n', iter, residual);
str =['Finished:','  ',model(mdel).tmod,'  ','with','  ',...
    model(mdel).fmod{flu},'  ','on case','  ',model(mdel).cases{cas},'  ',...
    'and radmod','  ',num2str(model(mdel).rmod(cas))];
fprintf('%s\n\n',str);

%fprintf('The residuals with this method is: %12.6e\n\n',norm(T-Tdns));

%% ------------------------------------------------------------------------
%
%       plot RANS results and compare with DNS
%
%  Plot function: plot_(...), plot_log() one side with log scale, and
%  plot_loglog(...) both side with log scale.
%  plotit_*(x,y,i,flagHold,style, xaxis, yaxis) 
%  x: x-axis variable, y: y-axis variable, i: figure numb, style: color-
% line type, xaxis: label x-axis, and yaxis: label y-axis

% For plots with latex style
set(groot, 'DefaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

y = MESH.y;

% %% ------------------------------------------------------------------------
% % plotting the velocity profiles
% 
% first calculate uplus (not really necessary, since utau = 1.0 already)
dudy   = MESH.ddy*u;
utau   = sqrt(mu(1)*dudy(1)/r(1));
upl    = u/utau;
utD    = Dm(:,4)./utau;

%% ------------------------------------------------------------------------
% plotting the DNS profiles
% 
if strcmp(turbPrT,'NO')
    alphat = mut./0.9./r;
end
alphatd = Df(:,2)./(-Df(:,3))/ReT; 

THF     = (r.*alphat.*(-MESH.ddy*T));
dT = MESH.ddy*T;
if strcmp(turbPrT,'V2T')
    THF = r.*uT;
    alphat = -uT./dT;
end
TST     = mut.*(-dudy)./r;


Pt = 2*r.*alphat.*dT.^2;

diffT = T - Tdns;
THFdns = (interp1(Df(:,1),Df(:,2),y,'pchip'));
close all
if(visual==1)
    figure(1); hold off
    plot(y,T); hold on
    plot(y,Tdns);
    plot(y,diffT);
    
    figure(2); hold off
    plot(y,THF); hold on
    plot(y,THFdns);
    
    figure(3); hold off
    plot(y,alphat); hold on
    plot(Df(:,1),alphatd);
    
    figure(4); hold off
    plot(y,u); hold on
    plot(Dm(:,1),Dm(:,4));
end


alphatdns = interp1(Df(:,1),alphatd,y,'pchip');


%% Calculation of conductive heat transfer

CHF = alpha.*(-MESH.ddy*T);
CHFD = alpha.*(-MESH.ddy*Tdns);
diff = (CHF - CHFD)./CHFD;

%fprintf('The HF in this case is: %12.6e %12.6e\n\n',diff(1),diff(end));

%% Calculation of the Nusselt number

udns = interp1(Dm(:,1),Dm(:,4),MESH.y,'pchip');
Tb = trapz(y,Tdns.*udns)/trapz(y,udns);

NuD = -1./(Tb - 1).*interp1(MESH.y,CHFD,1.0,'pchip');
Tb = trapz(y,T.*u)/trapz(y,u);
NuR = -1./(Tb - 1).*interp1(MESH.y,CHF,1.0,'pchip');
%fprintf('The Nusselt number in this case is: %12.6e while for RANS: %12.6e\n\n',NuD,NuR);

%% ------------------------------------------------------------------------
%saving down the files
if solveRad==1 || solveRad==3
    switch RadMod
        case 2;   string = strcat('solution/calcQ/',turbMod,'/',turbPrT,'/',radCase,'_','rad_',num2str(kPMod));
        case 1;   string = strcat('solution/calcQ/',turbMod,'/',turbPrT,'/',radCase,'_','rad_',num2str(kPMod));
        case 0;   string = strcat('solution/calcQ/',turbMod,'/',turbPrT,'/',radCase);
    end
else
    switch RadMod
        case 2;   string = strcat('solution/',turbMod,'/',turbPrT,'/',radCase,'_','rad_',num2str(kPMod));
        case 1;   string = strcat('solution/',turbMod,'/',turbPrT,'/',radCase,'_','rad_',num2str(kPMod));
        case 0;   string = strcat('solution/',turbMod,'/',turbPrT,'/',radCase);
    end
end
if strcmp(turbPrT,'NO')
    ls = zeros(n,1);
    fk = zeros(n,1);
    fe = zeros(n,1);
    fg = zeros(n,1);
end
fid = fopen(string,'w');
for i=1:n
    fprintf(fid,'%20.10f\t%20.10f\t%20.10f\t%20.10f\t%20.10f\t%20.10f\t%20.10f\t%20.10f\t%20.10f\t%20.10f\t%20.10f\t%20.10f\t%20.10f\t%20.10f\t%20.10f\t%20.10f\t%20.10f\t%20.10f\t%20.10f\t%20.10f\t%20.10f\t%20.10f\n',...
        y(i),u(i),T(i),diffT(i),mut(i),alphat(i),THF(i),k(i),e(i),v2(i),t2(i),et(i),r(i),CHF(i),CHFD(i),ls(i),kG(i),kP(i),fk(i),fe(i),fg(i),Pt(i));
end
fclose(fid);
fid = fopen(strcat(string,'_bulks'),'w');
fprintf(fid,'%20.10f\t%20.10f\t%20.10f\t%20.10f\t%20.10f\n',...
    norm(T-Tdns),diff(1),diff(end),NuD,NuR);
fclose(fid);

if visual==0 
    close all;
end
%% CLOSING THE LOOP ON THE MODELS
        end
    end
end

