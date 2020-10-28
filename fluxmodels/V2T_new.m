%flux model for turbulent heat transfer
% Lai & So "Near-wall modeling of turbulent heat transfer"
%   
%   uT-eq :
%   0 = Pt + Tt - et + Ft
%
%   Pt : production of turbulent heat flux 
%        = v2 dTdy
%
%   Tt : turbulent transport + viscous diffusion
%        = d/dy((cs0 2 kv2/e + 1/Re + 1/3(Pr+1)/RePr) duTdy )
%
%        = d/dy(a) duTdy + a d2uTdy
%   where  a =  cs0 2 kv2/e + 1/Re + 1/3(Pr+1)/RePr
%
%   et : viscous + molecular dissipation
%        = fw (1+1/Pr) e/k uT
%
%   Ft : pressure term (diffusion + strain)
%        = -c10 e/k uT + fw (c10 - 1) e/k uT
% Input:
%   u           velocity
%   k           turbulent kinetic energy, from previous time step
%   e           turbulent kinetic energy dissipation rate per unit volume,  
%               from previous time step
%   v2          wall normal velocity fluctuation, from previos time step
%   mu          molecular viscosity
%   mesh        mesh structure
%
% Output:
%   uT          turbulent heat flux
%   Prt         turbulent Prandtl number

function [ uT, lam ] = V2T_new(uT,k,e,v2,mu,alpha,ReT,Pr,Pl,T,kP,kG,r,mesh,RadMod,Em,G,Th,Tc,T0,cP,cP2,WV)


    lam = alpha;
    n = size(T,1);
    y = mesh.y;
    underrelaxuT = 0.8;
    cs0 = 0.2;
    c10 = 3.0;
    wallDist = min(y,2-y);

    Ret = r.*k.^2./(mu.*e);
    fw = exp(-(Ret./80).^2);
          
    %% -----------------------GENERAL RADIATIVE MODIFICATION TO alphat AND functions fem,fk,fg
    
    % Radiative Planck mean absorption coefficient
    Tr = T*(Th-Tc) + Tc;
    fk = (cP(6).*5./Tr.^6+cP(5)*4./Tr.^5+cP(4)*3./Tr.^4+cP(3)*2./Tr.^3 ...
        +cP(2)*1./Tr.^2 );
    fk = -(Th-Tc)*fk/(ReT*Pr*Pl);
    fkg= (cP2(6).*5./Tr.^6+cP2(5)*4./Tr.^5+cP2(4)*3./Tr.^4+cP2(3)*2./Tr.^3 ...
        +cP2(2)*1./Tr.^2 );
    fkg = -(Th-Tc)*fkg/(ReT*Pr*Pl);
        
    fem  = (16/T0^4*T.^3 + 48/T0^3*T.^2 + 48/T0^2*T + 16/T0)/(ReT*Pr*Pl);
  
    fg  = (fem.*(kG)+fkg.*(Em-G))./WV.*atan(WV./(kG));

    %%
    

    %Production and diffusion terms (diffusion implicit)
    Pt     = v2.* (mesh.ddy*T);
    Dturb  = cs0*2.0.*k.*v2./e;
    Dvisc  = mu + 1./(wallDist+2).*(alpha-mu); 
 
    a = Dturb + Dvisc;
    
    % diffusion matrix: a*d2()/dy2 + da/dy d()/dy
    A  =   bsxfun(@times, a, mesh.d2dy2) ... 
       +   bsxfun(@times, (mesh.ddy*a), mesh.ddy);
    
    % implementation of dissipation
    for i=2:n-1
      A(i,i) = A(i,i) - 0.5*(1 + 1./Pr) * e(i)/k(i) * fw(i);
    end
    % implementation of pressure term
    P1 = -c10*e./k;
    P2 = 0.0;
    Pr = 0.0;
    Pw = c10*e./k;
    
    phi = P1+P2+Pr+fw.*Pw;
    for i=2:n-1
        A(i,i) = A(i,i) + phi(i); %%(c10*(fw(i)-1)-fw(i)) * e(i)/k(i)*(1 - kP(i).*(fem(i)-fg(i)) - (Em(i)-G(i)).*fk(i));  %CHANGED: SIMONE!
    end
    % Addition of radiative terms
    if RadMod == 1
        for i=2:n-1
            A(i,i) = A(i,i) - kP(i).*(fem(i)-fg(i)) - (Em(i)-G(i)).*fk(i);
        end
    end
    b = Pt(2:n-1);
   
    % Solve
    uT = solveEq(uT,A,b,underrelaxuT);
    uT(2:n-1) = max(uT(2:n-1), 1.e-12);
    
end
