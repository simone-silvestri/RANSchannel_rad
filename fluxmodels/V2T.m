

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
%
%%
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

function [ uT, lam ] = V2T(uT,k,e,v2,mu,mut,ReT,Pr,T,mesh)


    n = size(T,1);
    lam = mu./Pr; % + mut./0.9;   

    underrelaxuT = 0.8;
    cs0 = 0.11;
    c10 = 4.0;
    
    Ret = k.^2./(mu.*e);
    fw = exp(-(Ret./80).^2);
    
    %Production and diffusion terms (diffusion implicit)
    Pt = v2.* (mesh.ddy*T);
    a  = cs0 * 2.0 * k.*v2./e + 1/ReT + 1./3.*(1-Pr)./ReT./Pr; 
    
    
    % diffusion matrix: a*d2()/dy2 + da/dy d()/dy
    A  =   bsxfun(@times, a, mesh.d2dy2) ... 
       +   bsxfun(@times, (mesh.ddy*a), mesh.ddy);
    
    % implementation of dissipation
    for i=2:n-1
        A(i,i) = A(i,i) - (1 + 1./Pr) * e(i)/k(i) * fw(i);
    end
    % implementation of pressure term
    for i=2:n-1
        A(i,i) = A(i,i) + (c10*(fw(i)-1)-fw(i)) * e(i)/k(i);
    end
    b = Pt(2:n-1);
   
    % Solve
    uT = solveEq(uT,A,b,underrelaxuT);
    uT(2:n-1) = max(uT(2:n-1), 1.e-12);
end

