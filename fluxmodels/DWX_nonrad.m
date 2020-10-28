%       Implementation of the dwx model
%       Reference,
%       M. Karcz and J. Badur., "A TURBULENT HEAT FLUX TWO–EQUATION
%       θ02 –εθ CLOSURE BASED ON THE V2F TURBULENCE MODEL"
%**************************************************************************
%
% Conventional models without radiation modifications:
%    t2-eq:  0 = 2 Pt - 2 et + ddy[(alpha+r.*alphat/sigmat2) dt2dy] 
%    et-eq:  0 = Cp1 fp1 sqrt(e et/ (k t2)) Pt - Cd1 fd1 et^2 / t2 
%              - Cd2 fd2 e et / t2 + ddy[(alpha+r.*alphat/sigmaet)detdy] 
%
% Models with radiative modifications:
%    t2-eq:  0 = 2 Pt - 2 et + ddy[(alpha+alphat/sigmat2) dt2dy] - 2kP Emt
%               +2kP Gt
%    et-eq:  0 = Cp1 fp1 sqrt(e et/ (k t2)) Pt - Cd1 fd1 et^2 / t2 
%              - Cd2 fd2 e et / k + ddy[(alpha+alphat/sigmaet)detdy] 
%              - 2*kP dEm dy dt dy + 2*kP dG dy dt dy
%
%    alphat = Cl v2 k^l/e t2^m/et
%    Pt     = - 2 r.*alphat (dTdy)^2
%
% Radiation modification: introduction of an emission dissipation and
% absorption source in both t2 and et equations:
%
% alphat = Cl v2 k^l/e t2^m/et
%
% Modeling of radiation fluctuations:
%
% Em   = fem t
% G    = fem fg t
%
% Where:
%
% fem  = (16/T0^4 T^3 + 48/T0^3 T^2 + 48/T0^2 T + 16/T0) / (ReT Pr Pl)
% fg  = kP/cr22 atan(cr22/kP)
% cr22 = 7.0;
%
% Radiative term modeling in t2 equation:
%
%  - 2 kP fem t2 (1 - fg)
%
% Radiative term modeling in et equation:
%
%  - kP ( fem et + 0.5 d fem dy dt2 dy ) (1 - fg)       (if fg != f(y))
%
%
% Additional terms due to k prime
%  c   = [-0.23093, -1.12390*A, 9.41530*A^2, -2.99880*A^3, 0.51382*A^4, -1.8684e-05*A^5];
%  A   = 1000
%
%  k prime = - t prime dT (c(6) 5/T^6 + c(5) 4/T^5 + c(4) 3/T^4 + c(3) 2/T^3 + c(2) 1/T^2 )


function [ lam,t2,et,alphat ] = DWX_nonrad(T,r,u,t2,et,k,e,alpha,mu,ReT,mesh)

    n = size(T,1);
    
    tv       = (mesh.ddy*u).*mu;
    ut       = sqrt(tv(1)./r(1)/2-tv(end)./r(end)/2);
    Retau    = ReT.*ut.*r;
    y        = mesh.y;
   	wallDist = min(y, 2-y);
    yplus    = wallDist.*Retau; 

    % Model constants
    Cl    = 0.1;
    Cp1   = 2.34; 
    % Cd1   = 2.0;  %%% ORIGINAL: DENG
    Cd1   = 1.5;  %%% CHANGED: SIMONE
    Cd2   = 0.9;  % NOT OF ANY INTEREST (GETS ELIMINATED BY fD2)
    Ce2   = 1.9;
    siget = 1.0;
    sigt2 = 1.0;
    m     = 0.5;
        
    % Relaxation factors
    underrelaxt2  = 0.8;
    underrelaxet  = 0.8;  

    % Time and length scales, eddy diffusivity and turbulent production
    Reps   = (mu./r.*e).^(1./4)./mu.*r.*wallDist;
    Rturb  = r.*(k.^2)./(mu.*e);
    
    % Model damping functions
    fd1    = 1 - exp(-(Reps./1.7)).^2;
    % feps   = (1 - 0.3*exp(-(Rturb/6.5).^2)).*(1-exp(-Reps/3.1)).^2;  %%% ORIGINAL: DENG
    feps   = (1 - 0.3*exp(-(Rturb/6.5).^2));                           %%% CHANGED: SIMONE
    fw0    = exp( -(Rturb./80) .^2);
    fd2    = (1/Cd2)*(Ce2*feps - 1).*(1 - exp(-Reps./5.8).^2);
    
    % turbulent diffusivity and production
    R = 0.5*(t2./et.*e./k);

    %fl  = (1-exp(-Reps/16)).^2.*(1+3./(Rturb.^(3/4.)));       %%% ORIGINAL: DENG
    fl  = (fw0.*0.1./(Rturb.^(1./4))) + (1-exp(-yplus/30)).^2; %%% CHANGED: SIMONE
    fl(1:n-1:n) = 0.0;
    
    alphat = max(Cl.*fl.*k.^2./e.*(2*R).^m,0.0);
    
    Pt  = r.*alphat.*(mesh.ddy*T).^2;

    % ---------------------------------------------------------------------
    % et-equation
    %    0 = Cp1 fp1 sqrt(e et/ (k t2)) Pt - Cd1 fd1 et^2 / t2 
    %        - Cd2 fd2 e et / t2 + ddy[(alpha+r.*alphat/sigmaet)detdy]    
    
    % effective diffusivity
    lam = alpha + r.*alphat./siget;
    
    % diffusion matrix: lam*d2()/dy2 + dlam/dy d()/dy
    A =   bsxfun(@times, lam, mesh.d2dy2) ... 
        + bsxfun(@times, (mesh.ddy*lam), mesh.ddy);

    % implicitly treated source term : -Cd2 fd2 e/t2 - Cd1 fd1 et^2 / t2 
    for i=2:n-1
        A(i,i) = A(i,i) - Cd2*fd2(i)*e(i)/k(i).*r(i) - Cd1*fd1(i)*et(i)/t2(i).*r(i);
    end
        
    % Right-hand-side: - Cp1 fp1 sqrt(e et/ (k t2)) Pt
    b = - Cp1*sqrt(e(2:n-1).*(et(2:n-1))./k(2:n-1)./t2(2:n-1)).*Pt(2:n-1);
    
    % Boundary conditions
    et(1) = alpha(1)*(sqrt(t2(2)  )/(mesh.y(2)-mesh.y(1  ))).^2;
    et(n) = alpha(n)*(sqrt(t2(n-1))/(mesh.y(n)-mesh.y(n-1))).^2;
    
    %solve
    et = solveEq(et, A, b, underrelaxet);
    et(2:n-1) = max(et(2:n-1), 1.e-12);

    % ---------------------------------------------------------------------
	% t2-equation
    %    0 = 2 Pt - 2 et + ddy[(alpha+r.*alphat/sigmat2) dt2dy] - Rad
    
    % effective diffusivity
    lam = alpha + r.*alphat./sigt2;
        
    % diffusion matrix: lam*d2()/dy2 + dlam/dy d()/dy
    A =   bsxfun(@times, lam, mesh.d2dy2) ... 
        + bsxfun(@times, (mesh.ddy*lam), mesh.ddy);
    
    % implicitly treated source term
    for i=2:n-1
        A(i,i) = A(i,i) - 2*et(i)/t2(i).*r(i);
    end
    
    % right-hand side
    b = - 2*Pt(2:n-1);
    
    % Boundary conditions
    t2(1) = 0;
    t2(n) = 0;
    
    % Solve
    t2 = solveEq(t2, A, b, underrelaxt2);
    
    lam = alpha + r.*alphat;

end

