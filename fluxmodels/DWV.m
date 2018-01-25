%       Implementation of the dwx model
%       Reference,
%       M. Karcz and J. Badur., "A TURBULENT HEAT FLUX TWO–EQUATION
%       θ02 –εθ CLOSURE BASED ON THE V2F
%       TURBULENCE MODEL"
%**************************************************************************
%
% Conventional models without radiation modifications:
%    t2-eq:  0 = 2 Pt - 2 et + ddy[(alpha+alphat/sigmat2) dt2dy] 
%    et-eq:  0 = Cp1 fp1 sqrt(e et/ (k t2)) Pt - Cd1 fd1 et^2 / t2 
%              - Cd2 fd2 e et / t2 + ddy[(alpha+alphat/sigmaet)detdy] 
%
% Models with radiative modifications:
%    t2-eq:  0 = 2 Pt - 2 et + ddy[(alpha+alphat/sigmat2) dt2dy] - 2kP Emt
%               +2kP Gt
%    et-eq:  0 = Cp1 fp1 sqrt(e et/ (k t2)) Pt - Cd1 fd1 et^2 / t2 
%              - Cd2 fd2 e et / t2 + ddy[(alpha+alphat/sigmaet)detdy] 
%              - 2*kP dEm dy dt dy + 2*kP dG dy dt dy
%
%    alphat = Cl v2 k^l/e t2^m/et
%    Pt     = - 2 alphat (dTdy)^2
%
% Radiation modification: introduction of an emission dissipation and
% absorption source in both t2 and et equations:
%
% alphat = Cl v2 k^l/e t2^m/et
%
% Modeling of radiation fluctuations:
%
% Em   = Cr1 t
% G    = Cr1 Cr2 t
%
% Where:
%
% Cr1  = (16/T0^4 T^3 + 48/T0^3 T^2 + 48/T0^2 T + 16/T0) / (ReT Pr Pl)
% Cr2  = kP/cr22 atan(cr22/kP)
% cr22 = 7.0;
%
% Radiative term modeling in t2 equation:
%
%  - 2 kP Cr1 t2 (1 - Cr2)
%
% Radiative term modeling in et equation:
%
%  - kP ( Cr1 et + 0.5 d Cr1 dy dt2 dy ) (1 - Cr2)       (if Cr2 != f(y))
%
%

function [ lam,t2,et,alphat ] = DWV( T,r,v2,t2,et,k,e,alpha,mu,kP,ReT,Pr,Pl,mesh,RadMod,kPMod)

    n = size(T,1);
    
    y        = mesh.y;
   	wallDist = min(y, 2-y);

    % Model constants
    Cl    = 0.28;
    Cp1   = 2.6;
    Cd1   = 2.0;
    Cd2   = 1.5;
    Ce2   = 1.9;
    siget = 1.0;
    sigt2 = 1.0;
    m     = 0.5;
    
    % Radiative model functions
    cr22 = 7.0;
    cr11 = 0.5;
    Cr1  = (16/1.5^4*T.^3 + 48/1.5^3*T.^2 + 48/1.5^2*T + 16/1.5)/(ReT*Pr*Pl); 
    Cr2  = kP./cr22.*atan(cr22./kP);
    
    % Relaxation factors
    underrelaxt2  = 0.8;
    underrelaxet  = 0.8;  
    
    % Time and length scales, eddy diffusivity and turbulent production
    Reps   = (mu.*e).^(1./4)./mu.*wallDist;
    Rturb  = (k.^2)./(mu.*e);
    
    % Model damping functions
    fd1    = ones(n,1);
    fd2    = sqrt(v2./k);
    
    % turbulent diffusivity and production
    if RadMod == 1
        R      = 0.5*(t2./(et + cr11*t2.*(Cr1).*(1-Cr2).*kP).*e./k); 
    else
        R      = 0.5*(t2./et.*e./k); 
    end
    fl     = 1;
    fl(1:n-1:n) = 0.0;
    
    alphat = max(r.*Cl.*v2.*k./e.*(2*R).^m,0.0) ;
    
      
    Pt  = alphat.*(mesh.ddy*T).^2;
    
    % ---------------------------------------------------------------------
    % et-equation
    %    0 = Cp1 fp1 sqrt(e et/ (k t2)) Pt - Cd1 fd1 et^2 / t2 
    %        - Cd2 fd2 e et / t2 + ddy[(alpha+alphat/sigmaet)detdy]    
    
    % effective diffusivity
    lam = alpha + alphat./siget;
    
    % diffusion matrix: lam*d2()/dy2 + dlam/dy d()/dy
    A =   bsxfun(@times, lam, mesh.d2dy2) ... 
        + bsxfun(@times, (mesh.ddy*lam), mesh.ddy);

    % implicitly treated source term : -Cd2 fd2 e/t2 - Cd1 fd1 et^2 / t2 
    for i=2:n-1
        A(i,i) = A(i,i) - Cd2*fd2(i)*e(i)/k(i) - Cd1*fd1(i)*et(i)/t2(i);
    end
        
    % Right-hand-side: - Cp1 fp1 sqrt(e et/ (k t2)) Pt
    b = - Cp1*sqrt(e(2:n-1).*et(2:n-1)./k(2:n-1)./t2(2:n-1)).*Pt(2:n-1);
    
    % Radiation implicit source term addition and RHS modification
    if(RadMod == 1)
        dCr1dy   = mesh.ddy  *Cr1;
        dt22dy2  = mesh.d2dy2*t2;
        for i=2:n-1
            A(i,i) = A(i,i) - 2*kP(i) * Cr1(i)*(1-Cr2(i)); %...
           %- 0.5*kP(i)*dCr1dy(i).*dt22dy2(i).*(1-Cr2(i))./et(i);
        end
    end
    
    % Boundary conditions
    et(1) = alpha(1)*(sqrt(t2(2)  )/(mesh.y(2)-mesh.y(1  ))).^2;
    et(n) = alpha(n)*(sqrt(t2(n-1))/(mesh.y(n)-mesh.y(n-1))).^2;
    
    %solve
    et = solveEq(et, A, b, underrelaxet);
    et(2:n-1) = max(et(2:n-1), 1.e-12);

    % ---------------------------------------------------------------------
	% t2-equation
    %    0 = 2 Pt - 2 et + ddy[(alpha+alphat/sigmat2) dt2dy]
    
    % effective diffusivity
    lam = alpha + alphat./sigt2;
        
    % diffusion matrix: lam*d2()/dy2 + dlam/dy d()/dy
    A =   bsxfun(@times, lam, mesh.d2dy2) ... 
        + bsxfun(@times, (mesh.ddy*lam), mesh.ddy);
    
    % implicitly treated source term
    for i=2:n-1
        A(i,i) = A(i,i) - 2*et(i)/t2(i);
    end
    
    % radiative implicit source modification
    if RadMod == 1
        for i=2:n-1
            A(i,i) = A(i,i) - 2*kP(i).*Cr1(i).*(1-Cr2(i));
        end
    end  
    
    % right-hand side
    b = - 2*Pt(2:n-1);
    
    % Boundary conditions
    t2(1) = 0;
    t2(n) = 0;
    
    % Solve
    t2 = solveEq(t2, A, b, underrelaxt2);
    
    lam = alpha + alphat;

end

