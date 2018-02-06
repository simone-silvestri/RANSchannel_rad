function [ lam, alphat ] = PRT( mu,mut,alpha,T,r,qy,ReT,MESH,RadMod)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    % model constants
    C1 = 0.5882; 
    C2 = 0.228;
    C3 = 0.0441;
    C4 = 5.165;
    
    gam = mut./mu;
    
    Prt0 = 1./(C1 + C2.*gam - C3.*gam.^2.*(1-exp(-C4./gam)));
    
    
    if RadMod == 2
        Prt = r./ReT.*abs(MESH.ddy*T)./abs(qy).*(Prt0+mut.*ReT);
    else
        Prt = Prt0;
    end
    
    lam    = alpha + mut./Prt;
    alphat = mut./Prt;
    
    
end

