function [ lam, Prt ] = PRT( mu,mut,Prt,Pr )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    C = 3.0; 
    PeT = mut./mu.*Pr;
    gam = 1.0./(Prt + 0.1*Pr.^0.83);
    A = sqrt(2*(1./Prt-gam));
    Prt = 1./(gam + C*PeT.*A - ((C*PeT).^2).*(1-exp(-A./(C*PeT))));
    
    lam = mu./Pr + mut./Prt;

end

