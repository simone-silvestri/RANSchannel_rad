function [ k,e,u,mut,r ] = DNS( Dk,Dm,mu,ReT,MESH,vardens )
%DNS Summary of this function goes here
%   Detailed explanation goes here

    k    = interp1(Dk(:,1),Dk(:,2),MESH.y,'spline');
    e    = interp1(Dk(:,1),-Dk(:,3),MESH.y,'spline');
    u    = interp1(Dm(:,1),Dm(:,4),MESH.y,'spline');
    if vardens==1 
        r    = interp1(Dm(:,1),Dm(:,14),MESH.y,'spline');
    else
        r = ones(length(MESH.y),1);
    end
    u    = abs(u); 
    
    n = length(u);
       
    tv       = (MESH.ddy*u)./ReT;
    ut       = sqrt(tv(1)./r(1)/2-tv(end)./r(end)/2);
    Retau    = ReT.*ut.*r;
    y        = MESH.y;
   	wallDist = min(y, 2-y);
    yplus    = wallDist.*Retau; 
    % Model constants
    cmu  = 0.09; 
    
    % ---------------------------------------------------------------------
    % eddy viscosity
    ReTurb = r.*(k.^2)./(mu.*e);
    fmue   = (1-exp(-yplus/70)).*(1.0+3.45./(ReTurb.^0.5));   
    fmue(1:n-1:n) = 0.0;
    
    mut  = cmu*fmue./e.*k.^2;

end

