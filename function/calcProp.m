%**************************************************************************
%       function, calculating transport properties 
%**************************************************************************
% Inputs:
%   T           temperature
%   ReT         friction Reynolds number ReT=utau r_wall h/ mu_wall
%   casename    fluid
%
% Output:
%   r           density
%   mu          molecular viscosity
%
function [r,mu,alpha] = calcProp(T,dns,ReT,Pr,casename,T0,y)

    n = size(T,1);
    
    if strcmp(casename,'cRets')
        r  = 1./T;
        mu = (T.^-0.5)./ReT;
        alpha = ones(n,1)./ReT./Pr;
    elseif strcmp(casename,'liquidLike')
        r  = ones(n,1);
        mu = (T.^-1.0)./ReT;
        alpha = ones(n,1)./ReT./Pr;
    elseif strcmp(casename,'gasLike')
        r  = 1./T;
        mu = (T.^0.7)./ReT;
        alpha = ones(n,1)./ReT./Pr;
    elseif strcmp(casename,'vardens')
        r     = T0./(T+T0);
        mu    = ones(n,1)./ReT;
        alpha = ones(n,1)./ReT./Pr;
    elseif strcmp(casename,'varvisc')
        r     = T0./(T+T0);
        mu    = ((T+T0)./(T0)).^1.15./ReT;
        alpha = ((T+T0)./(T0)).^1.35./ReT./Pr;
    elseif strcmp(casename,'nodensity')
        r     = interp1(dns(:,1),dns(:,14),y,'pchip');
        mu    = ones(n,1)./ReT;
        alpha = ones(n,1)./ReT./Pr;
    else
        r     = ones(n,1);
        mu    = ones(n,1)./ReT;
        alpha = ones(n,1)./ReT./Pr;
    end

end