%**************************************************************************
%       function, transform of wall units uplus to van Driest uvd
%**************************************************************************
%   duvd = sqrt(rho/rho_wall) duplus
%
% Inputs:
%   uplus   velocity with transformation using wall values (utau)
%   r       density
%
% Output:
%   uvd     velocity with van Driest transformation
%

function [uvs] = velTransVS(uplus,m)

    n = size(uplus,1);
    uvs = zeros(n,1);
   
    mw = m(1);

    uvs(1) = 0.0;

    for i=2:n
        uvs(i) = uvs(i-1) + 0.5*(m(i)+m(i-1))*(uplus(i)-uplus(i-1));
    end
end
