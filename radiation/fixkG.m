function kG = fixkG(T,radCase,ls,kP)


if strcmp(radCase,'H2OR') || strcmp(radCase,'Rb')
    Tnb = 500:20:1810;
else
    Tnb = 500:20:1010;
end

% Radiative model functions
switch radCase
    case 'H2OR'
        WV = interp1(T.*1200+600,ls,Tnb,'pchip');
    otherwise
        WV = interp1(T.*(955-573)+573,ls,Tnb,'pchip');
end

for i=1:length(Tnb)
    tmp      = sprintf('%d\n',Tnb(i));
    switch radCase
        case 'CO2'
            tmp      = strcat('radiation/hitran/spec',radCase,'/',tmp,'K.txt');
            spec     = importdata(tmp);
            kl       = spec(:,2)*10;
        case 'H2OP'
            tmp      = strcat('radiation/hitran/specH2O/',tmp,'K.txt');
            spec     = importdata(tmp);
            kl       = spec(:,2)*10+2;
        case 'H2OR'
            tmp      = strcat('radiation/hitran/specH2O/',tmp,'K.txt');
            spec     = importdata(tmp);
            kl       = spec(:,2)*100;
        otherwise
            tmp      = strcat('radiation/hitran/spec',radCase,'/',tmp,'K.txt');
            spec     = importdata(tmp);
            kl       = spec(:,2)*100;
    end
    l        = spec(:,1);
    bl       = blackbody(Tnb(i),l);
    Ib(i)    = trapz(l*100,bl);
    kP2(i)   = trapz(l*100,kl.*bl)./Ib(i);
    fn       = kl.*WV(i).*atan(1./(kl.*WV(i)));
    fn2      = kl.*bl;
    ft(i)    = trapz(l*100,fn.*fn2);
    final(i) = ft(i)./(kP2(i)*Ib(i));
    if mod(i,10)==0
        disp(i);
    end
    fun = @(x) x.*WV(i).*atan(1./(x.*WV(i)))-final(i);
    kGt(i) = fzero(fun,[0 10000]);
end
switch radCase
    case 'H2OR'
%         finInt = interp1(Tnb,final,T.*(1200)+600,'pchip');
          kG = interp1(Tnb,kGt,T.*(1200)+600,'pchip');
    otherwise
        kG = interp1(Tnb,kGt,T.*(955-573)+573,'pchip');
end
% switch radCase
%     case 'H2OR'
%         ls = interp1(Tnb,WV,T.*1200+600,'pchip');
%     otherwise
%         ls = interp1(Tnb,WV,T.*(955-573)+573,'pchip');
% end
% for i=1:length(T)
%     if ls(i)<=1e-4
%         ls(i) = 1e-4;
%     end
%     fun = @(x) x.*ls(i).*atan(1./(x.*ls(i)))-finInt(i);
%     kG(i) = fzero(fun,1.5);
% end
%
end