clear all
close all
clc
addpath('radiation');           % functions for the radiative calculations

radCase = 't100';
baseCase= 'b';
kPMod   =  0;

pathm    = strcat('solution/DNS/',radCase,'m');
patho    = strcat('solution/DNS/radiation/',radCase,'m');
Dm = importdata(pathm);
if (kPMod == 2) || (kPMod == 1)
    QR    = Dm(:,7);
    kP    = Dm(:,10);
    G     = Dm(:,11);
    if(isnan(G(10)))
        G = zeros(length(QR),1);
    end
    Em    = Dm(:,9 );
    ReT   = 3750;
    Pr    = 1.0;
    if strcmp(radCase,'H2OR')
        ReT = 16700;
        Pr  = 0.93;
    end
else
    QR    = Dm(:,7);
    Em    = Dm(:,8);
    G     = Dm(:,10);
    kP    = Dm(:,9);
    ReT   = 2900;
end

pathb    = strcat('~/Dropbox/RANSchannel_rad/solution/DNS/',baseCase,'m');
pathm    = strcat('~/Dropbox/RANSchannel_rad/solution/DNS/',radCase,'m');
patho    = strcat('~/Dropbox/RANSchannel_rad/solution/DNS/radiation/',radCase,'m');
Dm = importdata(pathm);
Db = importdata(pathb);

if strcmp(radCase,'H2OR') || strcmp(radCase,'Rb')
    Tnb = 500:20:1810;
else
    Tnb = 500:10:1010;
end
if kPMod == 2
    % Radiative model functions
    switch radCase
        case 'H2OR'
                        cr33 = 7*(ReT/2900).^(3./4).*Pr.^(1./2);
                        cr22 = 7*(ReT/2900).^(3./4).*Pr.^(1./2);
                        WVN   = ((cr33-cr22).*Dm(:,1).^2 - 2*(cr33-cr22).*Dm(:,1) +cr33);
                        WV = interp1(Dm(:,5).*1200+600,WVN,Tnb,'pchip');
%             load('lsR.mat');
%             WV = interp1(Db(:,5).*1200+600,ls,Tnb,'pchip');
        otherwise
                        cr33 = 7*(ReT/2900).^(3./4).*Pr.^(1./2);
                        cr22 = 7*(ReT/2900).^(3./4).*Pr.^(1./2);
                        WVN   = ((cr33-cr22).*Dm(:,1).^2 - 2*(cr33-cr22).*Dm(:,1) +cr33);
                        WV = interp1(Dm(:,5).*(955-573)+573,WVN,Tnb,'pchip');
%             load('lsr.mat');
%             WV = interp1(Db(:,5).*(955-573)+573,ls,Tnb,'pchip');
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
        fn       = kl./WV(i).*atan(WV(i)./(kl));
        fn2      = kl.*bl;
        ft(i)    = trapz(l*100,fn.*fn2);
        final(i) = ft(i)./(kP2(i)*Ib(i));
        i
    end
    switch radCase
        case 'H2OR'
            finInt = interp1(Tnb,final,Dm(:,5).*(1200)+600,'pchip');
        otherwise
            finInt = interp1(Tnb,final,Dm(:,5).*(955-573)+573,'pchip');
    end
else
    finInt = zeros(length(QR),1);
end
    

fid = fopen(patho,'w');
for i=1:length(QR)
    fprintf(fid,'%20.10f\t%20.10e\t%20.10e\t%20.10e\t%20.10e\t%20.10e\n',...
        Dm(i,1),kP(i),QR(i),Em(i),G(i),finInt(i));
end
fclose(fid);

