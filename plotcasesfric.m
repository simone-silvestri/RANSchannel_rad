function plotcasesfric(cas,gry,ReT,var)


fig=figure('Position',[0 0 1200 1200]);
axp = axes('Position',[0 0 1 1],'Visible','off');
[ha,pos] = tight_subplot(var,length(cas),[0.001 0.001],[0.15 0.1],[0.07 0.07]);
TickLength = [0.02 0.05];
sk    = 4;
set(0,'defaultLineMarkerSize',6)
set(0,'defaultLineLineWidth',1.5)
hot  = [1 0 0];
cold = [0 0 1];
grey = [0.6 0.6 0.6];
grey2= [0.2 0.2 0.2];
DNS  = [0.6 0.6 0.6];
v2fno= [0.4 0.4 0.4];
v2fdx= [1.0 0.0 0.0];
v2fdr= [0.0 0.0 1.0];
enhanced = [0. 0. 0.];
maxfig = zeros(3,length(cas));
width  = 0.8;

for c=1:length(cas)
    
    if strcmp(cas{c},'H2OR') || strcmp(cas{c},'Rb')
        ReT = 16700;
        sk  = 8;
    end
    f1 = c; f2 = c+length(cas); f3 = c+2*length(cas);
    string = strcat('solution/DNS/',cas{c},'m');
    dns(c).m = importdata(string);
    string = strcat('solution/DNS/',cas{c},'f');
    dns(c).f = importdata(string);
    dns(c).a = dns(c).f(:,2)./(-dns(c).f(:,3))/ReT;
    
    axes(ha(f1)); hold on
    set(gca,'FontSize',14);
    set(gca,'XTickLabel','');
    xlabel('')
    set(gca,'XMinorTick','on','YMinorTick','on')
    set(gca,'ticklength',TickLength);
    
    
    fr1 = (1-dns(c).m(1:sk:end/2,5)/
    plot(dns(c).m(1:sk:end/2,1),,'o','Color',hot,'LineWidth',width);
    plot(2-dns(c).m(end/2+1:sk:end,1),dns(c).m(end/2+1:sk:end,5),'o','Color',cold,'LineWidth',width);
    set(gca,'Xscale','log');
    ylim([0 0.8]);
    yticks([0.2 0.4 0.6])
    if c==1
        ylabel('$\overline{T}$','interpreter','latex');
    else
        ylabel(''); set(gca,'YTickLabel','');
    end
    box on
    
    if (var>=2)
        axes(ha(f2)); hold on
        set(gca,'FontSize',14);
        set(gca,'XTickLabel','');
        xlabel('')
        set(gca,'XMinorTick','on','YMinorTick','on')
        set(gca,'ticklength',TickLength);
        if cha==0 || cha==3
            plot(dns(c).f(1:sk:end,1),dns(c).f(1:sk:end,2),'o','Color',DNS,'LineWidth',width);
            xticks([0.5 1 1.5])
        elseif cha==1
            plot(dns(c).f(1:sk:end/2,1),dns(c).f(1:sk:end/2,2),'o','Color',hot,'LineWidth',width);
            plot(2-dns(c).f(end/2+1:sk:end,1),dns(c).f(end/2+1:sk:end,2),'o','Color',cold,'LineWidth',width);
        elseif cha==2
            plot(dns(c).f(1:sk:end/2,1),dns(c).f(1:sk:end/2,2),'o','Color',hot,'LineWidth',width);
            plot(2-dns(c).f(end/2+1:sk:end,1),dns(c).f(end/2+1:sk:end,2),'o','Color',cold,'LineWidth',width);
            set(gca,'Xscale','log');
        end
        if c==1
            ylabel('$\overline{\rho u^{\prime \prime} \theta^{\prime \prime}}$','interpreter','latex');
        else
            ylabel(''); set(gca,'YTickLabel','');
        end
        box on
    end
    
    if(var>=3)
        axes(ha(f3)); hold on
        set(gca,'FontSize',14);
        set(gca,'XMinorTick','on','YMinorTick','on')
        set(gca,'ticklength',TickLength);
        xlabel('$y/h$','Interpreter','latex');
        if cha==0 || cha==3
            plot(dns(c).m(1:sk:end,1),dns(c).a(1:sk:end),'o','Color',DNS,'LineWidth',width);        
            xticks([0.5 1 1.5])
        elseif cha==1
            plot(dns(c).m(1:sk:end/2,1),dns(c).a(1:sk:end/2),'o','Color',hot,'LineWidth',width);
            plot(2-dns(c).m(end/2+1:sk:end,1),dns(c).a(end/2+1:sk:end),'o','Color',cold,'LineWidth',width);
        elseif cha==2
            plot(dns(c).m(1:sk:end/2,1),dns(c).a(1:sk:end/2),'o','Color',hot,'LineWidth',width);
            plot(2-dns(c).m(end/2+1:sk:end,1),dns(c).a(end/2+1:sk:end),'o','Color',cold,'LineWidth',width);
            set(gca,'Xscale','log');
        end
        if c==1
            ylabel('$\alpha_t$','interpreter','latex');
        else
            ylabel(''); set(gca,'YTickLabel','');
        end
        box on
    end
    % y,u,T,diffT,mut,alphat,THF,k,e,v2,t2,et);
    
    string = strcat('solution/SA/NO/',cas{c});
    sa.no  = importdata(string);
    string = strcat('solution/V2F/NO/',cas{c});
    v2f.no = importdata(string);
    string = strcat('solution/V2F/DWX/',cas{c});
    v2f.dwx = importdata(string);
    string = strcat('solution/V2F/DWX/',cas{c},'_rad_0');
    v2f.dwr = importdata(string);
    string = strcat('solution/V2F/V2T/',cas{c});
    v2f.v2t = importdata(string);
    string = strcat('solution/V2F/V2T/',cas{c},'_rad_0');
    v2f.v2r = importdata(string);
    y = sa.no(:,1);
    
    if(gry==1)
        string = strcat('solution/V2F/DWX/',cas{c},'_rad_2');
        v2f.dwrng = importdata(string);
    end
    
    axes(ha(f1)); hold on
    if cha==0
        %plot(y(1:end),sa.no(1:end,3),':','Color',grey2);
        plot(y(1:end),v2f.no(1:end,3),'--','Color',v2fno);
        plot(y(1:end),v2f.dwx(1:end,3),'-.','Color',v2fdx);
        plot(y(1:end),v2f.dwr(1:end,3),'-','Color',v2fdr);
        if(gry==1)
            plot(y(1:end),v2f.dwrng(1:end,3),':','Color',enhanced);
        end
    elseif cha == 3
        %plot(y(1:end),sa.no(1:end,3),':','Color',grey2);
        plot(y(1:end),v2f.no(1:end,3),'--','Color',v2fno);
        plot(y(1:end),v2f.v2t(1:end,3),'-.','Color',v2fdx);
        plot(y(1:end),v2f.v2r(1:end,3),'-','Color',v2fdr);
        if(gry==1)
            plot(y(1:end),v2f.dwrng(1:end,3),':','Color',enhanced);
        end
    else
        %plot(y(1:end/2),1-sa.no(1:end/2,3),':','Color',hot);
        %plot(2-y(end/2+1:end),sa.no(end/2+1:end,3),':','Color',cold);
        plot(y(1:end/2),1-v2f.no(1:end/2,3),'--','Color',hot);
        plot(2-y(end/2+1:end),v2f.no(end/2+1:end,3),'--','Color',cold);
        plot(y(1:end/2),1-v2f.dwx(1:end/2,3),'-.','Color',hot);
        plot(2-y(end/2+1:end),v2f.dwx(end/2+1:end,3),'-.','Color',cold);
        plot(y(1:end/2),1-v2f.dwr(1:end/2,3),'-','Color',hot);
        plot(2-y(end/2+1:end),v2f.dwr(end/2+1:end,3),'-','Color',cold);
        if(gry==1)
            plot(y(1:end/2),1-v2f.dwrng(1:end/2,3),':','Color',enhanced);
            plot(2-y(end/2+1:end),v2f.dwrng(end/2+1:end,3),':','Color',enhanced);
        end
    end
    
    if(var>=2)
        axes(ha(f2)); hold on
        if cha==0
            %plot(y(1:end),sa.no(1:end,7),':','Color',grey2);
            plot(y(1:end),v2f.no(1:end,7),'--','Color',v2fno);
            plot(y(1:end),v2f.dwx(1:end,7),'-.','Color',v2fdx);
            plot(y(1:end),v2f.dwr(1:end,7),'-','Color',v2fdr);
            if(gry==1)
                plot(y(1:end),v2f.dwrng(1:end,7),':','Color',enhanced);
            end
        elseif cha==3
            %plot(y(1:end),sa.no(1:end,7),':','Color',grey2);
            plot(y(1:end),v2f.no(1:end,7),'--','Color',v2fno);
            plot(y(1:end),v2f.v2t(1:end,7),'-.','Color',v2fdx);
            plot(y(1:end),v2f.v2r(1:end,7),'-','Color',v2fdr);
            if(gry==1)
                plot(y(1:end),v2f.dwrng(1:end,7),':','Color',enhanced);
            end 
        else
            %plot(y(1:end/2),sa.no(1:end/2,7),':','Color',hot);
            %plot(2-y(end/2+1:end),sa.no(end/2+1:end,7),':','Color',cold);
            plot(y(1:end/2),v2f.no(1:end/2,7),'--','Color',hot);
            plot(2-y(end/2+1:end),v2f.no(end/2+1:end,7),'--','Color',cold);
            plot(y(1:end/2),v2f.dwx(1:end/2,7),'-.','Color',hot);
            plot(2-y(end/2+1:end),v2f.dwx(end/2+1:end,7),'-.','Color',cold);
            plot(y(1:end/2),v2f.dwr(1:end/2,7),'-','Color',hot);
            plot(2-y(end/2+1:end),v2f.dwr(end/2+1:end,7),'-','Color',cold);
            if(gry==1)
                plot(y(1:end/2),v2f.dwrng(1:end/2,7),':','Color',enhanced);
                plot(2-y(end/2+1:end),v2f.dwrng(end/2+1:end,7),':','Color',enhanced);
            end
        end
        maxfig(2,c) = max([sa.no(:,7);sa.no(:,7);v2f.no(:,7);v2f.dwx(:,7);v2f.dwr(:,7)]);
    end
    
    if(var>=3)
        axes(ha(f3)); hold on
        if cha==0
            %plot(y(1:end),sa.no(1:end,6).*sa.no(1:end,13),':','Color',grey2);
            plot(y(1:end),v2f.no(1:end,6).*v2f.no(1:end,13),'--','Color',v2fno);
            plot(y(1:end),v2f.dwx(1:end,6).*v2f.dwx(1:end,13),'-.','Color',v2fdx);
            plot(y(1:end),v2f.dwr(1:end,6).*v2f.dwr(1:end,13),'-','Color',v2fdr);
            if(gry==1)
                plot(y(1:end),v2f.dwrng(1:end,6).*v2f.dwrng(1:end,13),':','Color',enhanced);
            end
        elseif cha==3
            %plot(y(1:end),sa.no(1:end,6).*sa.no(1:end,13),':','Color',grey2);
            plot(y(1:end),v2f.no(1:end,6).*v2f.no(1:end,13),'--','Color',v2fno);
            plot(y(1:end),v2f.v2t(1:end,6).*v2f.v2t(1:end,13),'-.','Color',v2fdx);
            plot(y(1:end),v2f.v2t(1:end,6).*v2f.v2t(1:end,13),'-','Color',v2fdr);
            if(gry==1)
                plot(y(1:end),v2f.dwrng(1:end,6).*v2f.dwrng(1:end,13),':','Color',enhanced);
            end
        else
            %plot(y(1:end/2),sa.no(1:end/2,6).*sa.no(1:end/2,13),':','Color',hot);
            %plot(2-y(end/2+1:end),sa.no(end/2+1:end,6).*sa.no(end/2+1:end,13),':','Color',cold);
            plot(y(1:end/2),v2f.no(1:end/2,6).*v2f.no(1:end/2,13),'--','Color',hot);
            plot(2-y(end/2+1:end),v2f.no(end/2+1:end,6).*v2f.no(end/2+1:end,13),'--','Color',cold);
            plot(y(1:end/2),v2f.dwx(1:end/2,6).*v2f.dwx(1:end/2,13),'-.','Color',hot);
            plot(2-y(end/2+1:end),v2f.dwx(end/2+1:end,6).*v2f.dwx(end/2+1:end,13),'-.','Color',cold);
            plot(y(1:end/2),v2f.dwr(1:end/2,6).*v2f.dwr(1:end/2,13),'-','Color',hot);
            plot(2-y(end/2+1:end),v2f.dwr(end/2+1:end,6).*v2f.dwr(end/2+1:end,13),'-','Color',cold);
            if(gry==1)
                plot(y(1:end/2),v2f.dwrng(1:end/2,6).*v2f.dwrng(1:end/2,13),':','Color',enhanced);
                plot(2-y(end/2+1:end),v2f.dwrng(end/2+1:end,6).*v2f.dwrng(end/2+1:end,13),':','Color',enhanced);
            end
        end
        maxfig(3,c) = max([sa.no(:,6);sa.no(:,6);v2f.no(:,6);v2f.dwx(:,6);v2f.dwr(:,6)]);
    end
end

for i=2:var
    maxf = max(maxfig(i,:));
    for c=1:length(cas)
        f = c+(i-1)*length(cas);
        axes(ha(f)); hold on
        ylim([0 maxf*1.1])
    end
end



end