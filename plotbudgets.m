function plotbudgets(cas,gry,ReT,cha)


fig=figure('Position',[0 0 1500 1200]);
axp = axes('Position',[0 0 1 1],'Visible','off');
[ha,pos] = tight_subplot(4,length(cas),[0.001 0.001],[0.15 0.1],[0.07 0.07]);
TickLength = [0.02 0.05];
sk    = 3;
set(0,'defaultLineMarkerSize',5)
set(0,'defaultLineLineWidth',1.0)
hot  = [1 0 0];
cold = [0 0 1];
grey = [0.6 0.6 0.6];
grey2= [0.2 0.2 0.2];
enhanced = [0.5 0.5 0.5];
maxfig = zeros(3,length(cas));

for c=1:length(cas)
    
    f1 = c; f2 = c+length(cas); f3 = c+2*length(cas); f4 = c+3*length(cas);
    string = strcat('solution/DNS/',cas{c},'m');
    dns(c).m = importdata(string);
    string = strcat('solution/DNS/',cas{c},'f');
    dns(c).f = importdata(string);
    string = strcat('solution/DNS/',cas{c},'c');
    dns(c).c = importdata(string);
    
    axes(ha(f1)); hold on
    set(gca,'FontSize',14);
    set(gca,'XTickLabel','');
    xlabel('')
    set(gca,'XMinorTick','on','YMinorTick','on')
    set(gca,'ticklength',TickLength);
    if cha==0
        plot(dns(c).m(1:sk:end,1),dns(c).m(1:sk:end,5),'o','Color',grey);
        ylim([0 1.0]);
    elseif cha==1
        plot(dns(c).m(1:sk:end/2,1),1-dns(c).m(1:sk:end/2,5),'o','Color',hot);
        plot(2-dns(c).m(end/2+1:sk:end,1),dns(c).m(end/2+1:sk:end,5),'o','Color',cold);
        ylim([0 0.8]);
    elseif cha==2
        plot(dns(c).m(1:sk:end/2,1),1-dns(c).m(1:sk:end/2,5),'o','Color',hot);
        plot(2-dns(c).m(end/2+1:sk:end,1),dns(c).m(end/2+1:sk:end,5),'o','Color',cold);
        set(gca,'Xscale','log');
        ylim([0 0.8]);
    end
    yticks([0.2 0.4 0.6])
    if c==1
        ylabel('$\overline{T}$','interpreter','latex');
    else
        ylabel(''); set(gca,'YTickLabel','');
    end
    box on
    
    axes(ha(f2)); hold on
    set(gca,'FontSize',14);
    set(gca,'XTickLabel','');
    xlabel('')
    set(gca,'XMinorTick','on','YMinorTick','on')
    set(gca,'ticklength',TickLength);
    if cha==0
        plot(dns(c).f(1:sk:end,1),dns(c).f(1:sk:end,2),'o','Color',grey);
    elseif cha==1
        plot(dns(c).f(1:sk:end/2,1),dns(c).f(1:sk:end/2,2),'o','Color',hot);
        plot(2-dns(c).f(end/2+1:sk:end,1),dns(c).f(end/2+1:sk:end,2),'o','Color',cold);
    elseif cha==2
        plot(dns(c).f(1:sk:end/2,1),dns(c).f(1:sk:end/2,2),'o','Color',hot);
        plot(2-dns(c).f(end/2+1:sk:end,1),dns(c).f(end/2+1:sk:end,2),'o','Color',cold);
        set(gca,'Xscale','log');
    end
    if c==1
        ylabel('$\overline{\rho u^{\prime \prime} \theta^{\prime \prime}}$','interpreter','latex'); 
    else
        ylabel(''); set(gca,'YTickLabel','');
    end
    box on
    
    axes(ha(f3)); hold on
    set(gca,'FontSize',14);
    set(gca,'XMinorTick','on','YMinorTick','on')
    set(gca,'ticklength',TickLength);
    xlabel('$y/h$','Interpreter','latex');
    if cha==0
        plot(dns(c).c(1:sk:end,1),dns(c).c(1:sk:end,2),'o','Color',grey);
    elseif cha==1
        plot(dns(c).m(1:sk:end/2,1),dns(c).c(1:sk:end/2,2),'o','Color',hot);
        plot(2-dns(c).m(end/2+1:sk:end,1),dns(c).c(end/2+1:sk:end,2),'o','Color',cold);
    elseif cha==2
        plot(dns(c).m(1:sk:end/2,1),dns(c).c(1:sk:end/2,2),'o','Color',hot);
        plot(2-dns(c).m(end/2+1:sk:end,1),dns(c).c(end/2+1:sk:end,2),'o','Color',cold);
        set(gca,'Xscale','log');
    end
    if c==1
        ylabel('$\alpha_t$','interpreter','latex'); 
    else
        ylabel(''); set(gca,'YTickLabel','');
    end
    box on
    
    axes(ha(f4)); hold on
    set(gca,'FontSize',14);
    set(gca,'XMinorTick','on','YMinorTick','on')
    set(gca,'ticklength',TickLength);
    xlabel('$y/h$','Interpreter','latex');
    if cha==0
        plot(dns(c).c(1:sk:end,1),dns(c).c(1:sk:end,3),'o','Color',grey);
    elseif cha==1
        plot(dns(c).m(1:sk:end/2,1),dns(c).c(1:sk:end/2,3),'o','Color',hot);
        plot(2-dns(c).m(end/2+1:sk:end,1),dns(c).c(end/2+1:sk:end,3),'o','Color',cold);
    elseif cha==2
        plot(dns(c).m(1:sk:end/2,1),dns(c).c(1:sk:end/2,3),'o','Color',hot);
        plot(2-dns(c).m(end/2+1:sk:end,1),dns(c).c(end/2+1:sk:end,3),'o','Color',cold);
        set(gca,'Xscale','log');
    end
    if c==1
        ylabel('$\alpha_t$','interpreter','latex'); 
    else
        ylabel(''); set(gca,'YTickLabel','');
    end
    box on
    
    % y,u,T,diffT,mut,alphat,THF,k,e,v2,t2,et);
    
    string = strcat('solution/SA/NO/',cas{c});
    sa.no  = importdata(string);
    string = strcat('solution/V2F/NO/',cas{c});
    v2f.no = importdata(string);
    string = strcat('solution/V2F/DWX/',cas{c});
    v2f.dwx = importdata(string);
    string = strcat('solution/V2F/DWX/',cas{c},'_rad_0');
    v2f.dwr = importdata(string);
    y = sa.no(:,1);
    
    if(gry==1)
        string = strcat('solution/V2F/DWX/',cas{c},'_rad_2');
        v2f.dwrng = importdata(string);
    end
    
    axes(ha(f1)); hold on
    if cha==0
        plot(y(1:end),sa.no(1:end,3),':','Color',grey2);
        plot(y(1:end),v2f.no(1:end,3),'--','Color',grey2);
        plot(y(1:end),v2f.dwx(1:end,3),'-.','Color',grey2);
        plot(y(1:end),v2f.dwr(1:end,3),'-','Color',grey2);
        if(gry==1)
            plot(y(1:end),v2f.dwrng(1:end,3),'-','Color',enhanced);
        end
    else
        plot(y(1:end/2),1-sa.no(1:end/2,3),':','Color',hot);
        plot(2-y(end/2+1:end),sa.no(end/2+1:end,3),':','Color',cold);
        plot(y(1:end/2),1-v2f.no(1:end/2,3),'--','Color',hot);
        plot(2-y(end/2+1:end),v2f.no(end/2+1:end,3),'--','Color',cold);
        plot(y(1:end/2),1-v2f.dwx(1:end/2,3),'-.','Color',hot);
        plot(2-y(end/2+1:end),v2f.dwx(end/2+1:end,3),'-.','Color',cold);
        plot(y(1:end/2),1-v2f.dwr(1:end/2,3),'-','Color',hot);
        plot(2-y(end/2+1:end),v2f.dwr(end/2+1:end,3),'-','Color',cold);
        if(gry==1)
            plot(y(1:end/2),1-v2f.dwrng(1:end/2,3),'-','Color',enhanced);
            plot(2-y(end/2+1:end),v2f.dwrng(end/2+1:end,3),'-','Color',enhanced);
        end
    end

    axes(ha(f2)); hold on
    if cha==0
        plot(y(1:end),sa.no(1:end,7),':','Color',grey2);
        plot(y(1:end),v2f.no(1:end,7),'--','Color',grey2);
        plot(y(1:end),v2f.dwx(1:end,7),'-.','Color',grey2);
        plot(y(1:end),v2f.dwr(1:end,7),'-','Color',grey2);
        if(gry==1)
            plot(y(1:end),v2f.dwrng(1:end,7),'-','Color',enhanced);
        end
    else
        plot(y(1:end/2),sa.no(1:end/2,7),':','Color',hot);
        plot(2-y(end/2+1:end),sa.no(end/2+1:end,7),':','Color',cold);
        plot(y(1:end/2),v2f.no(1:end/2,7),'--','Color',hot);
        plot(2-y(end/2+1:end),v2f.no(end/2+1:end,7),'--','Color',cold);
        plot(y(1:end/2),v2f.dwx(1:end/2,7),'-.','Color',hot);
        plot(2-y(end/2+1:end),v2f.dwx(end/2+1:end,7),'-.','Color',cold);
        plot(y(1:end/2),v2f.dwr(1:end/2,7),'-','Color',hot);
        plot(2-y(end/2+1:end),v2f.dwr(end/2+1:end,7),'-','Color',cold);
        if(gry==1)
            plot(y(1:end/2),v2f.dwrng(1:end/2,7),'-','Color',enhanced);
            plot(2-y(end/2+1:end),v2f.dwrng(end/2+1:end,7),'-','Color',enhanced);
        end
    end
    maxfig(2,c) = max([sa.no(:,7);sa.no(:,7);v2f.no(:,7);v2f.dwx(:,7);v2f.dwr(:,7)]);
    
    axes(ha(f3)); hold on
    if cha==0
        plot(y(1:end),sa.no(1:end,11),':','Color',grey2);
        plot(y(1:end),v2f.no(1:end,11),'--','Color',grey2);
        plot(y(1:end),v2f.dwx(1:end,11),'-.','Color',grey2);
        plot(y(1:end),v2f.dwr(1:end,11),'-','Color',grey2);
        if(gry==1)
            plot(y(1:end),v2f.dwrng(1:end,11),'-','Color',enhanced);
        end
    else
        plot(y(1:end/2),sa.no(1:end/2,11),':','Color',hot);
        plot(2-y(end/2+1:end),sa.no(end/2+1:end,11),':','Color',cold);
        plot(y(1:end/2),v2f.no(1:end/2,11),'--','Color',hot);
        plot(2-y(end/2+1:end),v2f.no(end/2+1:end,11),'--','Color',cold);
        plot(y(1:end/2),v2f.dwx(1:end/2,11),'-.','Color',hot);
        plot(2-y(end/2+1:end),v2f.dwx(end/2+1:end,11),'-.','Color',cold);
        plot(y(1:end/2),v2f.dwr(1:end/2,11),'-','Color',hot);
        plot(2-y(end/2+1:end),v2f.dwr(end/2+1:end,11),'-','Color',cold);
        if(gry==1)
            plot(y(1:end/2),v2f.dwrng(1:end/2,11),'-','Color',enhanced);
            plot(2-y(end/2+1:end),v2f.dwrng(end/2+1:end,11),'-','Color',enhanced);
        end
    end
    maxfig(3,c) = max([sa.no(:,11);sa.no(:,11);v2f.no(:,11);v2f.dwx(:,11);v2f.dwr(:,11)]);

    axes(ha(f4)); hold on
    if cha==0
        plot(y(1:end),-sa.no(1:end,12),':','Color',grey2);
        plot(y(1:end),-v2f.no(1:end,12),'--','Color',grey2);
        plot(y(1:end),-v2f.dwx(1:end,12),'-.','Color',grey2);
        plot(y(1:end),-v2f.dwr(1:end,12),'-','Color',grey2);
        if(gry==1)
            plot(y(1:end),-v2f.dwrng(1:end,12),'-','Color',enhanced);
        end
    else
        plot(y(1:end/2),-sa.no(1:end/2,12),':','Color',hot);
        plot(2-y(end/2+1:end),-sa.no(end/2+1:end,12),':','Color',cold);
        plot(y(1:end/2),-v2f.no(1:end/2,12),'--','Color',hot);
        plot(2-y(end/2+1:end),-v2f.no(end/2+1:end,12),'--','Color',cold);
        plot(y(1:end/2),-v2f.dwx(1:end/2,12),'-.','Color',hot);
        plot(2-y(end/2+1:end),-v2f.dwx(end/2+1:end,12),'-.','Color',cold);
        plot(y(1:end/2),-v2f.dwr(1:end/2,12),'-','Color',hot);
        plot(2-y(end/2+1:end),-v2f.dwr(end/2+1:end,12),'-','Color',cold);
        if(gry==1)
            plot(y(1:end/2),v2f.dwrng(1:end/2,11),'-','Color',enhanced);
            plot(2-y(end/2+1:end),v2f.dwrng(end/2+1:end,11),'-','Color',enhanced);
        end
    end
    maxfig(4,c) = min([-sa.no(:,12);-sa.no(:,12);-v2f.no(:,12);-v2f.dwx(:,12);-v2f.dwr(:,12)]);

end

for i=2:4
    maxf = max(maxfig(i,:));
    maxf = max(0,maxf);
    minf = min(maxfig(i,:));
    minf = min(0,minf);
    for c=1:length(cas)
        f = c+(i-1)*length(cas);
        axes(ha(f)); hold on
        ylim([minf*1.1 maxf*1.1])
    end
end



end