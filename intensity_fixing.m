%% LED parlaklık optimizasyonu
clear
link_m10='m10\';
link_m9='m9\';
ngroup=25;LED={'UV';'Blue';'White';'Red';'Green'};LED=sort(LED); nLED=length(LED);
LEDtr={'Mavi';'Yeşil';'Kırmızı';'Morötesi';'Beyaz';};
num=1:5; nnum=length(num);
matfiles_new=dir(fullfile(link_m10, '*.txt'));
matfiles_old=dir(fullfile(link_m9, '*.txt'));
nfiles=length(matfiles_new);

Res_new=[148.4 170.6 191.6 211.3 238.2;...
    595.1 653.7 708 794 873;...
    27.9 36.4 43.8 54.8 66.8;...
    486.2 526.6 559.2 604.5 643;...
    315.9 338.4 377.3 419.2 467.8;];
Res_old=[518.3,543.2,568.5,604.7,643.6;...
    622.5,670,734,822,904;...
    130,139.6,155.8,171.5,184.7;...
    599.3,645.7,689,746,787;...
    613.6,637.4,673,714,755];
Res_old(1,:)=Res_old(1,:)-300;
Res_old(5,:)=Res_old(5,:)-300;

for i=1:ngroup
    for j=1:nnum
        for k=1:nLED
            names=matfiles_new(i+ngroup*(j-1)+nnum*ngroup*(k-1)).name;
            D=load([link_m10 names]);
            names=matfiles_old(i+ngroup*(j-1)+nnum*ngroup*(k-1)).name;
            D2=load([link_m9 names]);
            X(:,1)=D(:,1); Y_new(:,i,j,k)=D(:,2)/100; Y_old(:,i,j,k)=D2(:,2)/100;
        end
    end
end
N=length(X);
for i=1:N
    for j=1:nnum 
        for k=1:nLED
            Y_avgnew(i,j+nnum*(k-1))=mean(Y_new(i,:,j,k));
            Y_avgold(i,j+nnum*(k-1))=mean(Y_old(i,:,j,k));
        end
    end
end

for i=1:nnum*nLED
    Base1(i,1)=mean(Y_avgnew(100:350,i));
    Base2(i,1)=mean(Y_avgold(100:350,i));
end
for i=1:N
    for j=1:nnum*nLED
        if Y_avgnew(i,j)>=Base1(j)
            yavg_new(i,j)=Y_avgnew(i,j)-Base1(j);
        else
            yavg_new(i,j)=0;
        end
        if Y_avgnew(i,j)>=Base2(j)
            yavg_old(i,j)=Y_avgold(i,j)-Base2(j);
        else
            yavg_old(i,j)=0;
        end
    end
end

for i=1:nLED
    for j=1:nnum
        int_new(i,j)=sum(yavg_new(:,j+nnum*(i-1)))*(X(end)-X(1))/N;
        int_old(i,j)=sum(yavg_old(:,j+nnum*(i-1)))*(X(end)-X(1))/N;
    end
end
color={'+b';'+g';'+r';'+m';'+c'};
color2={'b';'g';'r';'m';'c'}; v2=[4 2 4 3 2];
%% Fit yapma
clc
for i=1:nLED
    color_o{i}=join(['İstasyon 1 ' LEDtr{i}]);
    color_n{i}=join(['İstasyon 2 ' LEDtr{i}]);
    f_new=fit(Res_new(i,:)',int_new(i,:)','poly1');
    ab_n(i,:)=coeffvalues(f_new);
    Rxn(:,i)=linspace(Res_new(i,1)*.9,Res_new(i,end)*1.1,N);
    intfn(:,i)=ab_n(i,1).*Rxn(:,i)+ab_n(i,2);
    
    
    f_old=fit(Res_old(i,:)',int_old(i,:)','poly1');
    ab_o(i,:)=coeffvalues(f_old);
    Rxo(:,i)=linspace(Res_old(i,1)*.9,Res_old(i,end)*1.1,N);
    intfo(:,i)=ab_o(i,1).*Rxo(:,i)+ab_o(i,2);
    
    hold on
    pn(i)=plot(Rxn(:,i),intfn(:,i),join([char(color2{i})],''),'linewidth',2);
end
for i=1:nLED
    po(i)=plot(Rxo(:,i),intfo(:,i),join(['--',char(color2{i})],''),'linewidth',2);
end
intt=7900;
p=plot(linspace(20,900,N),ones(N,1)*intt,'k','linewidth',2);

for i=1:nLED
    pasd=plot(Res_new(i,:),int_new(i,:),'--+k');
    plot(Res_old(i,:),int_old(i,:),'--+k')
end
xlim([0 1000]),ylim([6000 10000])
xlabel('Akım Sınırlayıcı Direnç (\Omega)'),ylabel('LED Spektrum İntegrali (\gamma /ms)')
l1=legend([p po pasd pn],[{'Hedef Parlaklık'} color_o {'Ölçüm değerleri'} color_n],'Location','n');grid on
l1.NumColumns=2; legend boxoff
l1.Position=l1.Position.*[1 1.05 1 1];
set(gca, 'YMinorTick','on', 'XMinorTick','on')
set(gca, 'YMinorGrid','on', 'XMinorGrid','on')
set(gcf,'position',[175,345,870,420])
exportgraphics(gcf,'IvsR.png','Resolution',600)
for i=1:nLED
    Rneedn(i,1)=(intt-ab_n(i,2))/ab_n(i,1);
    Rneedo(i,1)=(intt-ab_o(i,2))/ab_o(i,1);
end
T=table(char(LEDtr),num2str(round(ab_o,2)),num2str(round(Rneedo,2)),num2str(round(ab_n,2)),num2str(round(Rneedn)), 'VariableNames', {'LED Çeşidi',...
    'istasyon 1: a ve b','istasyon 1: direnç','istasyon 2: a ve b','istasyon 2: direnç'});

disp 'Resistant-Intensity equation [int=a*R + b]'
disp '=========================================='
disp 'Table of Old and New LED Stations'
disp(T)
%% R tasarımı
Rsetn=[150 20 0;560 100 3.9;39 2.2 0;390 120 0;220 68 68];
Rseto=[200 33 0;470 180 47;100 47 0;470 150 0;220 120 3.3];
Rn=sum(Rsetn')';
Ro=sum(Rseto')';
for i=1:nLED
    yeni(i,1)=ab_n(i,1)*Rn(i)+ab_n(i,2);
    eski(i,1)=ab_o(i,1)*Ro(i)+ab_o(i,2);
    ery(i,1)=abs(yeni(i)-intt)/intt*100;
    ere(i,1)=abs(eski(i)-intt)/intt*100;
end
disp(['Wanted int ==> ' num2str(intt)])
disp(' ')
T2=table(char(LED),[Rn Ro],[yeni eski],[ery ere],'VariableNames',{'LED color','Applied R(N&O)','Calculated Int(N&O)','Relative %err(N&O)'});
disp(T2)
%% Kırmızı veyvlengt şiftin
close
for i=1:nLED
    plot3(X,ones(length(X))*i,yavg_new(:,10+i),'LineWidth',1.5);
    ztxt=max(yavg_new(:,10+i));
    ytxt=i;
    xtxt=640;
    hold on
%     set(gca,'yscale','log')
    xlim([500 700])
    ttxt{i,1}=num2str(round(Res_new(3,i),1));
    ltxt{i,1}=['<==  ' num2str(round(int_new(3,i))) ' \gamma /ms'];
    text(xtxt,ytxt,ztxt,ltxt{i})
end,grid on
set(gca,'XMinorTick','on', 'ZMinorTick','on')
set(gca,'XMinorGrid','on', 'ZMinorGrid','on')
ylim([0.5 5.5])
yticks([1:5])
yticklabels(ttxt)
view([-35,35])
xlabel('Dalgaboyu (nm)')
ylabel('Uygulanan Direnç (\Omega)','Position',[451.9295,2.1335,-55.7400])
zlabel('Parlaklık (\gamma /ms.nm)')
set(gcf,'position',[175,345,870,420])
exportgraphics(gcf,'Istasyon1_kayma.png','Resolution',600),close
run intensity_control.m