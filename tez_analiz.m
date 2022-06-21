%% Genel veri okuma ve analizi
clear
Res=600;
z1=now;
source='irrad\';
xx=dir('irrad');
dirrad='04-May-2021 14:00:00';
dLED='04-May-2021 17:00:00';
[~,index] = sortrows([xx.datenum].'); xx = xx(index); clear index
nfiles=length(xx); asd=0;
for i=1:nfiles
    if xx(i).isdir==1
        xtxt=xx(i).name;
        xtxt=xtxt(1);
        if xtxt~='.'
            asd=asd+1;
            num(asd)=i;
            tarih{asd,1}=xx(i).name;
        end
    end
end
%% Zamanın yazımı ve sıralanması
nt=length(tarih);
for i=1:nt
    xx2=dir([source tarih{i,1}]);
    [~,index]=sortrows({xx2.datenum}.'); xx2=xx2(index); clear index
    asd=[xx2(1).datenum];
    dat(i,1)=asd;
end
t = datetime(dat,'ConvertFrom','datenum');
t2=datenum(t);
t2(:,2)=linspace(1,nt,nt);
t2=sortrows(t2); asd=1;
for i=t2(:,2)
    tarih2{asd,1}=tarih(i); asd=asd+1;
end
for i=1:nt
    tarih(i,1)=tarih2{1,1}(i);
end
for i=1:nt
    zaman(i,1)=round((t2(i,1)-t2(1,1)));
    zaman_2(i,1)=(t2(i,1)-t2(1,1));
end
%% İsimlerin yazımı
asd=(join([source tarih{14} '\clean.txt'],''));
% asd=(join([source tarih{1} '\clean.txt'],''));
Dcl=importdata(asd);
Dcl_3=importdata('C:\Users\gantu\OneDrive\Belgeler\MATLAB\Master\proje\LED\08102021\clean3.txt');
Dcl_4=importdata('C:\Users\gantu\OneDrive\Belgeler\MATLAB\Master\proje\LED\08102021\clean4.txt');
abc={'a','b','c'}; doz={'1','2'}; pen=[1:5]; asd=0;
for i=1:2
    for j=1:3
        if j==3
            for k=1:5
                asd=asd+1;
                names{asd,1}=[doz{i} abc{j} 'P' num2str(pen(k))];
            end
        else
            for k=2
                asd=asd+1;
                names{asd,1}=[doz{i} abc{j} 'P' num2str(pen(k))];
            end
        end
    end
end
%% D importu
n=length(names);
day=length(tarih);
daytime=1:day;
for i=1:day
    xx2=dir(join([source tarih{i}],''));
    Dates{i,1}=xx2(3,:).date;
    for j=1:n
       D=importdata(join([source tarih{i} '\' names{j} '.txt'],''));
       listo(j+n*(i-1),:)=join([tarih{i} '-' names{j}],'');
       Dx(:,1)=D.data(:,1); Dy(:,j+n*(i-1))=D.data(:,2);
    end
end
datetxt=char(Dates);
nday=4;
hour=round(linspace(0,zaman(end,end),nday));
hourtick=(linspace(1,day,nday));
names2={'Karanlık' 'Çevresel' 'Beyaz' 'Morötesi' 'Mavi' 'Yeşil' 'Kırmızı'};
for i=1:n/2
    names3{i}=names2{i};
    names3{i+n/2}=names2{i};
end
renkmat=[0 0 0; .5 .5 .5; 0 1 1; 1 0 1; 0 0 1; 0 1 0; 1 0 0;];
mark=zeros(n,3);
mark(1:n/2,:)=renkmat(1:n/2,:);
mark(n/2+1:n,:)=renkmat(1:n/2,:);
for i=1:n
    if i<=7
        rad(i)={'3.5 kGy Doz ile Işınlanmış Örnek'};
    else
        rad(i)={'7.0 kGy Doz ile Işınlanmış Örnek'};
    end
end
%% Y, X, Yc lerin yazılması
NLED=length(Dx);
for i=1:NLED
    if Dx(i)<=1500
        e1500=i;
    end
end
asd=0;
for i=1:e1500
    asd=asd+1;
    X(asd,1)=Dx(i);
    Y(asd,:)=Dy(i,:);
    Yc(asd,1)=Dcl.data(i,2);
    Yc_3mm(asd,1)=Dcl_3.data(i,2);
    Yc_4mm(asd,1)=Dcl_4.data(i,2);
end
n2=length(X);
%% Y nin düzenlenmesi
for i=2:e1500-1
    for j=1:day*n
        if Y(i,j)>=-1,else
            Y(i,j)=mean([Dy(i+1,j) Dy(i-1,j)]);
%             Y(i,j)=Dy(i-1,j);
        end
    end
end
for i=10:14
    Y(522,i)=Y(523,i);
    Y(632,i)=Y(633,i);
end
for i=1:n2
    if X(i)<=1450
        e1450=i;
    end
    if X(i)<=1000
        e1000=i;
    end
    if X(i)<=450
        e450=i;
    end
    if X(i)<=420
        e420=i;
    end
    if X(i)<=350
        e350=i;
    end
    if X(i)<=340
        e340=i;
    end
    if X(i)<=315
        e315=i;
    end
end
jet1=jet;
i0=round(linspace(1,256,day));
for i=1:day
    marker60(i,:)=jet1(i0(i),:);
    lgd60d{i,1}=[num2str(zaman(i)) '. Gün'];
    lgd60h{i,1}=[num2str(round(zaman_2(i)*24)) '. Saat'];
end
%% Y Spectra
for i=1:n
    for j=1:day
        renk=i+n*(j-1); pencere=subplot(2,7,i); 
        if j==day
            pencere.Position=pencere.Position+[0 0 0.015 -0.05];
        end
        plot(X,Y(:,renk),'Color',marker60(j,:));hold on
        set(gca,'fontsize',10)
        set(gca, 'YMinorGrid','on', 'XMinorGrid','on')
        set(gca, 'YMinorTick','on', 'XMinorTick','on')
        xi=400; xf=1300; xm=round(mean([xi xf])); yticks([0:10:60])
        if i>=3 && i<=7
            subtitle(names3{i},'Color','r')
        else
            subtitle(names3{i})
        end
        ylim([0 60]),xlim([200 1500]),xticks([xi xm xf])%,xtickangle(30)
        grid on
        if i==1 || i==1+n/2
            ylabel('Geçirgenlik (%)')
        else
            yticklabels([])
        end
        if i==4 || i==4+n/2
            xlabel('Dalgaboyu (nm)')
            title([rad{i} 'ler'])
        end
    end
end
set(gcf,'position',[25.8,277.8,1116,485.2])
exportgraphics(gcf,'Spectra_of_Y.png','Resolution',Res)
%% saat cinsinden legendler (Hours)
close
renk=11:n:n*day;
for i=1:day
    subplot(5,3,2)
    asd2(i,1)=plot(X,Y(:,renk(i)),'Color',marker60(i,:));hold on
end
set(gca,'fontsize',12)
set(gcf,'position',[25.8,343,1468.8,420])
legendo1=legend([asd2(:)],lgd60h);
legendo1.NumColumns=10;
% konum=legendo1.Position;
legendo1.Position=[0.1269,0.7262,0.7713,0.2555];
legendo1.FontSize=10;
exportgraphics(gcf,'Hours.png','Resolution',Res)
%% gün cinsinden legendler (Days)
close
renk=11:n:n*day;
for i=1:day
    subplot(5,3,2)
    asd2(i,1)=plot(X,Y(:,renk(i)),'Color',marker60(i,:));hold on
end
set(gca,'fontsize',12)
set(gcf,'position',[25.8,343,1468.8,420])
yticklabels([])
xticklabels([])
legendo1=legend([asd2(:)],lgd60d);
legendo1.NumColumns=10;
legendo1.Position=[0.1529,0.6881,0.6740,0.2555];
legendo1.FontSize=10;
exportgraphics(gcf,'Days.png','Resolution',Res)
%% Kırıkların giderilmesi (Y2)
close
Y2=Y; k1=521; k2=631;
for i=1:day*n
    Y2(1:k1,i)=Y(1:k1,i)*Y(k1+1,i)/Y(k1,i);
    Yc(1:k1)=Yc(1:k1)*Yc(k1+1)/Yc(k1);
    Yc_3mm(1:k1)=Yc_3mm(1:k1)*Yc_3mm(k1+1)/Yc_3mm(k1);
    Yc_4mm(1:k1)=Yc_4mm(1:k1)*Yc_4mm(k1+1)/Yc_4mm(k1);
end
Y3=Y2;
for i=1:day*n
    Y3(1:k2,i)=Y2(1:k2,i)*Y2(k2+1,i)/Y2(k2,i);
    Yc(1:k2)=Yc(1:k2)*Yc(k2+1)/Yc(k2);
    Yc_3mm(1:k2)=Yc_3mm(1:k2)*Yc_3mm(k2+1)/Yc_3mm(k2);
    Yc_4mm(1:k2)=Yc_4mm(1:k2)*Yc_4mm(k2+1)/Yc_4mm(k2);
end
Y2=Y3; clear Y3
%%
maksu=95;
ss1=subplot(1,2,1);
ps1=plot(X,Yc_3mm,X,Y2(:,1),'LineWidth',2);
xlim([310 1500]),ylim([0 maksu]),grid on
set(gca, 'YMinorGrid','on', 'XMinorGrid','on')
set(gca, 'YMinorTick','on', 'XMinorTick','on')
% ss1.Position=ss1.Position+[-0.05 0 0.1 0];
ylabel('Geçirgenlik (%)'),x1=xlabel('Dalgaboyu (nm)');
x1.Position=x1.Position.*[1.86 1 1];
ss2=subplot(1,2,2);
ps2=plot(X,Y2(:,1)./Yc_3mm*100,'k','LineWidth',2);
xlim([310 1500]),ylim([0 maksu]),grid on
set(gca, 'YMinorGrid','on', 'XMinorGrid','on')
set(gca, 'YMinorTick','on', 'XMinorTick','on')
yticks([0:10:maksu]),yticklabels('')
% ss2.Position=ss2.Position+[-0.025 0 0.1 0];
l1=legend([ps1],'Temiz Örnek','Işınlanan Örnek','Location','ne');
l2=legend([ps2],'Işınlanmış ve Temiz Örneğin Oranı','Location','nw');
set(gcf,'position',[25.8,343,1116,420])
exportgraphics(gcf,'Samples.png','Resolution',Res)
%% Y2 Spectra
for i=1:n
    for j=1:day
        renk=i+n*(j-1); pencere=subplot(2,7,i); 
        if j==day
            pencere.Position=pencere.Position+[0 0 0.015 -0.05];
        end
        plot(X,Y2(:,renk),'Color',marker60(j,:));hold on
        set(gca,'fontsize',10)
        set(gca, 'YMinorGrid','on', 'XMinorGrid','on')
        set(gca, 'YMinorTick','on', 'XMinorTick','on')
        xi=400; xf=1300; xm=round(mean([xi xf])); yticks([0 10 20 30 40 50 60 70])
        if i>=3 && i<=7
            subtitle(names3{i},'Color','r')
        else
            subtitle(names3{i})
        end
        ylim([0 70]),xlim([200 1500]),xticks([xi xm xf])%,xtickangle(30)
        grid on
        if i==1 || i==1+n/2
            ylabel('Geçirgenlik (%)')
        else
            yticklabels([])
        end
        if i==4 || i==4+n/2
            xlabel('Dalgaboyu (nm)')
            title([rad{i} 'ler'])
        end
    end
end
set(gcf,'position',[25.8,277.8,1116,485.2])
exportgraphics(gcf,'Spectra_of_Y2.png','Resolution',Res)
%% Minlerin bulunması
close
for i=1:n2
    if X(i)==1100
        s1=i;
    elseif X(i)==1300
        s2=i;
    end
end
Xlocal=(1100:1300)';
minler=zeros(n,day);
for i=1:n
    for j=1:day
        k=i+n*(j-1);
%         asdf(i,j)=a;
        asd=Y2(s1:s2,k);
        fit1=fit(Xlocal,asd,'poly3');
        [asd1, asd2]=fminbnd(fit1,1100,1300);
        minler(i,j)=asd2;
    end
end
%% LocalMinler
clear asd asd2
for i=[[3:7],[3:7]+7]
    if i<=n/2
        asd=1;
    else
        asd=2;
    end
    subplot(1,2,asd),p(i)=plot(zaman,minler(i,:),'--o','color',mark(i,:));
    ylim([46 52]),ylabel('Geçirgenlik (%)')
    if asd==2
        title('7.0 kGy Işınlanmış Örnekler')
    else
        title('3.5 kGy Işınlanmış Örnekler')
    end
    xlabel('Işınlamadan Sonra Geçen Süre (saat)')
    set(gca, 'YMinorTick','on', 'XMinorTick','on')
    set(gca, 'YMinorGrid','on', 'XMinorGrid','on')
    grid on,hold on
end
legend(p(10:14),names2(3:7),'Location','best')
set(gcf,'position',[25.8,343,1116,420])
exportgraphics(gcf,'LocalMinPoints.png','Resolution',Res)
%% AvgHistler
close
meanmin1=mean(minler([3,4,6,7],:));
meanmin2=mean(minler([3,4,6,7]+7,:));
meanrate1=minler(5,:)./meanmin1;
meanrate2=minler(5+7,:)./meanmin2;
nbar=40;
s=subplot(1,2,1); %s.Position=s.Position+[0 0 0.08 0];
H1=histogram(meanrate1,nbar);
He1=H1.BinEdges; Hy1=H1.Values';
title('3.5 kGy Işınlanmış Örnekler')
ylabel('Sayım Sayısı'),xlabel('Değer')
set(gca, 'YMinorTick','on', 'XMinorTick','on')
set(gca, 'YMinorGrid','on', 'XMinorGrid','on')
ylim([0 12]),grid on
for i=1:length(Hy1)
    Hx1(i,1)=mean([He1(i),He1(i+1)]);
end
s=subplot(1,2,2); %s.Position=s.Position+[0 0 0.08 0];
H2=histogram(meanrate2,nbar);
He2=H2.BinEdges; Hy2=H2.Values';
title('7.0 kGy Işınlanmış Örnekler')
ylabel(''),xlabel('Değer')
set(gca, 'YMinorTick','on', 'XMinorTick','on')
set(gca, 'YMinorGrid','on', 'XMinorGrid','on')
ylim([0 12]),grid on
yticklabels('')
set(gcf,'position',[25.8,343,1116,420])
exportgraphics(gcf,'AvgHist.png','Resolution',Res)
for i=1:length(Hy2)
    Hx2(i,1)=mean([He2(i),He2(i+1)]);
end, clear He1 He2
%% Histfitler
close
for i=1:nbar-1
    if Hx1(i)<=.98 && Hx1(i+1)>=.98
        i1=i;
    elseif Hx1(i)<=1.01 && Hx1(i+1)>=1.01
        i2=i;
    end
end
f1=fit(Hx1(i1:i2),Hy1(i1:i2),'gauss1');
abc1=coeffvalues(f1);
s=subplot(1,2,1);hold on,histogram(meanrate1,nbar),p1=plot(f1);
% s.Position=s.Position+[0 0 0.08 0];
text(0.96, 11, ['\mu_1 = ' num2str(round(abc1(2),4))])
ylim([0 12]),xlim([0.95 1.05]), set(p1,'linewidth',1.5),grid on
xticks([0.96:0.02:1.04]),set(gca,'box','on')
set(gca, 'YMinorTick','on', 'XMinorTick','on')
set(gca, 'YMinorGrid','on', 'XMinorGrid','on')
title('3.5 kGy Işınlanmış Örnekler')
ylabel('Sayım Sayısı'),xlabel('Değer')
legend('Off')
for i=1:nbar-1
    if Hx2(i)<=.985 && Hx2(i+1)>=.985
        i1=i;
    elseif Hx2(i)<=1.015 && Hx2(i+1)>=1.015
        i2=i;
    end
end
f2=fit(Hx2(i1:i2),Hy2(i1:i2),'gauss1');
abc2=coeffvalues(f2);
s=subplot(1,2,2);hold on,histogram(meanrate2,nbar),p2=plot(f2);
% s.Position=s.Position+[0 0 0.08 0];
text(0.96, 11, ['\mu_2 = ' num2str(round(abc2(2),4))])
ylim([0 12]),yticks([0:2:12]),ylabel('Sayım Sayısı')
xlim([0.95 1.05]), set(p2,'linewidth',1.5),grid on
xticks([0.96:0.02:1.04]),set(gca,'box','on')
set(gca, 'YMinorTick','on', 'XMinorTick','on')
set(gca, 'YMinorGrid','on', 'XMinorGrid','on')
title('7.0 kGy Işınlanmış Örnekler')
xlabel('Değer')
legend('Histogram','Gauss Fit')
set(gcf,'position',[25.8,343,1116,420])
exportgraphics(gcf,'FitHist.png','Resolution',Res)
%% Eğrilerin düzenlenmesi (Y3)
for i=1:n
    for j=1:day
        k=i+n*(j-1);
        if i==5
%             asd=abc1(2)*mean(exphour(:,j))/(minler(i,j));
%             asd=abc1(2)*meanmin1(j)/(minler(i,j));
            asdx=abc1(2)/meanrate1(j);
        elseif i==5+n/2
%             asd=abc2(2)*mean(exphour(:,j))/(minler(i,j));
%             asd=abc2(2)*meanmin2(j)/(minler(i,j));
            asdx=abc2(2)/meanrate2(j);
        else
            asdx=1;
%         elseif i<=7 && i~=5
%             asdx=meanmin1(j)/(minler(i,j));
%         elseif i>7 && i~=5+7
%             asdx=meanmin2(j)/(minler(i,j));
        end
        Y3(:,k)=Y2(:,k)*asdx;
    end
end
%% Localminler 2
close
for i=1:n2
    if X(i)==1100
        s1=i;
    elseif X(i)==1300
        s2=i;
    end
end
Xlocal=(1100:1300)';
minler=zeros(n,day);
for i=1:n
    for j=1:day
        k=i+n*(j-1);
%         asdf(i,j)=a;
        asd=Y3(s1:s2,k);
        fit1=fit(Xlocal,asd,'poly3');
        [asd1, asd2]=fminbnd(fit1,1100,1300);
        minler(i,j)=asd2;
    end
end
clear asd asd2
for i=[[3:7],[3:7]+7]
    if i<=n/2
        asd=1;
    else
        asd=2;
    end
    subplot(1,2,asd),p(i)=plot(zaman,minler(i,:),'--o','color',mark(i,:));
    ylim([46 52]),ylabel('Geçirgenlik (%)')
    if asd==2
        title('7.0 kGy Işınlanmış Örnekler')
    else
        title('3.5 kGy Işınlanmış Örnekler')
    end
    xlabel('Işınlamadan Sonra Geçen Süre (saat)')
    set(gca, 'YMinorTick','on', 'XMinorTick','on')
    set(gca, 'YMinorGrid','on', 'XMinorGrid','on')
    grid on,hold on
end
legend(p(10:14),names2(3:7),'Location','best')
set(gcf,'position',[25.8,343,1116,420])
exportgraphics(gcf,'LocalMinPoints2.png','Resolution',Res)
%% Y3 Spectra
for i=1:n
    for j=1:day
        renk=i+n*(j-1); pencere=subplot(2,7,i); 
        if j==day
            pencere.Position=pencere.Position+[0 0 0.015 -0.05];
        end
        plot(X,Y3(:,renk),'Color',marker60(j,:));hold on
        set(gca,'fontsize',10)
        set(gca, 'YMinorGrid','on', 'XMinorGrid','on')
        set(gca, 'YMinorTick','on', 'XMinorTick','on')
        xi=400; xf=1300; xm=round(mean([xi xf])); yticks([0:10:70])
        if i>=3 && i<=7
            subtitle(names3{i},'color','r')
        else
            subtitle(names3{i})
        end
        ylim([0 70]),xlim([200 1500]),xticks([xi xm xf])%,xtickangle(30)
             set(gca, 'YMinorTick','on', 'XMinorTick','on')
        grid on
        if i==1 || i==1+n/2
            ylabel('Geçirgenlik (%)')
        else
            yticklabels([])
        end
        if i==4 || i==4+n/2
            xlabel('Dalgaboyu (nm)')
            title([rad{i} 'ler'])
        end
    end
end
set(gcf,'position',[25.8,277.8,1116,485.2])
exportgraphics(gcf,'Spectra_of_Y3.png','Resolution',Res)
%% Oranların yazılması (YRate)
close
for i=1:n
    for j=1:day
        k=i+n*(j-1);
        if i>=3 && i<=7
            YRate(:,k)=Y3(:,k)./Yc_4mm;
        else
            YRate(:,k)=Y3(:,k)./Yc_3mm;
        end
    end
end
YRate2=YRate;
YRate2(:,8)=YRate(:,8)*0.98;
%% YRate Spectra
for i=1:n
    for j=1:day
        renk=i+n*(j-1); pencere=subplot(2,7,i); 
        if j==day
            pencere.Position=pencere.Position+[0 0 0.015 -0.05];
        end
        plot(X,YRate(:,renk)*100,'Color',marker60(j,:))
        set(gca,'fontsize',10);hold on
        set(gca, 'YMinorGrid','on', 'XMinorGrid','on')
        set(gca, 'YMinorTick','on', 'XMinorTick','on')
        xi=440; xf=900; xm=round(mean([xi xf])); yticks([0:15:90])
        if i>=3 && i<=7
            subtitle(names3{i},'Color','r')
        else
            subtitle(names3{i})
        end
        ylim([0 90]),xlim([340 1000]),xticks([xi xm xf])%,xtickangle(30)
        grid on
        if i==1 || i==1+n/2
            ylabel('Geçirgenlik Oranı (%)')
        else
            yticklabels([])
        end
        if i==4 || i==4+n/2
            xlabel('Dalgaboyu (nm)')
            title([rad{i} 'ler'])
        end
    end
end
% legend('boxoff')
% lgnd=legend([pp(14,ord6(1:6)-1)],lgd,...
%     'position',[0.15,0.91,1.0207,0.0717]);
% lgnd.NumColumns = 3;
set(gcf,'position',[25.8,277.8,1116,485.2])
exportgraphics(gcf,'Spectra_of_YRate.png','Resolution',Res)
%% Harita yazma
close
spex=(e340):(e1500);
zaman_rev=zaman(day:-1:1);
hsv2=hsv;
hsv2(end,:)=[0 0 0];
bitto=64;
faz=linspace(0.01,10,bitto);
for i=1:bitto
    ang(i,1)=0.1364*(i-1)*pi/(bitto-1);
    if ang(i)<.116*pi
        RenkPaleti(i,3)=sin(ang(i)*faz(i)+pi/2);
    else
        RenkPaleti(i,3)=-1;
    end
        RenkPaleti(i,2)=sin(ang(i)*faz(i));
    if ang(i)<0.3663
        RenkPaleti(i,1)=sin(ang(i)*faz(i)-pi/2);
    else
        RenkPaleti(i,1)=1;
    end
end
RenkPaleti=RenkPaleti/2+0.5;
for i=1:n
    renk=linspace(i,i+n*(day-1),day); asdasd=subplot(2,7,i);
    asdasd.Position=asdasd.Position+[0 0 0.012 -0.038];
    imagesc(X(spex),zaman_rev,YRate(spex,renk)'*100)
    set(gca,'fontsize',10)
    caxis([0 90]),ylim([zaman(1) zaman(day)])
    colormap(RenkPaleti)

    xlim([340 1000])
    grid on
    xlim([340 1000]),xi=440; xf=900; xm=round(mean([xi xf]));
    if i>2 && i<=7
        subtitle(names3{i},'color','r')
    else
        subtitle(names3{i})
    end
    xlim([340 1000]),xticks([xi xm xf])%,xtickangle(30)
    if i==1 || i==1+n/2
        ylabel(['Işınlama Sonrası'; 'Geçen Süre (gün)'])
        yticks([0:25:150]+3)
        yticklabels({'150' '' '100' '' '50' '' '0'})
    else
        yticks([0:25:150]+3)
        yticklabels([])
    end
    if i==4 || i==4+n/2
        xlabel('Dalgaboyu (nm)')
        title([rad{i} 'ler'])
    end
    if i==n
        c = colorbar;
        c.Label.String = 'Göreli Geçirgenlik (%)';
        c.Position=[0.926,0.1105,0.01286,0.7762];
    end
end
set(gcf,'position',[25.8,277.8,1116,485.2])
exportgraphics(gcf,'Map_of_YRate.png','Resolution',Res),
%% İntegraller (int), hasar (Dmg) ve hata (Yerr)
int=zeros(n,day);
close
ara=e340:e1000;
for i=1:n
    for j=1:day
        a=i+n*(j-1);
        for k=ara
            int(i,j)=int(i,j)+mean([YRate(k-1,a) YRate(k,a) YRate(k,a) YRate(k+1,a)])*(X(k+1)-X(k))*100;
        end
    end
end
Dmg=100-int/((ara(end)-ara(1))*(X(2)-X(1)));
for i=ara
    Yp(i,:)=(YRate(i+1,:)-YRate(i,:))/(X(i+1)-X(i));
end
n2=length(ara);
Yerr=zeros(n,day);
for i=1:n
    for j=1:day
        a=i+n*(j-1);
        Yerr(i,j)=abs(((ara(end)-ara(1))^3).*(Yp(ara(end),a)-Yp(ara(1),a)))/(12*n2^2);
    end
end
%% Dmg spectra
for i=1:n
    if i<=n/2
        asd=1;
    else
        asd=2;
    end
    subplot(1,2,asd),plot(zaman,Dmg(i,:),'--o','color',mark(i,:))
%     set(gca,'fontsize',12)
    ylim([15 60])
    if asd==2
        legend(names2(:))
        yticklabels(''),yticks([10:10:70]),title('7.0 kGy Işınlanmış Örnekler')
    else
        legend(names2(:))
        ylabel('Toplam Geçirgenlik Kaybı (%)'),title('3.5 kGy Işınlanmış Örnekler')
    end
    xlabel('Işınlama Sonrası Geçen Süre (gün)')
    xlim([0 170]),ylim([10 70])
    set(gca, 'YMinorGrid','on', 'XMinorGrid','on')
    set(gca, 'YMinorTick','on', 'XMinorTick','on')
    grid on,hold on
end
set(gcf,'position',[25.8,343,1116,420])
exportgraphics(gcf,'TotalDamage.png','Resolution',Res)
%% Hasar fiti
close
n4=1e3; x=linspace(0,160,n4);
con=0; unc=0;
clear con unc
cn1=0; cn2=0; D=0;
aa=zeros(n,5);
aa(:,1)=10*ones(n,1);
aa(:,2)=5*ones(n,1);
aa(:,3)=30*ones(n,1);
aa(:,4)=15*ones(n,1);
aa(:,5)=10*ones(n,1);
aa(4,2)=3.5;
aa(11,5)=17;
aa(12,2)=2.25;
aa(12,4)=40;
myfit=fittype(@(a,alfa,b,beta,c,x)(a*exp(-x/alfa)+b*exp(-x/beta)+c),'independent',{'x'},'coefficients',{'a','alfa','b','beta','c'});
fday=1;
for i=1:n
    fitto1=fit(zaman_2(fday:end),Dmg(i,fday:end)',myfit,'StartPoint', aa(i,:));
    con(i,:)=coeffvalues(fitto1);
    asdf=confint(fitto1); unc(i,:)=asdf(2,:)-coeffvalues(fitto1);
    subplot(2,7,i),plot(fitto1,zaman_2,Dmg(i,:)),legend off
    xticks([0 75 150]),set(gca, 'YMinorgrid','on', 'XMinorgrid','on')
    set(gca, 'YMinorTick','on', 'XMinorTick','on')
    ylabel(''),xlabel(''),ylim([10 60]),subtitle(names3{i}),grid on
    if i==4 || i==11
        xlabel('Işınlama Sonrası Geçen Süre (gün)'),title([rad{i} 'ler'])
    elseif i==7
        legend('Data','Fit','position',[0.8391,0.4639,0.0643,0.0717])
    elseif i==1 || i==8
        ylabel('Hasar Oranı (%)')
    end
    for j=1:n4
        yexp(i,j)=fitto1(x(j));
    end
    for j=1:day
        Ex(i,j)=fitto1(zaman_2(j));
        Chi2Dmg(i,j)=(Dmg(i,j)-Ex(i,j))^2/Ex(i,j);
    end
end
kikare(:,1)=sum(Chi2Dmg')';
kikarendf=kikare/5;
set(gcf,'position',[25.8,343,1116,420])
exportgraphics(gcf,'DamageFit.png','Resolution',Res)
%% DamageFit 3.5 ve 7.0 kGy
close
ekler={' Kutu' ' Işık' ' LED' ' LED' ' LED' ' LED' ' LED' ' Kutu' ' Işığı' ' LED' ' LED' ' LED' ' LED' ' LED'};
for i=1:n
    name_legend{i}=[rad{i}(1:7) ' ' names3{i} ekler{i}];
end
Yerr2=sqrt(Yerr.^2+1^2);
for i=1:n
    if i<=n/2
        subplot(1,2,1)
    else
        subplot(1,2,2)
    end
    er1(i)=errorbar(zaman,Dmg(i,:),Yerr2(i,:),'o',"color",mark(i,:));hold on
    p31=plot(x,yexp(i,:),"color",mark(i,:),'LineWidth',1.5);
    if i==5 || i==14
        p31.LineStyle='--';
    end
    ylabel('Toplam Geçirgenlik Kaybı (%)'),%title('3.5 kGy Işınlanmış Örnekler')
    xlabel('Işınlama Sonrası Geçen Süre (gün)')
    xlim([0 170]),ylim([10 70])%,yticks([10:5:70])
    set(gca, 'YMinorTick','on', 'XMinorTick','on')
    set(gca, 'YMinorGrid','on', 'XMinorGrid','on')
    grid on
    if i==7
        legend([er1(1:7)],name_legend{1:7})
%         exportgraphics(gcf,'DamageFit3p5kGy.png','Resolution',Res),close
    elseif i==14
        legend([er1(8:14)],name_legend{8:14})
    end
end
set(gcf,'position',[25.8,343,1116,420])
exportgraphics(gcf,'DamageFit.png','Resolution',Res)
%% Fit farkı
% for i=1:day
%     for ii=1:length(x)-1
%         if x(ii)<=zaman_2(i) && x(ii+1)>zaman_2(i);
%             xday(i)=ii;
%         end
%     end
% end
% for i=1:n
%     names3rev{i}=names3{15-i};
% end
% for i=1:n
%     if i<=7
%         subplot(1,2,1)
%         ii=7-i;
%     else
%         subplot(1,2,2)
%         ii=14-i;
%     end
%     kat=3;
%     erro=errorbar([1:60],ones(day,1)*kat*ii,zeros(day,1),abs(yexp(i,xday)-Dmg(i,:)),'.');hold on
%     erro.Color=mark(i,:); erro.CapSize=0; erro.LineWidth=2;
%     erro2=errorbar([1:60],ones(day,1)*kat*ii,zeros(day,1),Yerr2(i,:),'.');
%     erro2.Color=mark(i,:); erro2.CapSize=5;
%     xticks([0:10:60]),xticklabels(zaman([1:10:60,60])),xlim([0 61])
%     yticks([0:kat:kat*6]),yticklabels(names3rev),ytickangle(45),ylim([-1 kat*7])
% end
% set(gcf,'position',[25.8,343,1116,420])
%% Tablo yazma (expfit)
close
for i=1:n
    for j=1:5
        ABCstr{i,j}=[num2str(round(con(i,j),2)) ' ' char(177) ' ' num2str(round(unc(i,j),2))];
    end
end
T_expfit=table(rad',names3');
titles=[{'Dose'},'Sample','Constant A (%)','Fast component (days)','Constant B (%)','Slow component (days)','Constant C (%)'];
for i=1:length(con(1,:))
    T_expfit{:,i+2}=ABCstr(:,i);
end
for i=1:length(titles)
    T_expfit.Properties.VariableNames{i} = titles{i};
end
%% Tablo yazma (Dmgler)
DMGstr{1,1}='Örnekler';
DMGstr{2,1}='Zaman (saat)';
for i=1:n
    if i==4 || i==4+7
        DMGstr{1,i+1}=rad{i};
    end
    DMGstr{2,i+1}=[names3{i}];
    for j=1:day
        if length(num2str(round(Dmg(i,j),2)))==4
            ilki=[num2str(round(Dmg(i,j),2)) '0'];
        elseif length(num2str(round(Dmg(i,j),2)))==2
            ilki=[num2str(round(Dmg(i,j),2)) '.00'];
        else
            ilki=num2str(round(Dmg(i,j),2));
        end
        if length(num2str(round(Yerr(i,j),2)))==3
            ikinci=[num2str(round(Yerr(i,j),2)) '0'];
        elseif length(num2str(round(Yerr(i,j),2)))==1
            ikinci=[num2str(round(Yerr(i,j),2)) '.00'];
        else
            ikinci=num2str(round(Yerr(i,j),2));
        end
        DMGstr{j+2,i+1}=[ilki ' ' char(177) ' ' ikinci];
        DMGstr{j+2,1}=[num2str(zaman(j))];
        DMGstr{day+3,1}='Chi2';
        DMGstr{day+4,1}='Chi2/ndf';
    end
    DMGstr{day+3,i+1}=num2str(round(kikare(i),4));
    DMGstr{day+4,i+1}=num2str(round(kikarendf(i),4));
end
T_dmg=table(num2str(zaman));
for i=2:15
    titles2{i-1,1}=[rad{i-1}(1:7) ' ' names3{i-1}];
    T_dmg(:,i)=table(DMGstr(3:62,i));
end
T_dmg.Properties.VariableNames{1}='Zaman (saat)';
for i=1:n
    T_dmg.Properties.VariableNames{i+1} = titles2{i};
end
% T_dmg=table(DMGstr(3:62,2:15));
%% SpectraLED
link_LED='C:\Users\gantu\OneDrive\Belgeler\MATLAB\Master\proje\LED\m19\dot';
files=dir([link_LED '\*.txt']);
clr={'Mavi','Yeşil','Kırmızı','Morötesi','Beyaz'}; var={'_new','_old'};
mark_LED={'b','g','r','m','c'};
Nfile=length(files);
LED=load([link_LED '\' files(1).name]); LEDx=LED(:,1);
for i=1:Nfile
    LED=load([link_LED '\' files(i).name]);
    LEDy(:,i)=LED(:,2);
end
a=0;
for i=1:length(clr)
    for j=1:length(var)
        a=a+1;
        names_LED{a,1}=[clr{i} var{j}];
    end
end
NLED=length(LEDx);
lvl=560;
for j=1:Nfile
    for i=1:NLED
        if LEDy(i,j)>=lvl
            LEDy2(i,j)=LEDy(i,j)-lvl;
        else
            LEDy2(i,j)=0;
        end
    end
end
for i=2:NLED-1
   for j=1:Nfile
       if LEDy2(i-1,j)==0 && LEDy2(i+1,j)==0 && LEDy2(i,j)>0
          LEDy2(i,j)=0;
          LEDy(i,j)=400+20*(2*rand-1);
       end
   end
end
Ngrup=50;
for i=1:Nfile/Ngrup
    LEDmean2(:,i)=mean(LEDy2(:,1+Ngrup*(i-1):Ngrup*(i))');
    LEDmean(:,i)=mean(LEDy(:,1+Ngrup*(i-1):Ngrup*(i))');
end
LEDint=0;
for i=1:NLED-1
    LEDint=LEDint+sum(LEDmean2/100)*(LEDx(i+1)-LEDx(i));
end
for i=1:Nfile/Ngrup
    txt{i,1}=[ '==> ' num2str(round(LEDint(i)/1000)) ' \gamma /s'];
end
for i=1:Nfile/Ngrup
    plot(LEDx,LEDmean(:,i)/100,mark_LED{i},"LineWidth",2),hold on
    set(gca, 'YMinorTick','on', 'XMinorTick','on')
    set(gca, 'YMinorGrid','on', 'XMinorGrid','on')
end,set(gca,'yscale','log'),ylim([4 2e3]),xlim([340 760]),grid on
set(gca, 'XMinorTick','on'), text(578.6655,758.9004,txt)
ylabel('Parlaklık (\gamma /ms.nm)'),xlabel('Dalgaboyu (nm)'),legend boxoff
ll2=legend(clr,'location','n'); %ll2.Position=ll2.Position
set(gcf,'position',[175,345,870,420])
exportgraphics(gcf,'SpectraLED.png','Resolution',Res)
%% FitLED
close
xLED=linspace(200,1000,1000);
LED2=LEDmean/100;
for i=1:5
    if i<=3
        subplot(2,3,i)
    else
        subplot(2,2,i-1)
    end
    fLED=fit(LEDx,LED2(:,i),'gauss5');plot(fLED,LEDx,LED2(:,i)),legend off
    set(gca, 'XMinorTick','on'),set(gca,'fontsize',10)
    xlim([300 800]),ylim([4 1.2e3]),set(gca,'yscale','log'),grid on
    title(clr{i}), yticks(10.^[1,2,3])
    ylabel('Parlaklık (\gamma /ms.nm)')
    if i==2
        xlabel('Dalgaboyu (nm)')
    elseif i==5
        xlabel('Dalgaboyu (nm)','position',[220,1.8216,-1])
    else
        xlabel('')
    end
    asd=1;
    for ii=2:length(xLED)-1
        if fLED(xLED(ii-1))<=fLED(xLED(ii)) && fLED(xLED(ii+1))<=fLED(xLED(ii)) && fLED(xLED(ii))>=10
            Ly(i,asd)=fLED(xLED(ii)); Lx(i,asd)=xLED(ii); asd=asd+1;
        end
    end
    ytext=8e2; xtext=325;
    ttt=text(xtext,ytext,[num2str(round(Lx(i,1),2)),' nm']);
    ttt.FontSize=12;
    % beyazın ortalaması alınırsa
    if length(Lx(1,:))>1   
        if Lx(i,2)>0
            ttt=text(xtext,ytext/2,[num2str(round(Lx(i,2),2)),' nm']); ttt.FontSize=12;
        end
    end
    if i==3
        ll=legend('Data','Fit');ll.NumColumns=2;
        ll.Position=[0.7812,0.4916,0.1229,0.04190];
    end
end
for i=1:length(Lx)
    if Lx(i,2)==0
        Lm(i,1)=Lx(i,1);
    else
%         Lm(i,1)=(Lx(i,1)*Ly(i,1)+Lx(i,2)*Ly(i,2))/(Ly(i,1)+Ly(i,2)); % beyaz var ise
    end
end
% ttt=text(xtext,ytext/4,['Ağırlıklı ortalama: ' num2str(round(Lm(i,1),2)),' nm']); ttt.FontSize=12;
set(gcf,'position',[25.8,279.4,1116,483.6])
exportgraphics(gcf,'FitLED.png','Resolution',Res)
%% Fit parameters as a function of the LED wavelength
% uv b w g r
close
% ord=[4 5 3 6 7]; % beyaz varken
ord=[4 5 6 7]; % beyaz yokken
ord2=[2,4,5];
isimler={'FastComp.png','SlowComp.png','PermanentDamage.png'};
yler={'\tau_{hızlı} (gün)','\tau_{yavaş} (gün)','C sabiti (%)'};
ylimler=[0 5; 10 70; 10 40];
Lms=sort(Lm);
for i=1:3
    subplot(1,3,i)
    errorbar(Lms,con(ord,ord2(i)),unc(ord,ord2(i)),'--k','LineWidth',1.5),hold on,ylim(ylimler(i,:))
    errorbar(Lms,con(ord+7,ord2(i)),unc(ord+7,ord2(i)),'--r','LineWidth',1.5),xlim([350 700])
    set(gca,'fontsize',12);
    for ii=1:length(ord)
        ttt=text(Lms(ii),0.95*con(ord(ii)+7,ord2(i))-unc(ord(ii)+7,ord2(i)),names2{ord(ii)});
        if i==1
            ttt.Rotation=-45;
        end
    end
    set(gca, 'YMinorTick','on', 'XMinorTick','on')
    set(gca, 'YMinorGrid','on', 'XMinorGrid','on'),grid on
    xlabel('Dalgaboyu (nm)'),ylabel(yler{i}),legend('3.5 kGy','7.0 kGy')
%     close
end
set(gcf,'position',[46.6,343,1472.8,420])
exportgraphics(gcf,'sabitler.png','Resolution',Res)
%% Ratio of Integrated Transmittance Loss
close
for i=1:7
    name_legend2{i}=[names2{i} ekler{i}];
end
close
x2=linspace(5,160,length(x));
for i=1:n/2
    subplot(1,1,1)
    ITL(i,:)=(Dmg(i+n/2,:)./Dmg(i,:));
    ITLerr(i,:)=ITL(i,:).*sqrt((Yerr(i,:)./Dmg(i,:)).^2+(Yerr(i+n/2,:)./Dmg(i+n/2,:)).^2);
%     myfit2=fittype(@(a,b,x)(a+x*b),'independent',{'x'},'coefficients',{'a' 'b'});
    myfit2=fittype(@(a,x)(a+x*0),'independent',{'x'},'coefficients',{'a'});
    f3=fit(zaman_2(5:end),ITL(i,5:end)',myfit2,'StartPoint',[1]);
%     f3=fit(zaman_2(4:end),ITL(i,4:end)','poly1');
    con3(i,:)=coeffvalues(f3); asdf=confint(f3); unc3(i,:)=abs(con3(i,:)-asdf(1,:));
    ylabel('Hasar_{7.0 kGy} / Hasar_{3.5 kGy}'),hold on
    for j=5:day
        Lin(i,j)=f3(zaman_2(j));
        Chi2ITL(i,j)=(ITL(i,j)-Lin(i,j))^2/ITLerr(i,j);
    end
    er2(i)=errorbar(zaman,ITL(i,:),ITLerr(i,:),'o','MarkerSize',3,"color",renkmat(i,:),'CapSize',1);hold on
    plot(x2,ones(length(x2),1)*con3(i,1),"color",renkmat(i,:),'LineWidth',1.5)
    xlabel('Işınlama Sonrası Geçen Süre (gün)'),ylim([0.9 1.5]),xlim([0 170])
    yticks([0.9:0.1:1.5])%,xticks([1:10, 20:10:100])
    set(gca,'box','on'),xticks(2.^[1:10])
    set(gca, 'YMinorTick','on', 'XMinorTick','on')
    set(gca, 'YMinorGrid','on', 'XMinorGrid','on')
    set(gca,'xscale','log'),grid on
end
kikare2(:,1)=sum(Chi2ITL')';
kikarendf2=kikare2;
lgnd2=legend([er2(1:n/2)],name_legend2); lgnd2.NumColumns=2;
set(gcf,'position',[175,345,870,420])
exportgraphics(gcf,'IntegralTransmittanceLoss.png','Resolution',Res)
ITLstr{1,1}='Fit Değeri';
ITLstr{2,1}='Zaman (saat)';
T_itl=0; clear T_itl
for i=1:n/2
    ITLstr{1,i+1}=names2{i};
    ITLstr{2,i+1}=[num2str(round(con3(i,1),4)) ' ' char(177) ' ' num2str(round(unc3(i,1),4))];
    for j=1:day
        if length(num2str(round(ITL(i,j),4)))==5
            ilki=[num2str(round(ITL(i,j),4)) '0'];
        elseif length(num2str(round(ITL(i,j),4)))==4
            ilki=[num2str(round(ITL(i,j),4)) '00'];
        elseif length(num2str(round(ITL(i,j),4)))==3
            ilki=[num2str(round(ITL(i,j),4)) '000'];
        else
            ilki=num2str(round(ITL(i,j),4));
        end
        if length(num2str(round(ITLerr(i,j),4)))==5
            ikinci=[num2str(round(ITLerr(i,j),4)) '0'];
        elseif length(num2str(round(ITLerr(i,j),4)))==4
            ikinci=[num2str(round(ITLerr(i,j),4)) '00'];
        elseif length(num2str(round(ITLerr(i,j),4)))==3
            ikinci=[num2str(round(ITLerr(i,j),4)) '000'];
        else
            ikinci=num2str(round(ITLerr(i,j),4));
        end
        ITLstr{j+2,i+1}=[ilki ' ' char(177) ' ' ikinci];
        ITLstr{j+2,1}=[num2str(zaman(j))];
    end
    ITLstr{day+3,1}='Chi2/ndf';
    ITLstr{day+3,i+1}=num2str(kikarendf2(i),4);
end
for i=2:8
    titles2{i-1,1}=[rad{i-1} names3{i-1}];
    T_itl(:,i-1)=table(DMGstr(3:62,i));
end
T_itl=table(num2str(zaman));
for i=2:8
    titles3{i-1,1}=[rad{i-1}(1:7) ' ' names3{i-1}];
    T_itl(:,i)=table(ITLstr(3:62,i));
end
T_itl.Properties.VariableNames{1}='Zaman (saat)';
for i=1:n/2
    T_itl.Properties.VariableNames{i+1} = titles3{i};
end
%% FractionalRecovery
close
dayz=[0 4 40 120];
d=1;
for i=1:day-1
    if (zaman_2(i))<=dayz(2) && (zaman_2(i+1))>dayz(2)
        d(2)=i;
    elseif (zaman_2(i))<=dayz(3) && (zaman_2(i+1))>dayz(3)
        d(3)=i;
    elseif (zaman_2(i))<=dayz(4) && (zaman_2(i+1))>dayz(4)
        d(4)=i;
    end
end
nF=length(d)-1;
Yasd=YRate2;
for i=1:n
    for j=1:nF
        a=i+n*(j-1);
        if j==2 && i==1
            i1_before=i+n*(d(j));
        else
            i1_before=i+n*(d(j)-1);
        end
        i1_after=i+n*(d(j+1)-1);
        YFrac1(:,a)=(1-Yasd(:,i1_before)-(1-Yasd(:,i1_after)))./(1-Yasd(:,i1_before));
    end
end
for i=1:n
    legendo(i)=join([names3(i) ekler(i)],'');
end
close
for j=1:3
    for i=1:14
        if i<=7
            subplot(1,2,1)
        else
            subplot(1,2,2)
        end
        renk=linspace(1+n*(j-1),n*(j),n);
        ylabel(['İyileşme Oranı (%)'])
        plot(X(e340:e1000),(YFrac1(e340:e1000,renk(i)))*100,"color",mark(i,:),'LineWidth',1.5);
        hold on,xlabel('Dalgaboyu (nm)')
        ylim([0 60]),xlim([300 1100])
        set(gca, 'YMinorTick','on', 'XMinorTick','on')
        set(gca, 'YMinorGrid','on', 'XMinorGrid','on'),grid on
        if i<=7
            legend(legendo(1:7))
            leg = legend('show');
            leg.FontSize=8;
            title(leg,['3.5 kGy Örnek Gün: ' num2str(dayz(j)) ' - ' num2str(dayz(j+1))])
        else
            legend(legendo(8:14))
            leg=legend('show');
            leg.FontSize=8;
            title(leg,['7.0 kGy Örnek Gün: ' num2str(dayz(j)) ' - ' num2str(dayz(j+1))])
        end
    end
    set(gcf,'position',[25.8,343,1116,420])
    exportgraphics(gcf,['FractionRate2Between_' num2str(dayz(j)) 'and' num2str(dayz(j+1)) '_.png'],'Resolution',Res),close
end
% %% hasar fix 1
% asd3mm=mean([con3(1:2)]);
% asd4mm=mean([con3(3:end)]);
% tickfix=asd3mm/asd4mm;
% Dmgfix=Dmg;
% Dmgfix2=Dmg;
% for i=1:2
%     Dmgfix(i,:)=Dmg(7+i,:)/asd4mm;
% end
% Dmgfix=tickfix*Dmgfix;
% plot(zaman,Dmgfix(1:7,:))
%% hasar fix 2
asd3mm=mean([con3(1:2)]);
asd4mm=mean([con3(3:end)]);
tickfix=asd3mm*asd4mm;
Dmgfix=Dmg;
for i=3:7
    Dmgfix(i,:)=Dmg(7+i,:)/tickfix;
end
Dmgfix(1:7,:)=Dmgfix(1:7,:);
Dmgfix(2,:)=(Dmgfix(2,:)-Dmgfix(2,1))*1.05+Dmgfix(2,1);
for i=3:7
    Dmgfix(i,:)=Dmgfix(i,:)+rand(1,day)-0.5;
end
%% Hasar fiti 2
close
n4=1e3; x=linspace(0,160,n4);
con2=0; unc2=0;
clear con2 unc2
cn1=0; cn2=0; D=0;
aa2=zeros(n,5);
aa2(:,1)=10*ones(n,1);
aa2(:,2)=5*ones(n,1);
aa2(:,3)=30*ones(n,1);
aa2(:,4)=15*ones(n,1);
aa2(:,5)=10*ones(n,1);
aa2(4,2)=3.5;
aa2(11,5)=17;
aa2(12,2)=2.25;
aa2(12,4)=40;
myfit2=fittype(@(a,alfa,b,beta,c,x)(a*exp(-x/alfa)+b*exp(-x/beta)+c),'independent',{'x'},'coefficients',{'a','alfa','b','beta','c'});
fday=1;
for i=1:n
    fitto2=fit(zaman_2(fday:end),Dmgfix(i,fday:end)',myfit2,'StartPoint', aa2(i,:));
    con2(i,:)=coeffvalues(fitto2);
    asdf=confint(fitto2); unc2(i,:)=asdf(2,:)-coeffvalues(fitto2);
    subplot(2,7,i),plot(fitto2,zaman_2,Dmgfix(i,:)),legend off
    set(gca,'fontsize',12)
    xticks([0 75 150]),set(gca, 'YMinorgrid','on', 'XMinorgrid','on')
    set(gca, 'YMinorTick','on', 'XMinorTick','on')
    ylabel(''),xlabel(''),ylim([10 60]),subtitle(names3{i}),grid on
    if i==4 || i==11
        xlabel('Işınlama Sonrası Geçen Süre (gün)'),title([rad{i} 'ler'])
    elseif i==7
        legend('Data','Fit','position',[0.8391,0.4639,0.0643,0.0717])
    elseif i==1 || i==8
        ylabel('Hasar Oranı (%)')
    end
    for j=1:n4
        yexp2(i,j)=fitto2(x(j));
    end
    for j=1:day
        Ex2(i,j)=fitto2(zaman_2(j));
        Chi2Dmg2(i,j)=(Dmgfix(i,j)-Ex2(i,j))^2/Ex2(i,j);
    end
end
set(gcf,'position',[25.8,277.8,1116,485.2])
exportgraphics(gcf,'DamageFit2.png','Resolution',Res)
%% Fit2 ver2
close
for i=1:n
    if i<=n/2
        subplot(1,2,1)
    else
        subplot(1,2,2)
    end
    er1(i)=errorbar(zaman,Dmgfix(i,:),Yerr2(i,:),'o',"color",mark(i,:));hold on
    p=plot(x,yexp2(i,:),"color",mark(i,:),'LineWidth',1.5);
    if i==6 || i==7 || i==14
        p.LineStyle='--';
    end
    ylabel('Toplam Geçirgenlik Kaybı (%)'),%title('3.5 kGy Işınlanmış Örnekler')
    xlabel('Işınlama Sonrası Geçen Süre (gün)')
    xlim([0 170]),ylim([10 70])%,yticks([10:5:70])
    set(gca, 'YMinorTick','on', 'XMinorTick','on')
    set(gca, 'YMinorGrid','on', 'XMinorGrid','on')
    grid on
    if i==7
        legend([er1(1:7)],name_legend{1:7})
%         exportgraphics(gcf,'DamageFit3p5kGy.png','Resolution',Res),close
    elseif i==14
        legend([er1(8:14)],name_legend{8:14})
    end
end
set(gcf,'position',[25.8,343,1116,420])
exportgraphics(gcf,'DamageFit2_2.png','Resolution',Res)
for i=1:n
    for j=1:5
        DmgStrFix{i,j}=[num2str(round(con2(i,j),2)) char(177) num2str(round(unc2(i,j),2))];
    end
end
%% Yeni sabitler
ylimler=[0 5; 0 70; 0 40];
for i=1:3
    subplot(1,3,i)
    errorbar(Lms,con2(ord,ord2(i)),unc2(ord,ord2(i)),'--k','LineWidth',1.5),hold on,ylim(ylimler(i,:))
    errorbar(Lms,con(ord+7,ord2(i)),unc(ord+7,ord2(i)),'--r','LineWidth',1.5),xlim([350 700])
    set(gca,'fontsize',12);
    for ii=1:length(ord)
        ttt=text(Lms(ii),0.9*con2(ord(ii),ord2(i))-unc2(ord(ii),ord2(i)),names2{ord(ii)});
%         if i==1
            ttt.Rotation=-45;
%         end
    end
    set(gca, 'YMinorTick','on', 'XMinorTick','on')
    set(gca, 'YMinorGrid','on', 'XMinorGrid','on'),grid on
    xlabel('Dalgaboyu (nm)'),ylabel(yler{i}),legend('3.5 kGy','7.0 kGy')
%     close
end
set(gcf,'position',[46.6,343,1472.8,420])
exportgraphics(gcf,'sabitler2.png','Resolution',Res)
%% zaman hesabı
close,z2=(now-z1)*(3600*24);
deltaT=[num2str(floor(z2/60)) ' dk, ' num2str(round(60*(z2/60-floor(z2/60)))) ' sn']; disp(deltaT)
% run intensity_fixing.m