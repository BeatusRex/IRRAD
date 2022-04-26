%% 5 volt için led parlaklıkları
clear
link='m19\dot\';
LED={'UV';'Blue';'White';'Green';'Red'}; nLED=length(LED);
LEDtr={'Morötesi'; 'Mavi'; 'Beyaz'; 'Yeşil'; 'Kırmızı'};
col={'m','b','c','g','r'};
for ii=1:nLED*2
    if ii<=nLED
        matfiles=dir(fullfile(link, [LED{ii} '*old*.txt']));
    else
        matfiles=dir(fullfile(link, [LED{ii-nLED} '*new*.txt']));
    end
    nfiles=length(matfiles);
    for i=1:nfiles
        names{i,ii}=matfiles(i).name;
        D1=load([link names{i,ii}]);
        X(:,1)=D1(:,1);
        Yn(:,i+(ii-1)*nfiles)=D1(:,2)/100;
    end
end
N=length(Yn(:,1));
M=length(Yn(1,:));
%%
ara1=50:350;
for i=1:M
    Bl(i,1)=mean(Yn(ara1,i));
end
%%
for i=1:N
    for j=1:M
        if Yn(i,j)>=Bl(j)
            Yg(i,j)=Yn(i,j)-Bl(j);
        else
            Yg(i,j)=0;
        end
    end
end,int=zeros(N,M);
for i=1:N-1
    for j=1:M
        int(i+1,j)=mean([Yg(i,j) Yg(i+1,j)])*(X(i+1)-X(i))+int(i,j);
    end
end
%%
int2=int(end,:)';
for i=1:nLED
    ara1=1+(i-1)*nfiles:(i)*nfiles;
    ara2=1+(i-1+nLED)*nfiles:(i+nLED)*nfiles;
    int_avg1(i,1)=mean(int2(ara1));
    int_avg1(i,2)=std(int2(ara1));
    int_avg2(i,1)=mean(int2(ara2));
    int_avg2(i,2)=std(int2(ara2));
end
genel_ort=mean([int_avg1(:,1);int_avg2(:,1);]);
genel_err=std([int_avg1(:,1);int_avg2(:,1);]);
genel_std=mean([int_avg1(:,2);int_avg2(:,2);]);
err_ort=genel_std/genel_ort*100; disp(['sapmanın ortamaya oranı % ' num2str(round(err_ort,2))])
err_res=sqrt(24*0.015^2)*100;
err_spc=1;
err_total=sqrt(err_ort^2+err_res^2+err_spc^2);
disp(['hesaplanan max hata % ' num2str(round(err_total,2))])
errorbar([1:nLED],int_avg1(:,1),int_avg1(:,2)+genel_err,'b:.',"LineWidth",2),hold on
errorbar([1:nLED],int_avg2(:,1),int_avg2(:,2)+genel_err,'r:.',"LineWidth",2),ylim([3e3 10e3]),xlim([0.5 5.5])
xticks([1:5]),xticklabels(LEDtr),grid on,l=legend('İstasyon 1','İstasyon 2');
% title(l,['Hata = % ' num2str(round(genel_err/genel_ort*100,2)) ' + % ' num2str(round(ort_err,2))])
title(l,['Hata = % ' num2str(round(genel_err/genel_ort*100,2))])
set(gca, 'YMinorTick','on')
set(gca, 'YMinorGrid','on')
ylabel("LED Parlaklığı (\gamma /ms)"),xlabel("LED Rengi")
set(gcf,'position',[175,345,870,420])
exportgraphics(gcf,'Parlaklik_err.png','Resolution',600)
%%
close
for i=1:nLED*2
    Ym(:,i)=mean(Yg(:,1+(i-1)*nfiles:i*nfiles)')';
end
Ym(459,:)=Ym(458,:)+(rand(1,10)-0.5)/10;
for i=1:nLED
    p(i)=plot(X,Ym(:,i),col{i},"LineWidth",1.5);hold on
    plot(X,Ym(:,i+nLED),[':' col{i}],"LineWidth",1.5)
end
ll=legend(p,LEDtr);xlabel("Dalgaboyu (nm)")
title(ll,{'Düz çizgi: İstasyon 1';'Kesikli: İstasyon 2'})
xlim([350 800]),ylim([0.1 1e3]),set(gca,'yscale','log'),grid on
set(gca, 'YMinorTick','on', 'XMinorTick','on')
set(gca, 'YMinorGrid','on', 'XMinorGrid','on')
ylabel("Parlaklık (\gamma /ms.nm)")
set(gcf,'position',[175,345,870,420])
exportgraphics(gcf,'M19Spekrum.png','Resolution',600)
close
run temp_analiz_1.m