%% Sıcaklık kontrolü bulunmayan normal oda
clear
z1=now;
link='Pulsed\deney2\';
Tlog=importdata('Pulsed\deney2\TEMP\temp.txt');
A=dir(fullfile(link, '*.txt'));
nfile=length(A);
nlog=length(Tlog);
for i=1:nfile
    D=importdata([link A(i).name]);
    if i==1
        X=D(:,1);
    end
    Y(:,i)=D(:,2);
    zaman(i,1)=round((A(i).datenum-A(1).datenum)*24*60);
end
Res=600;
%% Özel konumları bulma
nx=length(X);
for i=1:nx
    if X(i)<=200
        i200=i;
    elseif X(i)<=350
        i350=i;
    elseif X(i)<=450
        i450=i;
    elseif X(i)<=1000
        i1000=i;
    end
end
%% Temp datası okuma
it=0; ip=0; ih=0;
for i=1:nlog
    J=length(Tlog{i});
    for j=1:J
        if Tlog{i}(j)=='>'
            i1=j+2;
        end
        if Tlog{i}(j)=='='
            i2=j+2;
            cc=1;
        else
            cc=0;
        end
        if cc==1
            if Tlog{i}(i1:i1+2)=='Tem'
                it=it+1;
                T(it,1)=str2num(Tlog{i}(i2:i2+4));
            elseif Tlog{i}(i1:i1+2)=='Pre'
                ip=ip+1;
                P(ip,1)=str2num(Tlog{i}(i2:i2+5));
            elseif Tlog{i}(i1:i1+2)=='Hum'
                ih=ih+1;
                H(ih,1)=str2num(Tlog{i}(i2:i2+4));
            end
        end
    end
                
end
%% BG max noktası ve BG yok etme
BGL=max(max(Y([i200:i350, i450:i1000])));
for i=1:nx
    for j=1:nfile
        if Y(i,j)>=BGL
            Y2(i,j)=Y(i,j)-BGL;
        else
            Y2(i,j)=0;
        end
    end
end
%% Max konumları ve değerleri bulmak
% Gx=linspace(350,450,1001);
% for i=1:nfile
%     f=fit(X(i350:i450),Y2(i350:i450,i),'gauss2');
%     cons(i,:)=coeffvalues(f);
%     fmax=max(f(Gx));
%     for j=1:length(Gx)
%         if f(Gx(j))==fmax
%             Xmax(i,1)=Gx(j);
%         end
%     end
% end
% =========================================================================
for i=1:nfile
    fmax=max(max(Y2));
    for j=1:nx
        if Y2(j,i)==max(Y2(:,i))
            Y2max(i,1)=Y2(j,i);
            Xmax(i,1)=X(j);
        end
    end
end
%% Sıcamlık ve Nem
close
yyaxis left
plot(zaman,T(1:nfile),'LineWidth',2)
ylim([20 30]),ylabel('Sıcaklık (^{o}C)')
yyaxis right
plot(zaman,H(1:nfile),'LineWidth',2)
ylim([10 45]),ylabel('Bağıl Nem (%)')
xlabel('Zaman (dk)')
set(gca, 'YMinorTick','on', 'XMinorTick','on')
set(gca, 'YMinorGrid','on', 'XMinorGrid','on')
grid on
set(gcf,'position',[175,345,870,420]) % tek grafik boyutu
exportgraphics(gcf,'TvsH_1.png','Resolution',Res)
%% Max çizgisi 3D
close
subplot(1,2,1)
surf(zaman,X',Y2/100,'EdgeColor','none'),hold on,%view(2)
colormap jet; %colorbar
ylim([mean(Xmax)-20 mean(Xmax)+20]),xlim([min(zaman),max(zaman)])
p2=plot3(zaman,Xmax,Y2max/100,'k','LineWidth',2);
set(gca, 'YMinorTick','on', 'XMinorTick','on')
set(gca, 'YMinorGrid','on', 'XMinorGrid','on')
set(gca, 'ZMinorTick','on', 'ZMinorGrid','on')
view([40 45])
% legend(p2,'Azami değer konumu','Location','best')
ylabel('Dalgaboyu (nm)'),xlabel('Zaman (dk)')
zlabel('Foton Sayımı (\gamma /ms.nm)')
% exportgraphics(gcf,'MaxLine_3D_1.png','Resolution',Res)
% %% Max çizgisi 2D
% close
subplot(1,2,2)
surf(zaman,X',Y2/100,'EdgeColor','none');hold on,view(2)
colormap jet; c1=colorbar;
ylim([mean(Xmax)-20 mean(Xmax)+20]),xlim([min(zaman),max(zaman)])
p2=plot3(zaman,Xmax,Y2max/100,'k','LineWidth',2);
legend(p2,'Azami değer konumu')
ylabel('Dalgaboyu (nm)'),xlabel('Zaman (dk)')
ylabel(c1,'Foton Sayımı (\gamma /ms.nm)')
set(gcf,'position',[25.8,343,1116,420]) % çoklu grafik boyutu
exportgraphics(gcf,'MaxLine1.png','Resolution',Res)
%% Sıcamlık ve Max
close
subplot(1,2,1)
yyaxis left
plot(zaman,T(1:nfile),'LineWidth',2)
ylim([20 30]),ylabel('Sıcaklık (^{o}C)')
yyaxis right
plot(zaman,100*(Y2max./mean(Y2max)),'LineWidth',1.5)
ylabel('Azami Parlaklığın Bağıl Oranı (%)')
xlabel('Zaman (dk)')
grid on
set(gca, 'YMinorTick','on', 'XMinorTick','on')
set(gca, 'YMinorGrid','on', 'XMinorGrid','on')
% exportgraphics(gcf,'TvsY2Max_1.png','Resolution',Res)
% %% integral hesabı
% close
subplot(1,2,2)
Yint=sum(Y2)*(X(2)-X(1))/(100);
yyaxis left
plot(zaman,T(1:nfile),'LineWidth',2)
ylim([20 30]),ylabel('Sıcaklık (^{o}C)')
yyaxis right
plot(zaman,Yint,'LineWidth',1.5)
ylabel('Spektrum İntegrali (\gamma /ms)')
xlabel('Zaman (dk)')
set(gca, 'YMinorTick','on', 'XMinorTick','on')
set(gca, 'YMinorGrid','on', 'XMinorGrid','on')
grid on
set(gcf,'position',[25.8,343,1116,420])
exportgraphics(gcf,'TvsYint_1.png','Resolution',Res)
%% geçen zaman
z2=(now-z1)*(3600*24);
deltaT=[num2str(floor(z2/60)) ' dk, ' num2str(round(60*(z2/60-floor(z2/60)))) ' sn']; disp(deltaT)
close
run temp_analiz_2.m