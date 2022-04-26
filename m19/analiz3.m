%% 5 volt için led parlaklıkları
clear
link='C:\Users\gantu\OneDrive\Belgeler\MATLAB\Master\proje\LED\m19\dot\';
LED={'UV';'Blue';'White';'Green';'Red'}; nLED=length(LED);
col={'m','b','c','g','r'};
for ii=1:nLED
matfiles=dir(fullfile(link, [LED{ii} '*new*.txt']));
nfiles=length(matfiles);
for i=1:nfiles
    names=matfiles(i).name;
    D1=load([link names]);
    X(:,1)=D1(:,1);
    Yn(:,i+(ii-1)*nfiles)=D1(:,2);
end
end
n=length(Yn(:,1));
m=length(Yn(1,:));
%%
ara=50:350;
for i=1:m
    Bl(i,1)=mean(Yn(ara,i));
end
%%
for i=1:n
    for j=1:m
        if Yn(i,j)>=Bl(j)
            Y_g(i,j)=Yn(i,j)-Bl(j);
        else
            Y_g(i,j)=0;
        end
    end
end,int=zeros(n,m);
for i=1:n-1
    for j=1:m
        int(i+1,j)=mean([Y_g(i,j) Y_g(i+1,j)])*(X(i+1)-X(i))+int(i,j);
    end
end
%%
int1=int(end,:)';
for i=1:nLED
    ara=1+(i-1)*nfiles:(i)*nfiles;
    int_avg(i,1)=mean(int1(ara));
    int_avg(i,2)=std(int1(ara))*10;
end
errorbar([1:nLED],int_avg(:,1),int_avg(:,2),'k:.',"LineWidth",2),ylim([3e5 10e5]),xlim([0 6])
xticks([1:5]),xticklabels(LED),set(gca,"fontsize",24),grid on
ylabel("Intensity of LEDs"),xlabel("Type of LED")
%%
for i=1:nLED
    Y_m(:,i)=mean(Y_g(:,1+(i-1)*nfiles:i*nfiles)')';
end
for i=1:nLED
    plot(X,Y_m(:,i),col{i},"LineWidth",2),hold on
end
legend(LED),set(gca,"Fontsize",24),xlabel("Wavelength [nm]")
xlim([350 800]),ylim([10 1e5]),set(gca,'yscale','log'),grid on
ylabel("Intensith [Photons per 100ms]")