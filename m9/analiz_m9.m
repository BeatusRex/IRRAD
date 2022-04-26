clear
target_link='C:\Users\gantu\OneDrive\Belgeler\MATLAB\Master\proje\LED\m9\';
ngroup=25;LED={'UV';'Blue';'White';'Red';'Green'};LED=sort(LED); nLED=length(LED);
num=1:5; nnum=length(num);
matfiles=dir(fullfile(target_link, '*.txt'));
nfiles=length(matfiles);
%%
Res=[518.3,543.2,568.5,604.7,643.6;...
    622.5,670,734,822,904;...
    130,139.6,155.8,171.5,184.7;...
    599.3,645.7,689,746,787;...
    613.6,637.4,673,714,755];
Res(1,:)=Res(1,:)-300;
Res(5,:)=Res(5,:)-300;
%%
for i=1:ngroup
    for j=1:nnum
        for k=1:nLED
            names=matfiles(i+ngroup*(j-1)+nnum*ngroup*(k-1)).name;
            D=load([target_link names]);
            X(:,1)=D(:,1); Y(:,i,j,k)=D(:,2);
        end
    end
end
N=length(X);
%%
B(:,:,:)=Y(:,:,:,1);
G(:,:,:)=Y(:,:,:,2);
R(:,:,:)=Y(:,:,:,3);
U(:,:,:)=Y(:,:,:,4);
W(:,:,:)=Y(:,:,:,5);
for i=1:N
    for j=1:nnum 
        for k=1:nLED
            Y_avg(i,j+nnum*(k-1))=mean(Y(i,:,j,k));
        end
    end
end
%%
for i=1:nnum*nLED
    Base(i,1)=mean(Y_avg(100:350,i));
end
for i=1:N
    for j=1:nnum*nLED
        if Y_avg(i,j)>=Base(j)
            yavg(i,j)=Y_avg(i,j)-Base(j);
        else
            yavg(i,j)=0;
        end
    end
end
%%
for i=1:nLED
    for j=1:nnum
        int(i,j)=sum(yavg(:,j+nnum*(i-1)))*(X(end)-X(1))/N;
        intful(i,j)=sum(Y_avg(:,j+nnum*(i-1)))*(X(end)-X(1))/N;
    end
end
color={'+b';'+g';'+r';'+m';'+c'};
color2={'b';'g';'r';'m';'c'}; v2=[4 2 4 3 2];
for i=1:nLED
    v(i)=v2(i)+nLED*(i-1);
end
% for i=1:nLED
%     hold on
%     plot(Res(i,:),int(i,:),join(['--' char(color{i})],''))
% end
%%
pos=[415,34745; 500,34745; 595,46085;...
    365,70405; 520,6250;];
yy=pos(:,2); xx=pos(:,1);
for i=1:nLED
    hold on
    plot(X,yavg(:,v(i)),join([ char(color2{i})],''),'linewidth',2)
    xlabel('Wavelength [nm]'),xlim([350 750]),ylabel('Intensity [Photons per 100ms]')
    set(gca,'yscale','log')
    txt = [LED{i} ' ==> ' char(num2str(int(i,v2(i)),7))];
    text(xx(i),yy(i),txt,'fontsize',18)
end
set(gca,'fontsize',22)
ylim([3e1 1e5])
legend(LED),grid on