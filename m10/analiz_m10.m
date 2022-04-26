clear
target_link='C:\Users\gantu\OneDrive\Belgeler\MATLAB\Master\proje\LED\m10\';
ngroup=25;LED={'UV';'Blue';'White';'Red';'Green'};LED=sort(LED); nLED=length(LED);
num=1:5; nnum=length(num);
matfiles=dir(fullfile(target_link, '*.txt'));
nfiles=length(matfiles);
%%
Res=[148.4 170.6 191.6 211.3 238.2;...
    595.1 653.7 708 794 873;...
    27.9 36.4 43.8 54.8 66.8;...
    486.2 526.6 559.2 604.5 643;...
    315.9 338.4 377.3 419.2 467.8;];
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
color2={'b';'g';'r';'m';'c'}; v2=[2 2 3 2 2];
for i=1:nLED
    v(i)=v2(i)+nLED*(i-1);
end
% for i=1:nLED
%     hold on
%     plot(Res(i,:),int(i,:),join(['--' char(color{i})],''))
% end
%%
yy=3e2*ones(nLED,1); xx=[450 500 600 375 525];
yy(2)=yy(2)*2;
for i=1:nLED
    hold on
    plot(X,yavg(:,v(i)),join([ char(color2{i})],''),'linewidth',2)
    xlabel('Wavelength [nm]'),xlim([350 750]),ylabel('Intensity')
    set(gca,'yscale','log')
    txt = [LED{i} ' ==> ' char(num2str(int(i,v2(i))))];
    text(xx(i),yy(i),txt)
end
ylim([3e1 1e5])
legend(LED),grid on
