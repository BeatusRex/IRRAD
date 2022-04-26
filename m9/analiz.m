clear
target_link='C:\Users\gantu\OneDrive\Belgeler\MATLAB\Master\proje\LED\m9\';
ngroup=25;LED={'UV';'Blue';'White';'Red'};LED=sort(LED); nLED=length(LED);
num=1:5; nnum=length(num);
matfiles=dir(fullfile(target_link, '*.txt'));
nfiles=length(matfiles);
nLED=length(LED); N=2048;
x=0;
%%
Res=[518.3,543.2,568.5,604.7,643.6;...
    130,139.6,155.8,171.5,184.7;...
    599.3,645.7,689,746,787;...
    613.6,637.4,673,714,755];
Res(1,:)=Res(1,:)-300;
Res(4,:)=Res(4,:)-300;
%%
x=0;
for i=1:ngroup
    for j=1:nnum
        for k=1:nLED
            names=matfiles(i+nnum*(j-1)+nLED*(k-1)).name;
            D=load([target_link names]);
            X(:,1)=D(:,1); Y(:,i,j,k)=D(:,2);
        end
    end
end
%%
B(:,:,:)=Y(:,:,:,1);
R(:,:,:)=Y(:,:,:,2);
U(:,:,:)=Y(:,:,:,3);
W(:,:,:)=Y(:,:,:,4);
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
color={'+b';'+r';'+m';'+c'};
for i=1:nLED
    hold on
    plot(Res(i,:),intful(i,:),join(['--' char(color{i})],''))
end
legend(LED)
%%