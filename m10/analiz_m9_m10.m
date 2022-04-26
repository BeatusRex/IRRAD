clear
link_m10='C:\Users\gantu\OneDrive\Belgeler\MATLAB\Master\proje\LED\m10\';
link_m9='C:\Users\gantu\OneDrive\Belgeler\MATLAB\Master\proje\LED\m9\';
ngroup=25;LED={'UV';'Blue';'White';'Red';'Green'};LED=sort(LED); nLED=length(LED);
num=1:5; nnum=length(num);
matfiles_new=dir(fullfile(link_m10, '*.txt'));
matfiles_old=dir(fullfile(link_m9, '*.txt'));
nfiles=length(matfiles_new);
%%
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
%%
for i=1:ngroup
    for j=1:nnum
        for k=1:nLED
            names=matfiles_new(i+ngroup*(j-1)+nnum*ngroup*(k-1)).name;
            D=load([link_m10 names]);
            names=matfiles_old(i+ngroup*(j-1)+nnum*ngroup*(k-1)).name;
            D2=load([link_m9 names]);
            X(:,1)=D(:,1); Y_new(:,i,j,k)=D(:,2); Y_old(:,i,j,k)=D2(:,2);
        end
    end
end
N=length(X);
%%
% B(:,:,:)=Y_new(:,:,:,1);
% G(:,:,:)=Y_new(:,:,:,2);
% R(:,:,:)=Y_new(:,:,:,3);
% U(:,:,:)=Y_new(:,:,:,4);
% W(:,:,:)=Y_new(:,:,:,5);
for i=1:N
    for j=1:nnum 
        for k=1:nLED
            Y_avgnew(i,j+nnum*(k-1))=mean(Y_new(i,:,j,k));
            Y_avgold(i,j+nnum*(k-1))=mean(Y_old(i,:,j,k));
        end
    end
end
%%
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
%%
for i=1:nLED
    for j=1:nnum
        int_new(i,j)=sum(yavg_new(:,j+nnum*(i-1)))*(X(end)-X(1))/N;
        int_old(i,j)=sum(yavg_old(:,j+nnum*(i-1)))*(X(end)-X(1))/N;
    end
end
color={'+b';'+g';'+r';'+m';'+c'};
color2={'b';'g';'r';'m';'c'}; v2=[4 2 4 3 2];