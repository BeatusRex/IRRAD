%% 5 volt için led parlaklıkları
clear
link='C:\Users\gantu\OneDrive\Belgeler\MATLAB\Master\proje\LED\m19\dot\';
ngroup=25; LED={'UV';'Blue';'White';'Red';'Green'};LED=sort(LED); nLED=length(LED);
S={'Old','New'};S=sort(S); nS=length(S);
matfiles=dir(fullfile(link, '*.txt'));
nfiles=length(matfiles);
for i=1:nfiles
    names=matfiles(i).name;
    D1=load([link names]);
    X(:,1)=D1(:,1);
    Yn(:,i)=D1(:,2);
end
%%
for i=1:ngroup
    Yavg(i,:)=mean(Yn(:,1+(i-1)*nLED*nS:(i)*nLED*nS)');
end,Yavg=Yavg';
ara=50:350;
n=length(Yavg(:,1));
m=length(Yavg(1,:));
for i=1:m
    Bl(i,1)=mean(Yavg(ara,i));
end
%%
for i=1:n
    for j=1:m
        if Yavg(i,j)>=Bl(j)
            yavg(i,j)=Yavg(i,j)-Bl(j);
        else
            yavg(i,j)=0;
        end
    end
end,intp=zeros(n,m);
for i=1:n-1
    for j=1:m
        intp(i+1,j)=mean([yavg(i,j) yavg(i+1,j)])*(X(i+1)-X(i))+intp(i,j);
    end
end
col={'ob','og','or','om','oy'};
for i=1:nLED
    for j=1:nS
        lgnd(j+nS*(i-1))=join([S(j) LED(i)]);
    end
end
for i=1:nLED
    a1(i)=2*i-1;
    a2(i)=2*i;
end
int=intp(end,a1)';
int(:,2)=intp(end,a2)';
for i=1:5
    hold on
    plot(int(i,:),col{i},'LineWidth',10)
end
x=1:2; xticks(x),xticklabels(S),xlim([-1 4])
legend(LED),grid on
ylabel('Intensity')
