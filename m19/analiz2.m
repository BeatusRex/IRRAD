%% 5 volt için led parlaklıkları
clear
link='C:\Users\gantu\OneDrive\Belgeler\MATLAB\Master\proje\LED\m19\dot\';
ngroup=25; LED={'UV';'Blue';'White';'Red';'Green'};LED=sort(LED); nLED=length(LED);
S={'Old','New'};S=sort(S); nS=length(S);
S2={'No. 1','No. 2'};
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
col={'ob','og','or','om','oc'};
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
err1=(max(max(int))/min(min(int))-1)*100;
intm=mean(mean(int));
asd=abs(int/intm-1)*100;
err2=asd(:,1); err2(6:10,1)=asd(:,2);
xx(1:5)=-.3;xx(6:10)=2;xx=xx+0.1;
yy(:,1)=int(:,1);yy(6:10,1)=int(:,2);
%%
for i=1:5
    hold on
    plot(int(i,:),col{i},'LineWidth',10)
end
for i=1:10
    txt{i}=['Error of LED: ' num2str(err2(i),2) '%'];
end
set(gca,'fontsize',22)
text(xx,yy,txt,'fontsize',22)
subtitle(['Error Btwn max and min Val.: ', num2str(err1,2) '%'])
x=1:2; xticks(x),xticklabels(S2),xlim([-1 4]),xlabel('LED station')
legend(LED),ylabel('Integral of LED Intensity'),grid on