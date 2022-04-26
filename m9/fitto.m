run analiz_m9
%%
for i=1:nLED
    f=fit(Res(i,:)',int(i,:)','poly1');
    ab(i,:)=coeffvalues(f);
    Rx(:,i)=linspace(Res(i,1)*.9,Res(i,end)*1.1,N);
    intf(:,i)=ab(i,1).*Rx(:,i)+ab(i,2); 
    hold on
    plot(Rx(:,i),intf(:,i),join([char(color2{i})],''),'linewidth',2)
end
intt=7.9e5;
plot(linspace(100,900,N),ones(N,1)*intt,'linewidth',2)
for i=1:nLED
    plot(Res(i,:),int(i,:),'--+k')
end
xlabel('Resistance [\Omega]'),ylabel('Integral of LED Intensity')
legend([LED' {'Wanted Intensity'}]),grid on
c=ones(nLED,1)*1.005; %c(1)=1/c(1);
d=ones(nLED,1); %d(1)=.8*d(1);
% title('a x + b')
for i=1:nLED
    Rneed(i,1)=(intt-ab(i,2))/ab(i,1);
    txt = [LED{i} ' => ' char(num2str(Rneed(i))) ' ' '\Omega'];
    txt2=[num2str(round(ab(i,1))) ' R + ' num2str(round(ab(i,2)))];
    y=intt*c(i);e=d(i)*Rneed(i);
    text(e,y,txt)
    text(e+5,y*.99,txt2)
end,xlim([100 900])