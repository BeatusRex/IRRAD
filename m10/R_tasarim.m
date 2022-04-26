Rsetn=[150 20 0;560 100 3.9;39 2.2 0;390 120 0;220 68 68];
Rseto=[200 33 0;470 180 47;100 47 0;470 150 0;220 120 3.3];
Rn=sum(Rsetn')';
Ro=sum(Rseto')';
for i=1:nLED
    yeni(i,1)=ab_n(i,1)*Rn(i)+ab_n(i,2);
    eski(i,1)=ab_o(i,1)*Ro(i)+ab_o(i,2);
    ery(i,1)=abs(yeni(i)-intt)/intt*100;
    ere(i,1)=abs(eski(i)-intt)/intt*100;
end
disp(['Wanted int ==> ' num2str(intt)])
disp(' ')
T2=table(char(LED),[Rn Ro],[yeni eski],[ery ere],'VariableNames',{'LED color','Applied R(N&O)','Calculated Int(N&O)','Relative %err(N&O)'});
disp(T2)