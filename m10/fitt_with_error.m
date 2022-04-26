
for i=1:nLED
    color_o{i}=join([LED{i} ' old']);
    color_n{i}=join([LED{i} ' new']);
    f_new=fit(Res_new(i,:)',int_new(i,:)','poly1');
    ab_n(i,:)=coeffvalues(f_new);
    Rxn(:,i)=linspace(Res_new(i,1),Res_new(i,end),N);
    intfn(:,i)=ab_n(i,1).*Rxn(:,i)+ab_n(i,2);
    errn(:,i)=abs(ab_n(i,1).*Rxn(:,i)*0.05);
    
    
    f_old=fit(Res_old(i,:)',int_old(i,:)','poly1');
    ab_o(i,:)=coeffvalues(f_old);
    Rxo(:,i)=linspace(Res_old(i,1),Res_old(i,end),N);
    intfo(:,i)=ab_o(i,1).*Rxo(:,i)+ab_o(i,2);
    erro(:,i)=abs(ab_o(i,1).*Rxo(:,i)*0.05);
    
    hold on
    errorbar(Rxn(:,i),intfn(:,i),errn(:,i),join([char(color2{i})],''),'CapSize',0)
end
for i=1:nLED
    errorbar(Rxo(:,i),intfo(:,i),erro(:,i),join(['--',char(color2{i})],''),'CapSize',0)
end
intt=7.9e5;
plot(linspace(20,900,N),ones(N,1)*intt,'k','linewidth',2)
%%
for i=1:nLED
    plot(Res_new(i,:),int_new(i,:),'--+k')
    plot(Res_old(i,:),int_old(i,:),'--+k')
end
xlabel('Resistance [\Omega]'),ylabel('Integral of LED Intensity')
legend([color_n color_o {'Wanted Intensity'}]),grid on