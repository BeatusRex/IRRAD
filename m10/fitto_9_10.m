run analiz_m9_m10
%%
clc
for i=1:nLED
    color_o{i}=join([LED{i} ' old']);
    color_n{i}=join([LED{i} ' new']);
    f_new=fit(Res_new(i,:)',int_new(i,:)','poly1');
    ab_n(i,:)=coeffvalues(f_new);
    Rxn(:,i)=linspace(Res_new(i,1)*.9,Res_new(i,end)*1.1,N);
    intfn(:,i)=ab_n(i,1).*Rxn(:,i)+ab_n(i,2);
    
    
    f_old=fit(Res_old(i,:)',int_old(i,:)','poly1');
    ab_o(i,:)=coeffvalues(f_old);
    Rxo(:,i)=linspace(Res_old(i,1)*.9,Res_old(i,end)*1.1,N);
    intfo(:,i)=ab_o(i,1).*Rxo(:,i)+ab_o(i,2);
    
    hold on
    plot(Rxn(:,i),intfn(:,i),join([char(color2{i})],''),'linewidth',2)
end
for i=1:nLED
    plot(Rxo(:,i),intfo(:,i),join(['--',char(color2{i})],''),'linewidth',2)
end
intt=7.9e5;
plot(linspace(20,900,N),ones(N,1)*intt,'k','linewidth',2)
%%
for i=1:nLED
    plot(Res_new(i,:),int_new(i,:),'--+k')
    plot(Res_old(i,:),int_old(i,:),'--+k')
end,set(gca,'fontsize',20)
xlabel('Resistance [\Omega]'),ylabel('Integral of LED Intensity')
legend([color_n color_o {'Wanted Intensity'}]),grid on
%%
%%
for i=1:nLED
    Rneedn(i,1)=(intt-ab_n(i,2))/ab_n(i,1);
    Rneedo(i,1)=(intt-ab_o(i,2))/ab_o(i,1);
end
T=table(char(LED),ab_n,ab_o,[Rneedn,Rneedo], 'VariableNames', {'LED color',...
    'a and b for new','a and b for old','R_n and R_o'});

disp 'Resistant-Intensity equation [int=a*R + b]'
disp '=========================================='
disp 'Table of Old and New LED Stations'
disp(T)
% for i=1:nLED
%     tn{i}=join([LED(i)' ' ==> ' num2str(Rneedn(i)) ' \Omega']);
%     to{i}=join([LED(i)' ' ==> ' num2str(Rneedo(i)) ' \Omega']);
% end
% 
% text(100,5.7e5,tn')
run R_tasarim.m