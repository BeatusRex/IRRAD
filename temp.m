%% Karacisim ışıması grafikleri için kod
clear
h=6.6234e-34; %J.s
c=2.998e8; %m/s
k=1.3807e-23; %J/K
N=1e4;
M=11;
L=linspace(200,100000,N)'*1e-9;
v=c./L;
T=linspace(300,5000,M)';
for i=1:N
    for j=1:M
        B(i,j)=8*pi*h*v(i)^2/(c^2*(exp(h*v(i)/(k*T(j)))-1));
    end
end
% plot(L,B(:,:))
jett=jet;
subplot(1,2,1)
for i=1:M
    renk(i,:)=jett(ceil(length(jett)*(i-1)/M)+1,:);
    plot3(L*1e9,T(i).*ones(N,1),B(:,i)*1e12,'Color',[renk(i,:)]),hold on
end
grid on
set(gca, 'YMinorTick','on', 'XMinorTick','on')
set(gca, 'YMinorGrid','on', 'XMinorGrid','on')
set(gca, 'ZMinorTick','on', 'ZMinorGrid','on')
xlabel('Dalgaboyu (nm)'),ylabel('Sıcaklık (K)'),zlabel({'Katıaçı Başı'; 'Güç (W / \Omega)'})
xlim([200 8000]),view([-25,25])
subplot(1,2,2)
m=1;
plot(L*1e9,B(:,m)*1e12,'Color',[renk(m,:)])
grid on
set(gca, 'YMinorTick','on', 'XMinorTick','on')
set(gca, 'YMinorGrid','on', 'XMinorGrid','on')
xlabel('Dalgaboyu (nm)'),ylabel('Katıaçı Başı Güç (W /\Omega)'),title(['T = ' num2str(round(T(m))) ' K'])
set(gcf,'position',[25.8,343,1116,420])
exportgraphics(gcf,'TempFrequency.png','Resolution',600)