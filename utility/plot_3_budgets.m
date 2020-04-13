function out=plot_3_budgets(toutmat,CONC,INT,Nyear,valflag,nsmth,xlimsin)
figure;
subplot(3,2,1),...
plot(toutmat./86400./360,smooth(CONC{3}.PROD.*86400,Nyear),...
toutmat./86400./360,smooth(CONC{3}.LATADV.*86400,Nyear),...
toutmat./86400./360,smooth(CONC{3}.LATMIX.*86400,Nyear),...
toutmat./86400./360,smooth(CONC{3}.MIX.*86400,Nyear),...
toutmat./86400./360,smooth(CONC{3}.RATE.*86400,Nyear),'k-'); hold on;
hold on;
title('(A) Two-layer avg. nutrient concentration budget')
xlabel('Year')
xlim(xlimsin{3})
ylabel('mmol/m^3/d')
grid on
set(gcf,'color','w')
set(gca,'fontsize',11)
set(gca,'linewidth',1);

subplot(3,2,2),...

plot(toutmat./86400./360,smooth(INT{3}.PROD.*86400,Nyear),...
toutmat./86400./360,smooth(INT{3}.LATADV.*86400,Nyear),...
toutmat./86400./360,smooth(INT{3}.LATMIX.*86400,Nyear),...
toutmat./86400./360,smooth(INT{3}.MIX.*86400,Nyear),...
toutmat./86400./360,smooth(INT{3}.SUBD.*86400,Nyear),...
toutmat./86400./360,smooth(INT{3}.RATE.*86400,Nyear),'k-'); hold on;
if valflag == 1
plot(toutmat./86400./360,smooth(INT{3}.REAC.*86400,Nyear),'magenta--'); 
end
title('(B) Two-layer integrated nutrient budget')
legend('BIO','LADV','LMIX','VMIX','SUBD','RATE','location','northwest','orientation','horizontal')
xlabel('Year')
xlim(xlimsin{3})
ylabel('mmol/m^2/d')
grid on
set(gca,'fontsize',11)
set(gca,'linewidth',1)

subplot(3,2,3),...
    
plot(toutmat./86400./360,smooth(CONC{1}.PROD.*86400,Nyear),...
toutmat./86400./360,smooth(CONC{1}.LATADV.*86400,Nyear),...
toutmat./86400./360,smooth(CONC{1}.LATMIX.*86400,Nyear),...
toutmat./86400./360,smooth(CONC{1}.MIX.*86400,Nyear),...
toutmat./86400./360,smooth(CONC{1}.ENT.*86400,Nyear),...
toutmat./86400./360,smooth(CONC{1}.RATE.*86400,Nyear),'k-'); hold on;
if valflag == 1
    plot(toutmat./86400./360,smooth(CONC{1}.REAC.*86400,Nyear),'magenta--'); 
end
title('(C) Surface-layer nutrient concentration budget')
xlabel('Year')
xlim(xlimsin{3})
ylabel('mmol/m^3/d')
grid on
set(gca,'fontsize',11)
set(gca,'linewidth',1)

subplot(3,2,4),...
    
plot(toutmat./86400./360,smooth(INT{1}.PROD.*86400,Nyear),...
toutmat./86400./360,smooth(INT{1}.LATADV.*86400,Nyear),...
toutmat./86400./360,smooth(INT{1}.LATMIX.*86400,Nyear),...
toutmat./86400./360,smooth(INT{1}.MIX.*86400,Nyear),...
toutmat./86400./360,smooth(INT{1}.ENT.*86400,Nyear),...
toutmat./86400./360,smooth(INT{1}.RATE.*86400,Nyear),'k-'); hold on;
if valflag == 1
plot(toutmat./86400./360,smooth(INT{1}.REAC.*86400,Nyear),'magenta--'); 
end
title('(D) Depth-integrated surface-layer nutrient budget')
legend('BIO','LADV','LMIX','VMIX','ENT','RATE','location','northwest','orientation','horizontal')
xlabel('Year')
xlim(xlimsin{3})
ylabel('mmol/m^2/d')
grid on
set(gca,'fontsize',11);
set(gca,'linewidth',1);
ylim([-2 2])


subplot(3,2,5),...
plot(toutmat./86400./360,smooth(CONC{2}.PROD.*86400,Nyear),...
toutmat./86400./360,smooth(CONC{2}.LATADV.*86400,Nyear),...
toutmat./86400./360,smooth(CONC{2}.LATMIX.*86400,Nyear),...
toutmat./86400./360,smooth(CONC{2}.MIX.*86400,Nyear),...
toutmat./86400./360,smooth(CONC{2}.ENT.*86400,Nyear),...
toutmat./86400./360,smooth(CONC{2}.RATE.*86400,Nyear),'k-'); hold on;
if valflag == 1
plot(toutmat./86400./360,smooth(CONC{2}.REAC.*86400,Nyear),'magenta--'); 
end
title('(E) Deep-layer nutrient concentration budget')
xlabel('Year')
xlim(xlimsin{3})
ylim([-.015 .015])

ylabel('mmol/m^3/d')
grid on
set(gcf,'color','w')
set(gca,'fontsize',11)
set(gca,'linewidth',1)


subplot(3,2,6),...
plot(toutmat./86400./360,smooth(INT{2}.PROD.*86400,Nyear),...
toutmat./86400./360,smooth(INT{2}.LATADV.*86400,Nyear),...
toutmat./86400./360,smooth(INT{2}.LATMIX.*86400,Nyear),...
toutmat./86400./360,smooth(INT{2}.MIX.*86400,Nyear),...
toutmat./86400./360,smooth(INT{2}.ENT.*86400,Nyear),...
toutmat./86400./360,smooth(INT{2}.RATE.*86400,Nyear),'k-'); hold on;
if valflag == 1
    plot(toutmat./86400./360,smooth(INT{2}.REAC.*86400,Nyear),'magenta--'); 
end
title('(F) Depth-integrated deep-layer nutrient budget')
xlabel('Year')
xlim(xlimsin{3})
ylabel('mmol/m^2/d')
ylim([-2 2])
grid on
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 1.0;
end
set(gca,'fontsize',11)
set(gca,'linewidth',1)

%% first year
for it = 1:2
    if it == 1
        xlims=xlimsin{1}
    elseif it == 2
        xlims=xlimsin{2}
    end
figure;
subplot(3,2,1),...
plot(toutmat./86400./360,smooth(CONC{3}.PROD.*86400,nsmth),...
toutmat./86400./360,smooth(CONC{3}.LATADV.*86400,nsmth),...
toutmat./86400./360,smooth(CONC{3}.LATMIX.*86400,nsmth),...
toutmat./86400./360,smooth(CONC{3}.MIX.*86400,nsmth),...
toutmat./86400./360,smooth(CONC{3}.RATE.*86400,nsmth),'k-'); hold on;
hold on;
title('(A) Two-layer avg. nutrient concentration budget')
xlabel('Year')
xlim(xlims)
ylabel('mmol/m^3/d')
grid on
set(gcf,'color','w')
set(gca,'fontsize',11)
set(gca,'linewidth',1);

subplot(3,2,2),...

plot(toutmat./86400./360,smooth(INT{3}.PROD.*86400,nsmth),...
toutmat./86400./360,smooth(INT{3}.LATADV.*86400,nsmth),...
toutmat./86400./360,smooth(INT{3}.LATMIX.*86400,nsmth),...
toutmat./86400./360,smooth(INT{3}.MIX.*86400,nsmth),...
toutmat./86400./360,smooth(INT{3}.SUBD.*86400,nsmth),...
toutmat./86400./360,smooth(INT{3}.RATE.*86400,nsmth),'k-'); hold on;
if valflag == 1
plot(toutmat./86400./360,smooth(INT{3}.REAC.*86400,nsmth),'magenta--'); 
end
title('(B) Two-layer integrated nutrient budget')
legend('BIO','LADV','LMIX','VMIX','SUBD','RATE','location','northwest','orientation','horizontal')
xlabel('Year')
xlim(xlims)
ylabel('mmol/m^2/d')
grid on
set(gca,'fontsize',11)
set(gca,'linewidth',1)

subplot(3,2,3),...
    
plot(toutmat./86400./360,smooth(CONC{1}.PROD.*86400,nsmth),...
toutmat./86400./360,smooth(CONC{1}.LATADV.*86400,nsmth),...
toutmat./86400./360,smooth(CONC{1}.LATMIX.*86400,nsmth),...
toutmat./86400./360,smooth(CONC{1}.MIX.*86400,nsmth),...
toutmat./86400./360,smooth(CONC{1}.ENT.*86400,nsmth),...
toutmat./86400./360,smooth(CONC{1}.RATE.*86400,nsmth),'k-'); hold on;
if valflag == 1
    plot(toutmat./86400./360,smooth(CONC{1}.REAC.*86400,nsmth),'magenta--'); 
end
title('(C) Surface-layer nutrient concentration budget')
xlabel('Year')
xlim(xlims)
ylabel('mmol/m^3/d')
grid on
set(gca,'fontsize',11)
set(gca,'linewidth',1)

subplot(3,2,4),...
    
plot(toutmat./86400./360,smooth(INT{1}.PROD.*86400,nsmth),...
toutmat./86400./360,smooth(INT{1}.LATADV.*86400,nsmth),...
toutmat./86400./360,smooth(INT{1}.LATMIX.*86400,nsmth),...
toutmat./86400./360,smooth(INT{1}.MIX.*86400,nsmth),...
toutmat./86400./360,smooth(INT{1}.ENT.*86400,nsmth),...
toutmat./86400./360,smooth(INT{1}.RATE.*86400,nsmth),'k-'); hold on;
if valflag == 1
plot(toutmat./86400./360,smooth(INT{1}.REAC.*86400,nsmth),'magenta--'); 
end
title('(D) Depth-integrated surface-layer nutrient budget')
legend('BIO','LADV','LMIX','VMIX','ENT','RATE','location','northwest','orientation','horizontal')
xlabel('Year')
xlim(xlims)
ylabel('mmol/m^2/d')
grid on
set(gca,'fontsize',11);
set(gca,'linewidth',1);


subplot(3,2,5),...
plot(toutmat./86400./360,smooth(CONC{2}.PROD.*86400,nsmth),...
toutmat./86400./360,smooth(CONC{2}.LATADV.*86400,nsmth),...
toutmat./86400./360,smooth(CONC{2}.LATMIX.*86400,nsmth),...
toutmat./86400./360,smooth(CONC{2}.MIX.*86400,nsmth),...
toutmat./86400./360,smooth(CONC{2}.ENT.*86400,nsmth),...
toutmat./86400./360,smooth(CONC{2}.RATE.*86400,nsmth),'k-'); hold on;
if valflag == 1
plot(toutmat./86400./360,smooth(CONC{2}.REAC.*86400,nsmth),'magenta--'); 
end
title('(E) Deep-layer nutrient concentration budget')
xlabel('Year')
xlim(xlims)
ylabel('mmol/m^3/d')
grid on
set(gcf,'color','w')
set(gca,'fontsize',11)
set(gca,'linewidth',1)


subplot(3,2,6),...
plot(toutmat./86400./360,smooth(INT{2}.PROD.*86400,nsmth),...
toutmat./86400./360,smooth(INT{2}.LATADV.*86400,nsmth),...
toutmat./86400./360,smooth(INT{2}.LATMIX.*86400,nsmth),...
toutmat./86400./360,smooth(INT{2}.MIX.*86400,nsmth),...
toutmat./86400./360,smooth(INT{2}.ENT.*86400,nsmth),...
toutmat./86400./360,smooth(INT{2}.RATE.*86400,nsmth),'k-'); hold on;
if valflag == 1
    plot(toutmat./86400./360,smooth(INT{2}.REAC.*86400,nsmth),'magenta--'); 
end
title('(F) Depth-integrated deep-layer nutrient budget')
xlabel('Year')
ylabel('mmol/m^2/d')
xlim(xlims)
grid on
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 1.0;
end
set(gca,'fontsize',11)
set(gca,'linewidth',1)
end
out='SUCCESS'
end