% fig S9
clear all
close all
restoredefaultpath;

addpath ./utility/cmocean/
addpath ./utility/
addpath ./1Dmodel/
addpath ./utility/export_fig/

figure('Position',[100 100 1000 350]);
subplot(1,2,1),...
cmap=colormap(cmocean('thermal',5));
load('sweep_trans_Nsub_default.mat','Smaster')
for ibn=1:2:10
    Soutmat=Smaster{ibn}.Soutmat;
toutmat=Smaster{ibn}.toutmat;
bgcparams=Smaster{ibn}.bgcparams;
NEW_PRODUCTION =...
    360./1000.*86400.*bgcparams.mum./(bgcparams.kz).*...
    log((bgcparams.kI+Ifn(bgcparams.latlight,mod(toutmat./86400,360)))./bgcparams.kI).* ...
    (squeeze(Soutmat(1,1,:))./(bgcparams.kN+squeeze(Soutmat(1,1,:))))';

Nyear=round(3.1104e7/(toutmat(2)-toutmat(1)));
plot(toutmat./(86400.*360),smooth(NEW_PRODUCTION,Nyear),'color',cmap(round(ibn/2),:),'linewidth',1);
hold on
xlim([2 78])
ylim([0 1])
grid on
xlabel('Year','fontsize',12)
ylabel('(mol N/m^2/yr)','fontsize',12);
set(gcf,'color','w')
end

load('sweep_trans_Nsub_MLDonly.mat','Smaster')
for ibn=1:2:10
    Soutmat=Smaster{ibn}.Soutmat;
    toutmat=Smaster{ibn}.toutmat;
    bgcparams=Smaster{ibn}.bgcparams;
    NEW_PRODUCTION =...
        360./1000.*86400.*bgcparams.mum./(bgcparams.kz).*...
        log((bgcparams.kI+Ifn(bgcparams.latlight,mod(toutmat./86400,360)))./bgcparams.kI).* ...
        (squeeze(Soutmat(1,1,:))./(bgcparams.kN+squeeze(Soutmat(1,1,:))))';
    
    Nyear=round(3.1104e7/(toutmat(2)-toutmat(1)));
    plot(toutmat./(86400.*360),smooth(NEW_PRODUCTION,Nyear),'color',cmap(round(ibn/2),:),'linewidth',1,'linestyle','-.');
    hold on
    xlim([2 78])
    ylim([0 1])
    grid on
    xlabel('Year','fontsize',12)
    ylabel('(mol N/m^2/yr)','fontsize',12);
    set(gcf,'color','w')
end



load('sweep_trans_Nsub_PSIonly.mat','Smaster')
for ibn=1:2:10
    Soutmat=Smaster{ibn}.Soutmat;
    toutmat=Smaster{ibn}.toutmat;
    bgcparams=Smaster{ibn}.bgcparams;
    NEW_PRODUCTION =...
        360./1000.*86400.*bgcparams.mum./(bgcparams.kz).*...
        log((bgcparams.kI+Ifn(bgcparams.latlight,mod(toutmat./86400,360)))./bgcparams.kI).* ...
        (squeeze(Soutmat(1,1,:))./(bgcparams.kN+squeeze(Soutmat(1,1,:))))';
    
    Nyear=round(3.1104e7/(toutmat(2)-toutmat(1)));
    plot(toutmat./(86400.*360),smooth(NEW_PRODUCTION,Nyear),'color',cmap(round(ibn/2),:),'linewidth',1,'linestyle','--');
    hold on
    xlim([2 78])
    ylim([0 1])
    grid on
    xlabel('Year','fontsize',12)
    ylabel('(mol N/m^2/yr)','fontsize',12);
    set(gcf,'color','w')
end
cbh=colorbar;
set(cbh,'YTick',0:.2:.8)
set(cbh,'YTickLabel',{'0','-1.25','-2.5','-3.75','-5'})
ylabel(cbh,'d N_{sub}/dt (mmol/m^3 per century)','fontsize',14);
title('(a) Warming scenarios in Fig. 2a with declining N_{sub}')
set(gca,'Fontsize',14);
load('CESM-LENSmean-RCP8.5-boxmodelvars_wno3_vapr112020.mat')
subplot(1,2,2),...
    plot(2000+yrlist(16:end),NO3_AMOC_max45deg_ensmean(16:end)./AMOC_max45deg_ensmean(16:end),'k-','linewidth',1);
ylabel('(mmol/m^3)')
xlabel('year')
title('(b) N_{sub} in CESM')
grid on
annotation(gcf,'textbox',...
    [0.73 0.867 0.19 0.05],...
    'String','-2.2 mmol/m3 per century',...
    'LineStyle','none',...
    'FontSize',14,...
    'FontName','Arial',...
    'FitBoxToText','off');
set(gca,'fontsize',14)

set(gca,'Fontsize',13);
set(gcf,'color','w')
    
