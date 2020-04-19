% fig S9
clear all
close all
restoredefaultpath;

addpath ./utility/cmocean/
addpath ./utility/
addpath ./1Dmodel/
addpath ./utility/export_fig/

figure('Position',[100 100 1000 350]);
subplot(1,2,2),...

load('sweep_trans_Nsub_mfj_no_phys.mat','Smaster'); ibn=2;
    Soutmat=Smaster{ibn}.Soutmat;
    toutmat=Smaster{ibn}.toutmat;
    bgcparams=Smaster{ibn}.bgcparams;
    NEW_PRODUCTION =...
        360./1000.*86400.*bgcparams.mum./(bgcparams.kz).*...
        log((bgcparams.kI+Ifn(bgcparams.latlight,mod(toutmat./86400,360)))./bgcparams.kI).* ...
        (squeeze(Soutmat(1,1,:))./(bgcparams.kN+squeeze(Soutmat(1,1,:))))';
    
    Nyear=round(3.1104e7/(toutmat(2)-toutmat(1)));
    plot(toutmat./(86400.*360),smooth(NEW_PRODUCTION,Nyear),'color',[0.8*(ibn-1),0,0],'linewidth',1,'linestyle',':');
    hold on
    xlim([2 78])
    ylim([0 1])
    grid on
    xlabel('Year','fontsize',12)
    ylabel('(mol N/m^2/yr)','fontsize',12);
    set(gcf,'color','w')

load('sweep_trans_Nsub_mfj_MLDonly.mat','Smaster'); ibn=1;
    Soutmat=Smaster{ibn}.Soutmat;
    toutmat=Smaster{ibn}.toutmat;
    bgcparams=Smaster{ibn}.bgcparams;
    NEW_PRODUCTION =...
        360./1000.*86400.*bgcparams.mum./(bgcparams.kz).*...
        log((bgcparams.kI+Ifn(bgcparams.latlight,mod(toutmat./86400,360)))./bgcparams.kI).* ...
        (squeeze(Soutmat(1,1,:))./(bgcparams.kN+squeeze(Soutmat(1,1,:))))';
    
    Nyear=round(3.1104e7/(toutmat(2)-toutmat(1)));
    plot(toutmat./(86400.*360),smooth(NEW_PRODUCTION,Nyear),'color',[0.8*(ibn-1),0,0],'linewidth',1,'linestyle','-.');

load('sweep_trans_Nsub_mfj_PSIonly.mat','Smaster'); ibn=1;
    Soutmat=Smaster{ibn}.Soutmat;
    toutmat=Smaster{ibn}.toutmat;
    bgcparams=Smaster{ibn}.bgcparams;
    NEW_PRODUCTION =...
        360./1000.*86400.*bgcparams.mum./(bgcparams.kz).*...
        log((bgcparams.kI+Ifn(bgcparams.latlight,mod(toutmat./86400,360)))./bgcparams.kI).* ...
        (squeeze(Soutmat(1,1,:))./(bgcparams.kN+squeeze(Soutmat(1,1,:))))';
    
    Nyear=round(3.1104e7/(toutmat(2)-toutmat(1)));
    plot(toutmat./(86400.*360),smooth(NEW_PRODUCTION,Nyear),'color',[0.8*(ibn-1),0,0],'linewidth',1,'linestyle','--');

load('sweep_trans_Nsub_mfj_default.mat','Smaster')
for ibn=1:2
    Soutmat=Smaster{ibn}.Soutmat;
    toutmat=Smaster{ibn}.toutmat;
    bgcparams=Smaster{ibn}.bgcparams;
    NEW_PRODUCTION =...
        360./1000.*86400.*bgcparams.mum./(bgcparams.kz).*...
        log((bgcparams.kI+Ifn(bgcparams.latlight,mod(toutmat./86400,360)))./bgcparams.kI).* ...
        (squeeze(Soutmat(1,1,:))./(bgcparams.kN+squeeze(Soutmat(1,1,:))))';
    
    Nyear=round(3.1104e7/(toutmat(2)-toutmat(1)));
    plot(toutmat./(86400.*360),smooth(NEW_PRODUCTION,Nyear),'color',[0.8*(ibn-1),0,0],'linewidth',1);
end

legend({'Declining N_{sub} (dN_{sub}/dt = -2.2 mmol/m^3/century)','Declining winter time mixed layer depth (MLD)',...
       'Declining overturning','Declining MLD and overturning',...
       'Declining MLD, overturning, and N_{sub}'},'location','Southwest');
legend boxoff
title('(b) Warming scenarios of Fig. 2a and effect of declining N_{sub}')
set(gca,'Fontsize',14);
load('CESM-LENSmean-RCP8.5-boxmodelvars_wno3_vapr112020.mat')
subplot(1,2,1)
plot(2000+yrlist(16:end),NO3_AMOC_max45deg_ensmean(16:end)./AMOC_max45deg_ensmean(16:end),'k-','linewidth',1);
hold on
plot(2000+yrlist(16:end),16.8-0.022*(yrlist(16:end)-20),':k','linewidth',1);
ylabel('(mmol/m^3)')
xlabel('Year')
title('(a) N_{sub} in CESM')
grid on
annotation(gcf,'textbox',...
    [0.29 0.63 0.19 0.05],...
    'String','-2.2 mmol/m^3 per century',...
    'LineStyle','none',...
    'FontSize',14,...
    'FontName','Arial',...
    'FitBoxToText','off');
set(gca,'fontsize',14)

    
