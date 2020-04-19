% fig S4
clear all
close all
addpath ./utility/cmocean/
addpath ./utility/
addpath ./1Dmodel/
addpath ./utility/export_fig/
figure('Position',[100 100 1000 1000]);
load('sweep_trans_default.mat','Smaster')
bnrge=1.2E-8.*(1:10);
cmap=colormap(cmocean('thermal',10));

ibn=10;
Soutmat=Smaster{ibn}.Soutmat;
toutmat=Smaster{ibn}.toutmat;
bgcparams=Smaster{ibn}.bgcparams;
NEW_PRODUCTION =...
    360./1000.*86400.*bgcparams.mum./(bgcparams.kz).*...
    log((bgcparams.kI+Ifn(bgcparams.latlight,mod(toutmat./86400,360)))./bgcparams.kI).* ...
    (squeeze(Soutmat(1,1,:))./(bgcparams.kN+squeeze(Soutmat(1,1,:))))';

Nyear=round(3.1104e7/(toutmat(2)-toutmat(1)));
smthprod = smooth(NEW_PRODUCTION,Nyear);
smthNd=smooth(squeeze(Soutmat(1,2,:)),Nyear);
Ffactor=venfn(bgcparams,toutmat).*Dwinterfn(bgcparams,toutmat)./bgcparams.Atllength./exp(-Dwinterfn(bgcparams,toutmat)./bgcparams.delta).*86400;
PsioACTL= venfn(bgcparams,toutmat).*Dwinterfn(bgcparams,toutmat)./bgcparams.Atllength;
DwodeltaCTL= Dwinterfn(bgcparams,toutmat)./bgcparams.delta;
smthprodCTL=smthprod;

        subplot(2,2,4),...
            hold on;
        plot(Ffactor(Nyear:end-Nyear),...
            1000./360.*smthprod(Nyear:end-Nyear)./(19-smthNd(Nyear:end-Nyear)'),'k-','linewidth',1);
        hold on
        grid on
title('(D) Nutrient consumption vs physical forcing velocity scales','fontsize',12,'fontweight','normal');
xlabel('\Psi_{max}/A e^{D_w/\delta} (m/d)','fontsize',12)
ylabel('PROD/\Delta N (m/d)','fontsize',12)
grid on
set(gca,'FontName','Arial')


%%
subplot(2,2,2),...
    hold on
        plot((19-smthNd(Nyear:end-Nyear)),smthprod(Nyear:end-Nyear),'k-','linewidth',1);

title('(B) Nutrient consumption vs deficit','fontsize',12,'fontweight','normal');
xlabel('\Delta N (mmol/m^3)','fontsize',12)
ylabel('(mol N/m^2/yr)','fontsize',12)
grid on
set(gca,'FontName','Arial')

subplot(2,2,1),...
    hold on
        plot((smthNd(Nyear:end-Nyear)),smthprod(Nyear:end-Nyear),'k-','linewidth',1);

title('(A) Nutrient consumption vs concentration','fontsize',12,'fontweight','normal');
xlabel('N_d (mmol/m^3)','fontsize',16)
ylabel('(mol N/m^2/yr)','fontsize',16)
grid on
set(gca,'FontName','Arial')



%%



clear Smaster PROD NEW_PRODUCTION toutmat Soutmat bgcparams Ffactor smthNd smthprod
load('sweep_trans_PSIonly.mat','Smaster')
ibn=10;
Soutmat=Smaster{ibn}.Soutmat;
toutmat=Smaster{ibn}.toutmat;
bgcparams=Smaster{ibn}.bgcparams;
bgcparams.bndmocflag=0;
NEW_PRODUCTION =...
    360./1000.*86400.*bgcparams.mum./(bgcparams.kz).*...
    log((bgcparams.kI+Ifn(bgcparams.latlight,mod(toutmat./86400,360)))./bgcparams.kI).* ...
    (squeeze(Soutmat(1,1,:))./(bgcparams.kN+squeeze(Soutmat(1,1,:))))';
Nyear=round(3.1104e7/(toutmat(2)-toutmat(1)));
smthprod = smooth(NEW_PRODUCTION,Nyear);
smthNd=smooth(squeeze(Soutmat(1,2,:)),Nyear);
Ffactor=venfn(bgcparams,toutmat).*Dwinterfn(bgcparams,toutmat)./bgcparams.Atllength./exp(-Dwinterfn(bgcparams,toutmat)./bgcparams.delta).*86400;


PsioAPSIonly= venfn(bgcparams,toutmat).*Dwinterfn(bgcparams,toutmat)./bgcparams.Atllength;
DwodeltaPSIonly= Dwinterfn(bgcparams,toutmat)./bgcparams.delta;
smthprodPSIonly=smthprod;


subplot(2,2,2),...
    hold on
        plot((19-smthNd(Nyear:end-Nyear)),smthprod(Nyear:end-Nyear),'k--','linewidth',1);


subplot(2,2,1),...
    hold on
        plot(smthNd(Nyear:end-Nyear),smthprod(Nyear:end-Nyear),'k--','linewidth',1);

subplot(2,2,4),...
        plot(Ffactor(Nyear:end-Nyear),...
            1000./360.*smthprod(Nyear:end-Nyear)./(19-smthNd(Nyear:end-Nyear)'),'k--','linewidth',1);




clear Smaster PROD NEW_PRODUCTION toutmat Soutmat bgcparams Ffactor smthNd smthprod
load('sweep_trans_MLDonly.mat','Smaster')
ibn=10;
Soutmat=Smaster{ibn}.Soutmat;
toutmat=Smaster{ibn}.toutmat;
bgcparams=Smaster{ibn}.bgcparams;
bgcparams.bnpsiflag=0;
NEW_PRODUCTION =...
    360./1000.*86400.*bgcparams.mum./(bgcparams.kz).*...
    log((bgcparams.kI+Ifn(bgcparams.latlight,mod(toutmat./86400,360)))./bgcparams.kI).* ...
    (squeeze(Soutmat(1,1,:))./(bgcparams.kN+squeeze(Soutmat(1,1,:))))';
PROD(ibn) = nanmean(NEW_PRODUCTION(end-40:end));
Nyear=round(3.1104e7/(toutmat(2)-toutmat(1)));
smthprod = smooth(NEW_PRODUCTION,Nyear);
smthNd=smooth(squeeze(Soutmat(1,2,:)),Nyear);
Ffactor=venfn(bgcparams,toutmat).*Dwinterfn(bgcparams,toutmat)./bgcparams.Atllength./exp(-Dwinterfn(bgcparams,toutmat)./bgcparams.delta).*86400;

PsioAMLDonly= venfn(bgcparams,toutmat).*Dwinterfn(bgcparams,toutmat)./bgcparams.Atllength;
DwodeltaMLDonly= Dwinterfn(bgcparams,toutmat)./bgcparams.delta;
smthprodMLDonly=smthprod;




subplot(2,2,2),...
    hold on
        plot((19-smthNd(Nyear:end-Nyear)),smthprod(Nyear:end-Nyear),'k-.','linewidth',1);

set(gca,'fontsize',16)
set(gcf,'color','w')
set(gca,'linewidth',1)
ylim([0 1])
hold on

subplot(2,2,1),...
    hold on
        plot(smthNd(Nyear:end-Nyear),smthprod(Nyear:end-Nyear),'k-.','linewidth',1);

set(gca,'fontsize',16)
set(gcf,'color','w')
set(gca,'linewidth',1)
ylim([0 1])
hold on

subplot(2,2,4),...
        plot(Ffactor(Nyear:end-Nyear),...
            1000./360.*smthprod(Nyear:end-Nyear)./(19-smthNd(Nyear:end-Nyear)'),'k-.','linewidth',1);


set(gca,'FontName','Arial')
set(gca,'linewidth',1)




load('CESM-LENSmean-RCP8.5-boxmodelvars_wno3_vapr112020.mat')
hold on
NO3AMOCsmth=smooth(NO3_AMOC_max45deg_ensmean,5);

AMOCsmth=smooth(AMOC_max45deg_ensmean,5);
HMXLsmth=smooth(HMXL_mean_ensmean,5);
JNO3502020=mean(J_NO350msummed_mean_ensmean(15:17))
JNO3502099=mean(J_NO350msummed_mean_ensmean(93:95))
nanmean(NO3_AMOC_max45deg_ensmean(15:17))
nanmean(NO3_AMOC_max45deg_ensmean(93:95))

NdCESM=sum(repmat(diff(depth_ponflux(5:11)),[1 95]).*NO3_augmean_ensmean(5:10,:),1)./sum(repmat(diff(depth_ponflux(5:11)),[1 95]),1);

subplot(2,2,2),...
    hold on;
    plot((NO3AMOCsmth(16:end-2)./AMOCsmth(16:end-2)-NdCESM(16:end-2)),-J_NO350msummed_mean_ensmean(16:end-2),':','color',[.5 .5 .5],'linewidth',2);
    plot(exp(1/2).*(NO3AMOCsmth(16:end-2)./AMOCsmth(16:end-2)-NdCESM(16:end-2)),-J_NO350msummed_mean_ensmean(16:end-2),'b:','linewidth',2);


subplot(2,2,1),...
    hold on;
    plot(NdCESM(16:end-2),-J_NO350msummed_mean_ensmean(16:end-2),':','color',[.5 .5 .5],'linewidth',2);

   subplot(2,2,4),...
   hold on;
   plot(AMOCsmth(16:end-2)./11e6.*86400.*exp(HMXLsmth(16:end-2)./250),-1000./365.*J_NO350msummed_mean_ensmean(16:end-2)./...
       (NO3AMOCsmth(16:end-2)./AMOCsmth(16:end-2)-NdCESM(16:end-2))',':','color',[.5 .5 .5],'linewidth',2);

   plot(exp(-1./2).*AMOCsmth(16:end-2)./11e6.*86400.*exp(HMXLsmth(16:end-2)./250),-exp(-1./2).*1000./365.*J_NO350msummed_mean_ensmean(16:end-2)./...
       (NO3AMOCsmth(16:end-2)./AMOCsmth(16:end-2)-NdCESM(16:end-2))','b:','linewidth',2);

   
   set(gca,'fontsize',16)
set(gcf,'color','w')
set(gca,'linewidth',1)
xlim([0 .5])
ylim([0 .5])

   

subplot(2,2,3),...
plot(PON_FLUXz_mean_ensmean(:,16),-depth_ponflux,':','color',[.5 .5 .5],'linewidth',2);
hold on
plot(exp((-depth_ponflux+50)./500).*.5,-depth_ponflux,'k-','linewidth',2);
plot(exp((-depth_ponflux+50)./250).*.5,-depth_ponflux,'k--','linewidth',2);


ylim([-500 -50])
title('(C) Sinking particle flux profiles','fontsize',14,'fontweight','normal')

xlabel('(mol N/m^2/yr)','fontsize',14)
ylabel('Depth (m)')
legend('CESM-LE/2020','Exponential; \delta=500 m','Exponential; \delta=250 m','location','northwest')
grid on
set(gca,'fontsize',16)
set(gcf,'color','w')
set(gca,'linewidth',1)
