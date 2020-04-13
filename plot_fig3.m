% fig 3
clear all
close all
addpath ./utility/cmocean/
addpath ./utility/
addpath ./1Dmodel/
addpath ./utility/export_fig/

figure;
bnrge=1.2E-8.*(1:10);
cmap=colormap(cmocean('thermal',70));

% load data for background contours
load('sweep_steady_fig3.mat')
Smasterfig3=Smaster;
Dgs=zeros(10);
Pgs=zeros(10);
PRODgs=zeros(10);
for ibm = 1:2:20
    for ibn=1:2:20
        clear bgcnow tonow Sonow
        bgcnow=Smasterfig3{ibm,ibn}.bgcparams;
        tonow=Smasterfig3{ibm,ibn}.toutmat;
        Sonow=Smasterfig3{ibm,ibn}.Soutmat;
        Dgs(ceil(ibm/2),ceil(ibn/2))=nanmean(Dwinterfn(bgcnow,tonow));
        Pgs(ceil(ibm/2),ceil(ibn/2))=nanmean(venfn(bgcnow,tonow)...
            .*Dwinterfn(bgcnow,tonow)./bgcnow.Atllength);
        PRODgs(ceil(ibm/2),ceil(ibn/2)) = ...
        nanmean(360./1000.*86400.*bgcnow.mum./(bgcnow.kz).*...
        log((bgcnow.kI+Ifn(bgcnow.latlight,mod(tonow./86400,360)))./bgcnow.kI).* ...
        (squeeze(Sonow(1,1,:))./(bgcnow.kN+squeeze(Sonow(1,1,:))))');
    end
end


load('sweep_trans_default.mat','Smaster')
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
PsioACTL= venfn(bgcparams,toutmat).*Dwinterfn(bgcparams,toutmat)./bgcparams.Atllength;
DwodeltaCTL= Dwinterfn(bgcparams,toutmat)./bgcparams.delta;
smthprodCTL=smthprod;



clear Smaster PROD NEW_PRODUCTION toutmat Soutmat bgcparams smthprod
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
PsioAPSIonly= venfn(bgcparams,toutmat).*Dwinterfn(bgcparams,toutmat)./bgcparams.Atllength;
DwodeltaPSIonly= Dwinterfn(bgcparams,toutmat)./bgcparams.delta;
smthprodPSIonly=smthprod;



clear Smaster PROD NEW_PRODUCTION toutmat Soutmat bgcparams smthprod
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
Nyear=round(3.1104e7/(toutmat(2)-toutmat(1)));
smthprod = smooth(NEW_PRODUCTION,Nyear);
PsioAMLDonly= venfn(bgcparams,toutmat).*Dwinterfn(bgcparams,toutmat)./bgcparams.Atllength;
DwodeltaMLDonly= Dwinterfn(bgcparams,toutmat)./bgcparams.delta;
smthprodMLDonly=smthprod;


Dwodeltanow=(130:5:360)./500;
PsioAnow=(.5:.02:1.8).*1e-6;
[Pg,Dg]=meshgrid(PsioAnow,Dwodeltanow);

contourf(Pgs.*bgcparams.Atllength.*bgcparams.Atlwidth./1e6,Dgs,...
    PRODgs,.3:.01:1.3,'linestyle','none');
caxis([.3 .9]);
hold on;
caxis(caxis);
cbh=colorbar;
ylabel(cbh,'(mol N/m^2/yr)','fontsize',18)

[c,h]= contour(Pgs.*bgcparams.Atllength.*bgcparams.Atlwidth./1e6,Dgs,...
    PRODgs,.3:.1:1.3,'color',[1 1 1],'linewidth',2);
clabel(c,h,'fontsize',16);

xlabel('\Psi_{max} (10^6 m^3/s)')
ylabel('D_w (m)')

ix1yr=Nyear:205:(3291-Nyear);
plot(PsioACTL(ix1yr).*bgcparams.Atllength.*bgcparams.Atlwidth./1e6,DwodeltaCTL(ix1yr).*bgcparams.delta,'k-','linewidth',1)
plot(PsioAPSIonly(ix1yr).*bgcparams.Atllength.*bgcparams.Atlwidth./1e6,DwodeltaPSIonly(ix1yr).*bgcparams.delta,'k--','linewidth',1)
plot(PsioAMLDonly(ix1yr).*bgcparams.Atllength.*bgcparams.Atlwidth./1e6,DwodeltaMLDonly(ix1yr).*bgcparams.delta,'k-.','linewidth',1)

scatter(PsioACTL(ix1yr).*bgcparams.Atllength.*bgcparams.Atlwidth./1e6,...
    DwodeltaCTL(ix1yr).*bgcparams.delta,120,smthprodCTL(ix1yr),'filled','MarkerEdgeColor','k');
scatter(PsioAPSIonly(ix1yr).*bgcparams.Atllength.*bgcparams.Atlwidth./1e6,...
    DwodeltaPSIonly(ix1yr).*bgcparams.delta,120,smthprodPSIonly(ix1yr),'filled','MarkerEdgeColor','k');
scatter(PsioAMLDonly(ix1yr).*bgcparams.Atllength.*bgcparams.Atlwidth./1e6,...
    DwodeltaMLDonly(ix1yr).*bgcparams.delta,120,smthprodMLDonly(ix1yr),'filled','MarkerEdgeColor','k');





set(gca,'fontsize',16)
set(gca,'linewidth',1)
title('Annual new production vs D_w and \Psi_{max}','fontweight','normal')
load('CESM-LENSmean-RCP8.5-boxmodelvars_wno3_vapr112020.mat')
hold on
NO3AMOCsmth=smooth(NO3_AMOC_max45deg_ensmean,5);

AMOCsmth=smooth(AMOC_max45deg_ensmean,5);
HMXLsmth=smooth(HMXL_mean_ensmean,5);
PRODsmth=smooth(-J_NO350msummed_mean_ensmean,5);


NdCESM=sum(repmat(diff(depth_ponflux(5:11)),[1 95]).*NO3_mean_ensmean(5:10,:),1)./sum(repmat(diff(depth_ponflux(5:11)),[1 95]),1);
plot(AMOCsmth(16:5:end-2),HMXLsmth(16:5:end-2),':','color',[1 1 1],'linewidth',3);
scatter(AMOCsmth(16:5:end-2),HMXLsmth(16:5:end-2),120,-J_NO350msummed_mean_ensmean(16:5:end-2),'filled','MarkerEdgeColor','k');


annotation(gcf,'textarrow',[0.698888888888889 0.698888888888889],...
    [0.7325 0.77],'String',{'Transient box: declining overturning'},'LineWidth',1,...
    'FontSize',14,...
    'FontName','Helvetica Neue');
annotation(gcf,'textarrow',[0.714444444444444 0.738888888888889],...
    [0.6225 0.5125],'String',{'Transient box: declining overturning & mixing'},...
    'LineWidth',1,...
    'FontSize',14,...
    'FontName','Helvetica Neue');
annotation(gcf,'textarrow',[0.764444444444444 0.756666666666666],...
    [0.38 0.42],'Color',[1 1 1],...
    'String',{'CESM'},...
    'LineWidth',1,...
    'FontSize',14,...
    'FontName','Helvetica Neue');
annotation(gcf,'textarrow',[0.846666666666667 0.86],[0.2675 0.24],...
    'String',{'Transient box: declining mixing'},...
    'LineWidth',1,...
    'FontSize',14,...
    'FontName','Helvetica Neue');

set(gcf,'color','w')

%export_fig ./figures/fig_sweep_trans_physics_withCESM.pdf -painters

