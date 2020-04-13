%% Script to run example scenarios and generate some supplemental figs S1-S3
clear all;
close all;
restoredefaultpath;
addpath ./1Dmodel
addpath ./utility/


%% Start with climatechangeflag =0; fig s1 a-b; fig s2
%% then run with
%% climatechangeflag = 1 (MLD and PSI changing together); fig S3
%% climatechangeflag = 2 (only MLD changing, fixed PSI)
%% climatechangeflag = 3 (only PSI changing, fixed MLD)
climatechangeflag = 0;
%%

bgcparams=set_bgcparams_fn();
if climatechangeflag == 1
bgcparams.bnrate=1.2e-7;
yrs=80;
load('Sequibtest.mat');
elseif climatechangeflag == 2
bgcparams.bnrate=1.2e-7;
yrs=80;
load('Sequibtest.mat');
bgcparams.bnpsiflag=0;  
elseif climatechangeflag == 3
bgcparams.bnrate=1.2e-7;
yrs=80;
load('Sequibtest.mat');
bgcparams.bndmocflag=0;    
elseif climatechangeflag==0
yrs=200;
Sequib=[10 10];
else
    disp('CASE NOT SUPPORTED')
    pause
end

% run model to equilibrim
noutperyr=360;
dovis=1;
 [Soutmat,toutmat] = run_bio_NATL_box(bgcparams,Sequib,yrs,noutperyr,dovis);
 set(gcf,'color','w')
 Sequib=Soutmat(:,:,end);
 
 if climatechangeflag == 1
save 'Sequib2100.mat' Sequib
 elseif climatechangeflag == 2
save 'Sequib2100MLDonly.mat' Sequib
  elseif climatechangeflag == 3
save 'Sequib2100PSIonly.mat' Sequib
 elseif climatechangeflag == 0
 save 'Sequibtest.mat' Sequib
 end
 

% plot the budgets from years 1-2, years 78-79, 
%and years 2-78 (with 1-year smoother applied)
Nyear=round((86400.*360)/(toutmat(2)-toutmat(1)));

[REAC,CONC,INT,Dnow,Dwinter,wen,a,b,prod] = construct_reac_NATL_box_budget(Soutmat,1,toutmat,bgcparams);
valflag =0;
nsmth=round((86400.*30)/(toutmat(2)-toutmat(1))); % 2 months
xlims{1}=[1 2]
xlims{2}=[78 79]
xlims{3}=[2 78]
out=plot_3_budgets(toutmat,CONC,INT,Nyear,valflag,nsmth,xlims)



%  run model for another 1 year period to save and plot "equilibrium" year
if climatechangeflag==0
clear toutmat
 yrs=1;
 oneyearflag=1;
 noutperyr=360;
 dovis=0;
 [Soutmat,toutmat] = run_bio_NATL_box(bgcparams,Soutmat(:,:,end),yrs,noutperyr,dovis);
save('one_year_default_params.mat','Soutmat','toutmat','bgcparams')
Nyear=round((86400.*30)/(toutmat(2)-toutmat(1)));
valflag =0;
nsmth=round((86400.*30)/(toutmat(2)-toutmat(1))); % 2 months
xlims{1}=[1/24 23/24]
xlims{2}=[1/24 23/24]
xlims{3}=[1/24 23/24]
[REAC,CONC,INT,Dnow,Dwinter,wen,a,b,prod] = construct_reac_NATL_box_budget(Soutmat,1,toutmat,bgcparams);
out=plot_3_budgets(toutmat,CONC,INT,Nyear,valflag,nsmth,xlims)

    
Naout=bgcparams.Na;
depth=linspace(-500,0,251);
depthg=repmat(depth',[1 noutperyr]);
Ds=-Dfn(bgcparams.Db,(bgcparams.Dmldfrac).*Dwinterfn(bgcparams,toutmat),mod(toutmat./86400,360));
Dsg=repmat(Ds,[size(depthg,1) 1]);
maskDs = depthg>Dsg;
maskDd = depthg<=Dsg;
maskDa = depthg<=-Dwinterfn(bgcparams,toutmat);
Ns=squeeze(Soutmat(1,1,:));
Nd=squeeze(Soutmat(1,2,:));
Nsg=repmat(Ns',[size(depthg,1) 1]);
Ndg=repmat(Nd',[size(depthg,1) 1]);
Ng=Naout.*ones(size(depthg));
Ng(maskDs)=Nsg(maskDs);
Ng(maskDd)=Ndg(maskDd);
Ng(maskDa)=Naout;
figure;
contourf((((360/noutperyr)/2):(360/noutperyr):(360-360/noutperyr/2))./30,depthg(:,1),Ng,linspace(0,20,21),'linestyle','none');
hold on;
xlim([0 12])
caxis([0 18])
ylim([-500 0])
caxis(caxis)
hold on;
plot((((360/noutperyr)/2):(360/noutperyr):(360-360/noutperyr/2))./30,Ds,'k-','linewidth',1);
ylabel('Depth (m)');
title('Area averaged nitrate profile')
addpath ./utility/cmocean/
colormap(cmocean('dense',32))


cbh=colorbar('southoutside');
xlabel(cbh,'nitrate (mmol/m^3)')
hold on

set(gcf,'color','w')
xlabel('Months')
ylabel('Depth (m)');



annual_cycle_2pan_surface_layer_budget_plot;

end

