% generate data for fig 3: sweep_steady_fig3.mat
clear all;
close all;
restoredefaultpath;
addpath ./1Dmodel
addpath ./utility/
bgcparams=set_bgcparams_fn();

for ibm=1:2:20
for ibn=1:2:20
clearvars -except bgcparams Smaster ibn ibm
bgcparams.bn0=0.0003*ibn;
bndmoctarget=0.0003*ibm;
bgcparams.bndmocscale=bndmoctarget./bgcparams.bn0;
bgcparams.bndmocflag = 0;

% run model to equilibrim
dovis=0;
yrs=115;
noutperyr=12;
Sequib=[10 10];
[Soutmat,toutmat] = run_bio_NATL_box(bgcparams,Sequib,yrs,noutperyr,dovis);
% run model for another year with high output frequency for plots

Soutmat=Soutmat(:,:,end-11:end);
toutmat = toutmat(end-11:end);
Ds=Dfn(bgcparams.Db,(bgcparams.Dmldfrac).*Dwinterfn(bgcparams,toutmat),mod(toutmat./86400,360));
Dd=Dwinterfn(bgcparams,toutmat)-Ds;

OUTvar(1) = max((squeeze(Ds').*squeeze(Soutmat(1,1,:))+squeeze(Dd').*squeeze(Soutmat(1,2,:)))./(squeeze(Ds'+Dd')));
OUTvar(2) = min((squeeze(Ds').*squeeze(Soutmat(1,1,:))+squeeze(Dd').*squeeze(Soutmat(1,2,:)))./(squeeze(Ds'+Dd')));
OUTvar(3) = mean((squeeze(Ds').*squeeze(Soutmat(1,1,:))+squeeze(Dd').*squeeze(Soutmat(1,2,:)))./(squeeze(Ds'+Dd')));
% concentration in the surface layer
OUTvar(4) = max(squeeze(Soutmat(1,1,:)));
OUTvar(5) = min(squeeze(Soutmat(1,1,:)));
OUTvar(6) = mean(squeeze(Soutmat(1,1,:)));
OUTvar(9) = max(squeeze(Soutmat(1,2,:)));
OUTvar(10) = min(squeeze(Soutmat(1,2,:)));
OUTvar(11) = mean(squeeze(Soutmat(1,2,:)));

% new/export production
NEW_PRODUCTION =...
    bgcparams.mum./(bgcparams.kz).*...
    log((bgcparams.kI+Ifn(bgcparams.latlight,mod(toutmat./86400,360)))./bgcparams.kI).* ...
    (squeeze(Soutmat(1,1,:))./(bgcparams.kN+squeeze(Soutmat(1,1,:))))';
OUTvar(7) = reshape(mean(NEW_PRODUCTION),[1 1 1 1 1 1 1]);

% export flux to depths below the winter mld
EXPORT_BELOW_WINTERMLD = (aflxfn(Dfn(bgcparams.Db,(bgcparams.Dmldfrac).*Dwinterfn(bgcparams,toutmat),mod(toutmat./86400,360)),bgcparams.delta) ...
    -bflxfn(Dfn(bgcparams.Db,(bgcparams.Dmldfrac).*Dwinterfn(bgcparams,toutmat),mod(toutmat./86400,360)),...
    (Dwinterfn(bgcparams,toutmat)-Dfn(bgcparams.Db,(bgcparams.Dmldfrac).*Dwinterfn(bgcparams,toutmat),mod(toutmat./86400,360))),...
    bgcparams.delta)).*...
    NEW_PRODUCTION;
OUTvar(8)=mean(EXPORT_BELOW_WINTERMLD);
dovis=0;
yrs=1;
oneyearflag=1;
noutperyr=41;
[Soutmat,toutmat] = run_bio_NATL_box(bgcparams,Soutmat(:,:,end),yrs,noutperyr,dovis);
Smaster{ibm,ibn}.Soutmat=Soutmat;
Smaster{ibm,ibn}.toutmat=toutmat;
Smaster{ibm,ibn}.bgcparams=bgcparams;
Smaster{ibm,ibn}.OUTvar=OUTvar;
end
end
save('sweep_steady_fig3.mat','Smaster')
