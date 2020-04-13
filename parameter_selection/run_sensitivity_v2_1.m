clear all;
addpath ../1Dmodel
% to add a smooth function (missing toolbox)
addpath ../utility/
% bgc sensitivity

OUT=zeros(1728,12,12,11);
deltarge = linspace(50,750,12)
mumrge = linspace(.1,1.5,12)./86400
Dwrge = linspace(.05,.5,12)
kNrge = linspace(.5,8,12)
Nrge = linspace(12,22,12)
parfor it = 1:1728
    it1 = ceil(it./144)
    tempidx=mod(it,144);
    tempidx(tempidx==0)=144;
    it2 = ceil(tempidx./12)
    it3=mod(it,12);
    it3(it3==0)=12
    OUTvar=zeros(12,12,11);
   
            for it4 = 1:12
                for it5 = 1:12
                        bgcparams=set_bgcparams_fn();
                        bgcparams.delta=deltarge(it1);
                        bgcparams.mum=mumrge(it2);
                        bgcparams.Dwinterfrac=Dwrge(it3)
                        bgcparams.kN =kNrge(it4);
                        bgcparams.Nsub=Nrge(it5);
                        bgcparams.Na= 18;
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
                        % vertically-averaged concentration above the winter MLD
                        
                        OUTvar(it4,it5,1) = max((squeeze(Ds').*squeeze(Soutmat(1,1,:))+squeeze(Dd').*squeeze(Soutmat(1,2,:)))./(squeeze(Ds'+Dd')));
                        OUTvar(it4,it5,2) = min((squeeze(Ds').*squeeze(Soutmat(1,1,:))+squeeze(Dd').*squeeze(Soutmat(1,2,:)))./(squeeze(Ds'+Dd')));
                        OUTvar(it4,it5,3) = mean((squeeze(Ds').*squeeze(Soutmat(1,1,:))+squeeze(Dd').*squeeze(Soutmat(1,2,:)))./(squeeze(Ds'+Dd')));
                        % concentration in the surface layer
                        OUTvar(it4,it5,4) = max(squeeze(Soutmat(1,1,:)));
                        OUTvar(it4,it5,5) = min(squeeze(Soutmat(1,1,:)));
                        OUTvar(it4,it5,6) = mean(squeeze(Soutmat(1,1,:)));
                        % concentration in the surface layer
                        OUTvar(it4,it5,9) = max(squeeze(Soutmat(1,2,:)));
                        OUTvar(it4,it5,10) = min(squeeze(Soutmat(1,2,:)));
                        OUTvar(it4,it5,11) = mean(squeeze(Soutmat(1,2,:)));
                        
                        % new/export production
                        NEW_PRODUCTION =...
                            bgcparams.mum./(bgcparams.kz).*...
                            log((bgcparams.kI+Ifn(bgcparams.latlight,mod(toutmat./86400,360)))./bgcparams.kI).* ...
                            (squeeze(Soutmat(1,1,:))./(bgcparams.kN+squeeze(Soutmat(1,1,:))))';
                        OUTvar(it4,it5,7) = reshape(mean(NEW_PRODUCTION),[1 1 1 1 1 1 1]);
                        
                        % export flux to depths below the winter mld
                        EXPORT_BELOW_WINTERMLD = (aflxfn(Dfn(bgcparams.Db,(bgcparams.Dmldfrac).*Dwinterfn(bgcparams,toutmat),mod(toutmat./86400,360)),bgcparams.delta) ...
                            -bflxfn(Dfn(bgcparams.Db,(bgcparams.Dmldfrac).*Dwinterfn(bgcparams,toutmat),mod(toutmat./86400,360)),...
                            (Dwinterfn(bgcparams,toutmat)-Dfn(bgcparams.Db,(bgcparams.Dmldfrac).*Dwinterfn(bgcparams,toutmat),mod(toutmat./86400,360))),...
                            bgcparams.delta)).*...
                            NEW_PRODUCTION;
                        OUTvar(it4,it5,8)=mean(EXPORT_BELOW_WINTERMLD); 
                end
            end
    OUT(it,:,:,:)=OUTvar;
end
OUT2=zeros(12,12,12,12,12,11);
for it = 1:1728
    it1 = ceil(it./144)
    tempidx=mod(it,144);
    tempidx(tempidx==0)=144;
    it2 = ceil(tempidx./12)
    it3=mod(it,12);
    it3(it3==0)=12
    OUT2(it1,it2,it3,:,:,:)=OUT(it,:,:,:);
end
OUT=OUT2;
clear OUTvar OUT2
save(strcat('sensitivity_v2.mat'))



