% script to generate mat file data for Fig S5 
clear all;
close all;
addpath ./1Dmodel
addpath ./utility/
for scen=1:3
    clearvars -except scen
    if scen == 1
        strscenario='default'
    elseif scen == 2
        strscenario='PSIonly'
    elseif scen == 3
        strscenario='MLDonly'
    end
    for it = 1:8
        clear bgcparams
        bgcparams=set_bgcparams_fn();
        if  scen == 2
            bgcparams.bndmocflag = 0;
        elseif scen == 3
            bgcparams.bnpsiflag =0;
        end
        if it == 1
        elseif it == 2
            % delta/Dw fixed  = 1.56
            bgcparams.delta = 700;
            bgcparams.Dwinterfrac=0.24;
        elseif it == 3
            bgcparams.delta = 300;
            bgcparams.Dwinterfrac=0.10;
        elseif it == 4
            bgcparams.Dwinterfrac = 0.11;
            bgcparams.Nsub=22;
        elseif it == 5
            bgcparams.Dwinterfrac=0.27;
            bgcparams.Nsub=17;
        elseif it == 6
            bgcparams.Dwinterfrac=0.35;
            bgcparams.Nsub=15.5;
        elseif it == 7
            bgcparams.mum=.65./86400;
            bgcparams.kN=5.3;
        elseif it == 8
            bgcparams.mum=.35./86400;
            bgcparams.kN=1.2;
        end
        
        for ibn=1:2:40
            clearvars -except bgcparams Smaster ibn it scen strscenario
            bgcparams.bn0=0.0002+0.0003*ibn;
            if scen == 2
                bndmoctarget=0.0005;
                bgcparams.bndmocscale=bndmoctarget./bgcparams.bn0;
            end
            
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
            Smaster{it,ibn}.Soutmat=Soutmat;
            Smaster{it,ibn}.toutmat=toutmat;
            Smaster{it,ibn}.bgcparams=bgcparams;
            Smaster{it,ibn}.OUTvar=OUTvar;
        end
    end
    save(strcat('sweep_steady_c5toc7_',strscenario,'.mat'),'Smaster')
end
