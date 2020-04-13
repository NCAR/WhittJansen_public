% run for fig S8
clear all;
close all;
restoredefaultpath;
addpath ./1Dmodel
addpath ./utility/

for scen = 2:3
    clearvars -except scen
    if scen == 1
        strscenario='default'
    elseif scen == 2
        strscenario='PSIonly'
    elseif scen == 3
        strscenario='MLDonly'
    end
    
    bgcparams=set_bgcparams_fn();
    if  scen == 2
        bgcparams.bndmocflag = 0;
    elseif scen == 3
        bgcparams.bnpsiflag =0;
    end
    
    for ibn=1:10
        clearvars -except bgcparams Smaster ibn strscenario scen
        bgcparams.bnrate=1.2E-7;
        bgcparams.Nsubrate=-(ibn-1)*(.5./80./360); % decline at a rate between 0 and 5 mmol/m3 per 80 years
        
        dovis=0;
        yrs=80;
        noutperyr=41;
        load('Sequibtest.mat')
        [Soutmat,toutmat] = run_bio_NATL_box(bgcparams,Sequib,yrs,noutperyr,dovis);
        
        Smaster{ibn}.Soutmat=Soutmat;
        Smaster{ibn}.toutmat=toutmat;
        Smaster{ibn}.bgcparams=bgcparams;
    end
    save(strcat('sweep_trans_Nsub_',strscenario,'.mat'),'Smaster')
end
