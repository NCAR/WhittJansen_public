% run for figs 2-3
clear all;
close all;
restoredefaultpath;
addpath ./1Dmodel
addpath ./utility/

for scen = 1:3
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
    
    for ibn=10:10
        clearvars -except bgcparams Smaster ibn strscenario scen
        bgcparams.bnrate=1.2E-8.*ibn;
        
        dovis=0;
        yrs=80;
        noutperyr=41;
        load('Sequibtest.mat')
        [Soutmat,toutmat] = run_bio_NATL_box(bgcparams,Sequib,yrs,noutperyr,dovis);
        
        Smaster{ibn}.Soutmat=Soutmat;
        Smaster{ibn}.toutmat=toutmat;
        Smaster{ibn}.bgcparams=bgcparams;
    end
    save(strcat('sweep_trans_',strscenario,'.mat'),'Smaster')
end
