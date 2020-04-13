% fig S5
clear all
close all
restoredefaultpath;
addpath ./utility/cmocean/
addpath ./utility/
addpath ./1Dmodel/
addpath ./utility/export_fig/
figure('Position',[100 100 1500 525]);
cmap=colormap(cmocean('thermal',8));

idx=[1,7:8] % (top)
%idx=[1,4:6] % (middle)
%idx=[1,2:3]  % (bottom)

for it = [idx]
load('sweep_steady_c5toc7_default.mat','Smaster')
Smaster=Smaster(it,:);
bnrge=0.0002+0.0003.*(1:2:40);
cmap=colormap(cmocean('thermal',8));
    for ibn = 1:2:40
        clear  NEW_PRODUCTION OUTvar bgcparams toutmat Soutmat
        Soutmat=Smaster{ibn}.Soutmat;
        toutmat=Smaster{ibn}.toutmat;
        bgcparams=Smaster{ibn}.bgcparams;
        OUTvar=Smaster{ibn}.OUTvar;
        PROD(ceil(ibn./2))=OUTvar(7);
        SINK(ceil(ibn./2))=OUTvar(8);
        NEW_PRODUCTION =...
    bgcparams.mum./(bgcparams.kz).*...
    log((bgcparams.kI+Ifn(bgcparams.latlight,mod(toutmat./86400,360)))./bgcparams.kI).* ...
    (squeeze(Soutmat(1,1,:))./(bgcparams.kN+squeeze(Soutmat(1,1,:))))';
    end

subplot(1,2,[1 2]),...
    hold on
    plot(abs(0.0005-bnrge)./.0096,360./1000.*86400.*PROD(:),'-','linewidth',1,'color',cmap(it,:));
hold on;
title('Sensitivity of steady-state new production to surface density','fontsize',12,'fontweight','normal');
xlabel('Surface density anomaly (kg/m^3)','fontsize',12)
ylabel('Annual mean nutrient consumption rate (mol N/m^2/yr)','fontsize',12)
grid on

%%

%RATIO1=SINK(:,:)./PROD(:,:);
clear Smaster PROD SINK
pause(0.1)
load('sweep_steady_c5toc7_PSIonly.mat','Smaster')
pause(0.1)
Smaster=Smaster(it,:);
for ibn = 1:2:40
    clear  NEW_PRODUCTION OUTvar bgcparams toutmat Soutmat
    Soutmat=Smaster{ibn}.Soutmat;
    toutmat=Smaster{ibn}.toutmat;
    bgcparams=Smaster{ibn}.bgcparams;
    OUTvar=Smaster{ibn}.OUTvar;
    PROD(ceil(ibn./2))=OUTvar(7);
    SINK(ceil(ibn./2))=OUTvar(8);
    NEW_PRODUCTION =...
    bgcparams.mum./(bgcparams.kz).*...
    log((bgcparams.kI+Ifn(bgcparams.latlight,mod(toutmat./86400,360)))./bgcparams.kI).* ...
    (squeeze(Soutmat(1,1,:))./(bgcparams.kN+squeeze(Soutmat(1,1,:))))';
end

subplot(1,2,[1 2]),...
hold on
plot(abs(0.0005-bnrge)./.0096,360./1000.*86400.*PROD(:),'--','linewidth',1,'color',cmap(it,:));
hold on;







clear Smaster PROD SINK
pause(0.1)
load('sweep_steady_c5toc7_MLDonly.mat','Smaster')
Smaster=Smaster(it,:);
for ibn = 1:2:40
    clear  NEW_PRODUCTION OUTvar bgcparams toutmat Soutmat
    Soutmat=Smaster{ibn}.Soutmat;
    toutmat=Smaster{ibn}.toutmat;
    bgcparams=Smaster{ibn}.bgcparams;
    OUTvar=Smaster{ibn}.OUTvar;
    PROD(ceil(ibn./2))=OUTvar(7);
    SINK(ceil(ibn./2))=OUTvar(8);
    NEW_PRODUCTION =...
    bgcparams.mum./(bgcparams.kz).*...
    log((bgcparams.kI+Ifn(bgcparams.latlight,mod(toutmat./86400,360)))./bgcparams.kI).* ...
    (squeeze(Soutmat(1,1,:))./(bgcparams.kN+squeeze(Soutmat(1,1,:))))';
    PROD(ceil(ibn./2))=nanmean(NEW_PRODUCTION);

end



subplot(1,2,[1 2]),...
plot(abs(0.0005-bnrge)./0.0096,360./1000.*86400.*PROD(:),'-.','linewidth',1,'color',cmap(it,:));
ylim([0 1])
set(gca,'fontsize',16)
set(gcf,'color','w')
set(gca,'linewidth',1)
hold on
end

xlim([0 1])

set(gca,'linewidth',1)
