 % only plot every 8th day:
 noutperyr=noutperyr/8;
Soutmat=Soutmat(:,:,1:8:end);
toutmat=toutmat(1:8:end);

fig=figure('Position', [10 10 800 1200]);
subplot(4,2,1),...
plot(toutmat./86400,squeeze(Soutmat(:,1,:)),...
    toutmat./86400,squeeze(Soutmat(:,2,:)));
title('N concentrations (mmol N/m^3)')
xlabel('days since start')
legend('N_s','N_d')
ylim([0 15])

PROD_SURF_SUM=(...
    bgcparams.mum./(bgcparams.kz).*...
                      log((bgcparams.kI+Ifn(bgcparams.latlight,mod(toutmat./86400,360)))./bgcparams.kI).* ...
                     (squeeze(Soutmat(:,1,:))./(bgcparams.kN+squeeze(Soutmat(:,1,:))))');
                 
EXPORT_LOSS=((aflxfn(Dfn(bgcparams.Db,(bgcparams.Dmldfrac).*Dwinterfn(bgcparams,toutmat),mod(toutmat./86400,360)),bgcparams.delta))...
    .*bgcparams.mum./(bgcparams.kz).*...
                      log((bgcparams.kI+Ifn(bgcparams.latlight,mod(toutmat./86400,360)))./bgcparams.kI).* ...
                     (squeeze(Soutmat(:,1,:))./(bgcparams.kN+squeeze(Soutmat(:,1,:))))') - ...
                     bflxfn(Dfn(bgcparams.Db,(bgcparams.Dmldfrac).*Dwinterfn(bgcparams,toutmat),mod(toutmat./86400,360)),(Dwinterfn(bgcparams,toutmat)-Dfn(bgcparams.Db,(bgcparams.Dmldfrac).*Dwinterfn(bgcparams,toutmat),mod(toutmat./86400,360))),bgcparams.delta).*...
                     bgcparams.mum./(bgcparams.kz).*...
                      log((bgcparams.kI+Ifn(bgcparams.latlight,mod(toutmat./86400,360)))./bgcparams.kI).* ...
                     (squeeze(Soutmat(:,1,:))./(bgcparams.kN+squeeze(Soutmat(:,1,:))))';
PROD_SURF_LOSS=aflxfn(Dfn(bgcparams.Db,(bgcparams.Dmldfrac).*Dwinterfn(bgcparams,toutmat),mod(toutmat./86400,360)),bgcparams.delta).*...
    PROD_SURF_SUM;

subplot(4,2,2),...
plot(toutmat./86400,PROD_SURF_SUM.*bgcparams.Atlwidth.*bgcparams.Atllength./1e6,toutmat./86400,PROD_SURF_LOSS.*bgcparams.Atlwidth.*bgcparams.Atllength./1e6,...
toutmat./86400,EXPORT_LOSS.*bgcparams.Atlwidth.*bgcparams.Atllength./1e6); hold on;
hold on;
if oneyearflag~=1
Nyear=round((360.*86400)/(toutmat(2)-toutmat(1)));
plot(toutmat./86400,smooth(PROD_SURF_SUM.*bgcparams.Atlwidth.*bgcparams.Atllength./1e6,Nyear),'r-',...
    toutmat./86400,smooth(PROD_SURF_LOSS.*bgcparams.Atlwidth.*bgcparams.Atllength./1e6,Nyear),'r-',...
    toutmat./86400,smooth(EXPORT_LOSS.*bgcparams.Atlwidth.*bgcparams.Atllength./1e6,Nyear),'r-');
else
plot(toutmat./86400,nanmean(PROD_SURF_SUM.*bgcparams.Atlwidth.*bgcparams.Atllength./1e6).*ones(size(toutmat)),'r-',...
    toutmat./86400,nanmean(PROD_SURF_LOSS*bgcparams.Atlwidth.*bgcparams.Atllength./1e6).*ones(size(toutmat)),'r-',...
    toutmat./86400,nanmean(EXPORT_LOSS.*bgcparams.Atlwidth.*bgcparams.Atllength./1e6).*ones(size(toutmat)),'r-')
end
title('Surface nutrient consumption and export below winter MLD (kmol/s)')
legend('PROD','aPROD','(a-b)PROD')
xlabel('days since start')
%legend('D_s','D_d')
%xlim([3.6E4-720 3.6E4-360])

Dwinter=Dwinterfn(bgcparams,toutmat);
Dnow=Dfn(bgcparams.Db,(bgcparams.Dmldfrac).*Dwinter,mod(toutmat./86400,360));
ven=venfn(bgcparams,toutmat);
wen=dDdtfn_disc(bgcparams.Db,(bgcparams.Dmldfrac).*Dwinter,mod(toutmat./86400,360))...
    -ven.*Dnow./bgcparams.Atllength;
tau_lateral= (ven./(bgcparams.Atllength).*heaviside(ven) ... 
                  +bgcparams.kappah./(bgcparams.Atllength).^2 );
ENT=-wen.*heaviside(wen)./Dnow.*(squeeze(Soutmat(1,1,:)-Soutmat(1,2,:))').*max(bgcparams.Db./Dnow,1);
ENT2=-wen.*heaviside(wen)./Dnow.*(squeeze(Soutmat(1,1,:)-Soutmat(1,2,:))');
ENT3=-wen.*heaviside(wen).*(squeeze(-Soutmat(1,2,:))');
SUBDUCT3=wen.*heaviside(-wen).*(squeeze(Soutmat(1,1,:))');

MERID2=-tau_lateral.*(squeeze(Soutmat(1,1,:))'-bgcparams.Nsub);
MERID3=-tau_lateral.*(-bgcparams.Nsub);
MERID=-tau_lateral.*(-bgcparams.Nsub);
subplot(4,2,5),...
    if oneyearflag==1
  plot(toutmat./86400,86400.*ENT.*Dnow./1000.*365,...
    toutmat./86400,86400.*nanmean(ENT.*Dnow).*ones(1,length(toutmat))./1000.*365);
    else
          plot(toutmat./86400,86400.*ENT.*Dnow./1000.*365,...
    toutmat./86400,86400.*smooth(ENT.*Dnow,Nyear)./1000.*365);
    end
title('entrainment NO3 flux to D_b ~ comparable to Williams et al. 2000')
xlabel('days since start')
ylabel('mol N/m^2 yr')

subplot(4,2,6),...
    if oneyearflag == 1
  plot(toutmat./86400,86400.*MERID.*Dwinter,...
    toutmat./86400,86400.*nanmean(MERID.*Dwinter).*ones(1,length(toutmat)).*bgcparams.Atllength.*bgcparams.Atlwidth./1e6./86400);
    else
        if length(Dwinter)==1
            Dwinter=Dwinter.*ones(size(toutmat));
        end
        if length(MERID)== 1
            MERID=MERID.*ones(size(toutmat));
        end
   plot(toutmat./86400,86400.*MERID.*Dwinter.*bgcparams.Atllength.*bgcparams.Atlwidth./1e6./86400,...
    toutmat./86400,86400.*smooth(MERID.*Dwinter,Nyear).*bgcparams.Atllength.*bgcparams.Atlwidth./1e6./86400);   
    end
title('AMOC flux to Dwinter, comparable with Whitt 2019')
xlabel('days since start')
ylabel('kmol/s')

subplot(4,2,7),...
    if oneyearflag == 1
 % plot(toutmat./86400,86400.*ENT2.*Dnow,...
 %   toutmat./86400,86400.*nanmean(ENT2.*Dnow).*ones(1,length(toutmat)));
%hold on
  plot(toutmat./86400,86400.*SUBDUCT3,...
    toutmat./86400,86400.*nanmean(SUBDUCT3).*ones(1,length(toutmat)),'linewidth',1);
hold on
  plot(toutmat./86400,86400.*(ENT3),...
    toutmat./86400,86400.*nanmean(ENT3).*ones(1,length(toutmat)),'linewidth',1);
  plot(toutmat./86400,86400.*(ENT3+SUBDUCT3),'--',...
    toutmat./86400,86400.*nanmean(ENT3+SUBDUCT3).*ones(1,length(toutmat)),'--','linewidth',1);
    else    
  %plot(toutmat./86400,86400.*ENT2.*Dnow,...
 %   toutmat./86400,86400.*smooth(ENT2.*Dnow,Nyear));
%hold on
  plot(toutmat./86400,86400.*SUBDUCT3,...
    toutmat./86400,86400.*smooth(SUBDUCT3,Nyear),'linewidth',1);
hold on
  plot(toutmat./86400,86400.*(ENT3),...
    toutmat./86400,86400.*smooth(ENT3,Nyear),'linewidth',1);
  plot(toutmat./86400,86400.*(ENT3+SUBDUCT3),'--',...
    toutmat./86400,86400.*smooth(ENT3+SUBDUCT3,Nyear),'--','linewidth',1);

    end
title('total entrainment/subduction flux to D_s')
xlabel('days since start')
ylabel('mmol N/m^2 d')

subplot(4,2,8),...
    if oneyearflag == 1
%  plot(toutmat./86400,86400.*MERID2.*Dnow,...
%    toutmat./86400,86400.*nanmean(MERID2.*Dnow).*ones(1,length(toutmat)));
%hold on;
%  plot(toutmat./86400,86400.*nanmean(ENT2.*Dnow).*ones(1,length(toutmat)));
%hold on
  plot(toutmat./86400,86400.*nanmean(ENT3+SUBDUCT3).*ones(1,length(toutmat)),...
       toutmat./86400,86400.*nanmean(MERID3.*Dnow).*ones(1,length(toutmat)),...
       toutmat./86400,86400.*nanmean(MERID3.*Dnow+ENT3+SUBDUCT3).*ones(1,length(toutmat)),'--','linewidth',1);
    hold on
legend('MEAN TOTAL ENT + SUBD','MEAN MERID. FLUX at boundary','SUM SHOULD BALANCE aPROD')
    else
   plot(toutmat./86400,86400.*MERID2.*Dnow,...
    toutmat./86400,86400.*smooth(MERID2.*Dnow,Nyear));       
    end
title('total AMOC flux to D_s')
xlabel('days since start')
ylabel('mmol N/m^2 d')
%legend('D_s','D_d')
%legend('D_s','D_d')
  
subplot(4,2,3),...
plot(toutmat./86400,Dfn(bgcparams.Db,(bgcparams.Dmldfrac).*Dwinterfn(bgcparams,toutmat),mod(toutmat./86400,360)),...
    toutmat./86400,Dwinterfn(bgcparams,toutmat)-Dfn(bgcparams.Db,(bgcparams.Dmldfrac).*Dwinterfn(bgcparams,toutmat),mod(toutmat./86400,360)));
title('Depths of surface and deep layers (m)')
xlabel('days since start')
legend('D_s','D_d')
%xlim([3.6E4-720 3.6E4-360])

subplot(4,2,4),...
plot(toutmat./86400,Dfn(bgcparams.Db,(bgcparams.Dmldfrac).*Dwinterfn(bgcparams,toutmat),mod(toutmat./86400,360)).*squeeze(Soutmat(:,1,:))',...
    toutmat./86400,squeeze(Soutmat(:,2,:)).*(Dwinterfn(bgcparams,toutmat)-Dfn(bgcparams.Db,(bgcparams.Dmldfrac).*Dwinterfn(bgcparams,toutmat),mod(toutmat./86400,360)))');
title('Depth-integrals of N mmol N/m^2')
xlabel('days since start')
legend('N_s*D_1','N_d*D_2')
%xlim([3.6E4-720 3.6E4-360])
 if climatechangeflag==1
print(fig,'figures/annual_cycle_2100','-dpng')
 else
print(fig,'figures/annual_cycle','-dpng')
 end
