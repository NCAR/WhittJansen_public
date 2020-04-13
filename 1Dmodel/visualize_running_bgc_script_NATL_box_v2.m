clf; cla;
Nyear=round(3.1104e7/(toutmat(2)-toutmat(1)));
aplt=Dfn(bgcparams.Db,(bgcparams.Dmldfrac).*Dwinterfn(bgcparams,toutmat(1:ct)),mod(toutmat(1:ct)./86400,360)).*squeeze(Soutmat(:,1,1:ct))';
bplt=squeeze(Soutmat(:,2,1:ct)).*(Dwinterfn(bgcparams,toutmat(1:ct))-Dfn(bgcparams.Db,(bgcparams.Dmldfrac).*Dwinterfn(bgcparams,toutmat(1:ct)),mod(toutmat(1:ct)./86400,360)))';
dintNdt=zeros(size(bplt));
dintNdt(2:end-1)=aplt(3:ct)'+bplt(3:ct) - aplt(1:(ct-2))'-bplt(1:(ct-2));
dintNdt=dintNdt';
dintNdt(2:end-1)=dintNdt(2:end-1)./(toutmat(3:ct)-toutmat(1:(ct-2)));
dintNdt(1)=dintNdt(2);
dintNdt(end)=dintNdt(end-1);

subplot(4,1,1),...
plot(toutmat(1:ct)./86400./360,smooth(aplt'+bplt,Nyear)./1000);
title('Depth-integral of nutrient in both layers')
xlim([1 (yrs-1)])
%legend('N_s*D_s+D_d(D_w-D_s)')
ylim([0 4])
ylabel('(mol /m^2)')
set(gca,'fontsize',15);
grid on

% subplot(4,1,1),...
% plot(toutmat(1:ct)./86400./360,smooth((aplt'+bplt)./((Dwinterfn(bgcparams,toutmat(1:ct)))'),Nyear));
% title('Average N concentration, including both two boxes (mmol N/m^3)')
% xlabel('years since start')
% xlim([1 (yrs-1)])
% %legend('')
% ylim([0 11])

PROD_SURF_SUM=bgcparams.mum./(bgcparams.kz).*...
                      log((bgcparams.kI+Ifn(bgcparams.latlight,mod(toutmat(1:ct)./86400,360)))./bgcparams.kI).* ...
                     (squeeze(Soutmat(:,1,1:ct))./(bgcparams.kN+squeeze(Soutmat(:,1,1:ct))))';
                 EXPORT_LOSS=((aflxfn(Dfn(bgcparams.Db,(bgcparams.Dmldfrac).*Dwinterfn(bgcparams,toutmat(1:ct)),mod(toutmat(1:ct)./86400,360)),bgcparams.delta))...
                     .*bgcparams.mum./(bgcparams.kz).*...
                     log((bgcparams.kI+Ifn(bgcparams.latlight,mod(toutmat(1:ct)./86400,360)))./bgcparams.kI).* ...
                     (squeeze(Soutmat(:,1,1:ct))./(bgcparams.kN+squeeze(Soutmat(:,1,1:ct))))') - ...
                     bflxfn(Dfn(bgcparams.Db,(bgcparams.Dmldfrac).*Dwinterfn(bgcparams,toutmat(1:ct)),mod(toutmat(1:ct)./86400,360)),(Dwinterfn(bgcparams,toutmat(1:ct))-Dfn(bgcparams.Db,(bgcparams.Dmldfrac).*Dwinterfn(bgcparams,toutmat(1:ct)),mod(toutmat(1:ct)./86400,360))),bgcparams.delta).*...
                     bgcparams.mum./(bgcparams.kz).*...
                     log((bgcparams.kI+Ifn(bgcparams.latlight,mod(toutmat(1:ct)./86400,360)))./bgcparams.kI).* ...
                     (squeeze(Soutmat(:,1,1:ct))./(bgcparams.kN+squeeze(Soutmat(:,1,1:ct))))';

Dwinter=Dwinterfn(bgcparams,toutmat(1:ct));
Dnow=Dfn(bgcparams.Db,(bgcparams.Dmldfrac).*Dwinter,mod(toutmat(1:ct)./86400,360));
ven=venfn(bgcparams,toutmat(1:ct));
tau_lateral= (ven./(bgcparams.Atllength).*heaviside(ven) ... 
                  +bgcparams.kappah./(bgcparams.Atllength).^2);
              
wen=dDdtfn_disc(bgcparams.Db,(bgcparams.Dmldfrac).*Dwinter,mod(toutmat(1:ct)./86400,360))...
    -ven.*Dnow./bgcparams.Atllength;
w=-ven.*Dnow./bgcparams.Atllength;
ENTs=-wen.*heaviside(wen).*(squeeze(Soutmat(1,1,1:ct)-Soutmat(1,2,1:ct))');
ENTd=-wen.*heaviside(-wen).*(squeeze(Soutmat(1,1,1:ct)-Soutmat(1,2,1:ct))');
SUBtend=-w.*(squeeze(Soutmat(1,1,1:ct)-Soutmat(1,2,1:ct))');
SUBtrans=zeros(size(SUBtend));
SUBtrans(2:end-1)=(Dwinter(3:ct)-Dwinter(1:(ct-2)))./(toutmat(3:ct)-toutmat(1:(ct-2))).*(squeeze(Soutmat(1,2,2:(ct-1)))');
SUBtrans(1)=SUBtrans(2);
SUBtrans(end)=SUBtrans(end-1);
%ENT3=-wen.*heaviside(wen).*(squeeze(-Soutmat(1,2,:))');
%SUBDUCT3=wen.*heaviside(-wen).*(squeeze(Soutmat(1,1,:))');

MERIDs=-tau_lateral.*(squeeze(Soutmat(1,1,1:ct))'-bgcparams.Nsub).*Dnow;
MERIDd=-tau_lateral.*(squeeze(Soutmat(1,2,1:ct))'-bgcparams.Nsub).*(Dwinter-Dnow);

VMIX= -bgcparams.kappaz./((Dwinter)).*(squeeze(Soutmat(1,2,1:ct))'-bgcparams.Na);

 subplot(4,1,[3 4]),...
plot(toutmat(1:ct)./86400./360,smooth(360./1000.*86400.*(MERIDs+MERIDd+SUBtend-nanmean(MERIDs(1:noutperyr)+MERIDd(1:noutperyr)+SUBtend(1:noutperyr))),Nyear)',...
toutmat(1:ct)./86400./360,-smooth(360./1000.*86400.*(EXPORT_LOSS-nanmean(EXPORT_LOSS(1:noutperyr))),Nyear)',...
toutmat(1:ct)./86400./360,smooth(360./1000.*86400.*(VMIX),Nyear)','--',...
toutmat(1:ct)./86400./360,smooth(360./1000.*86400.*SUBtrans,Nyear)','--','linewidth',2);
hold on;
plot(toutmat(1:ct)./86400./360,smooth(360./1000.*86400.*(MERIDs+MERIDd+SUBtend-EXPORT_LOSS),Nyear)','--','linewidth',2)
%plot(toutmat(1:ct)./86400./360,smooth(360./1000.*86400.*(MERIDs+MERIDd+SUBtend-EXPORT_LOSS+VMIX+SUBtrans),Nyear)','k-','linewidth',2);
plot(toutmat(1:ct)./86400./360,smooth(360./1000.*86400.*dintNdt,Nyear)','--','linewidth',2);
legend;
legend('ADV pert.','EXPORT pert.','VMIX','SUBDUCT','ADV-EXPORT','dintN/dt','location','northwest')
title('Terms in the budget of depth-integrated N, including both layers')
xlabel('years since start')
xlim([1 (yrs-1)])
%legend('ADV.','VMIX','EXPORT')
%ylim([-.0001 .0001])
ylim([-.175 .175])
grid on;
ylabel('(mol /m^2 per year)')
set(gca,'fontsize',15);

subplot(4,1,2),...
plot(toutmat(1:ct)./86400./360,Dwinterfn(bgcparams,toutmat(1:ct)));
title('Winter mixed layer depth')
xlim([1 (yrs-1)])
ylim([0 400])
ylabel('m')
set(gca,'fontsize',15);
grid on

