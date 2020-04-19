% fig S6
close all;
% parameter selection analysis script 1 - without parameter constraints
clear all;
restoredefaultpath;
load('sensitivity_v2.mat');
load('../obs_nitrate_woa13_NATL_climatology.mat','depthN','nitrate_natl_profile');
addpath ../1Dmodel/
addpath ../utility/
addpath ../utility/cmocean/
Dwin=D_from_b(500,.0205,.0005);

% mean N_d target: (from fitML - case 2)
target11 = 11.33; 
weight11 = 1./2;

% maximum N_s to observed max surf. NO3 (fitML -case 3):
target4=9.66;
weight4=1./2;

% minimum N_s to observed August ML NO3 (fitML -case 2):
target5=1.32;
weight5=1./2;

% PROD target:
target7 = 0.7; % mol N/m^2/yr
weight7=1./0.7; 
% scaling factor to convert model mmol N/m^2/s -> mol N/m^2/yr
fac7=86400*365./1000;

% delta target:
target12=500;
weight12=1./300;

% Nsub target:
target13=19.5;
weight13=1./2.5;

% kN target:
target14=3;
weight14= 1./3;

delrgeg=repmat(deltarge',[1 12 12 12 12]);
Nrgeg = permute(repmat(Nrge',[1 12 12 12 12]),[2 3 4 5 1]);
kNrgeg = permute(repmat(kNrge',[1 12 12 12 12]),[2 3 4 1 5]);

sqerr=(weight11*abs(OUT(:,:,:,:,:,11)-target11)).^2 ...
     +(weight4*abs(OUT(:,:,:,:,:,4)-target4)).^2 ...
     +(weight5*abs(OUT(:,:,:,:,:,5)-target5)).^2 ...
     +(weight7*abs(OUT(:,:,:,:,:,7).*fac7-target7)).^2;
 
good=sqerr<.2653; % chosen to leave 1000 solutions
sum(good(:))
save('sensitivity_v2_good.mat','good','sqerr')


Ndelta=zeros(1,12); NDwinter=zeros(1,12);
Nmum=zeros(1,12);NkN=zeros(1,12);NN=zeros(1,12);
for ii=1:12      
   temp=good(ii,:,:,:,:);
   Ndelta(ii)=sum(temp(:));
   temp=good(:,ii,:,:,:);
   Nmum(ii)=sum(temp(:));
   temp=good(:,:,ii,:,:);
   NDwinter(ii)=sum(temp(:));
      temp=good(:,:,:,ii,:);
   NkN(ii)=sum(temp(:));
  temp=good(:,:,:,:,ii);
     NN(ii)=sum(temp(:));
end

figure;
clf
subplot(3,2,1)
bar(deltarge,Ndelta)
xlabel('\delta')
set(gca,'fontsize',18)
subplot(3,2,2)
bar(Dwrge,NDwinter)
xlabel('D_{winter}')
set(gca,'fontsize',18)
subplot(3,2,3)
bar(mumrge.*86400,Nmum)
xlabel('\mu_m')
set(gca,'fontsize',18)
subplot(3,2,4)
bar(kNrge,NkN)
xlabel('k_N')
set(gca,'fontsize',18)
subplot(3,2,5)
bar(Nrge,NN)
xlabel('N_{s}')
set(gca,'fontsize',18)
set(gcf,'color','w')


figure;
    image(kNrge,86400.*mumrge,squeeze(sum(sum(sum(good(:,:,:,:,:),1),3),5)),'CDataMapping','scaled');
xlabel('k_N')
ylabel('\mu_m')
xlim([0 9])
ylim([0 1.6])
colorbar;
caxis(caxis)
hold on
plot3(kNrge(3).*ones(3,1),linspace(0,1.6,3),1000.*ones(3,1),'color','magenta','linewidth',1)
plot3(kNrge(5).*ones(3,1),linspace(0,1.6,3),1000.*ones(3,1),'color','magenta','linewidth',1)
set(gcf,'color','w')
colormap(cmocean('deep',32))
set(gca,'fontsize',18)

crge=[0 16];
Dwrge = Dwin.*Dwrge;
goodkn=good(:,:,:,3:5,:);
sum(goodkn(:))
kn=3:5
%for kn = 3:3
figure
subplot(3,2,1),...
    image(Dwrge,deltarge,squeeze(sum(sum(sum(good(:,:,:,kn,:),2),4),5)),'CDataMapping','scaled');
xlabel('D_w')
ylabel('\delta')
xlim([0 1000])
ylim([0 800])
colorbar;

subplot(3,2,2),...
    image(Nrge,deltarge,squeeze(sum(sum(sum(good(:,:,:,kn,:),2),3),4)),'CDataMapping','scaled');
xlabel('N_s')
ylabel('\delta')
colorbar;
xlim([11 21])
ylim([0 800])

subplot(3,2,3),...
    image(mumrge.*86400,deltarge,squeeze(sum(sum(sum(good(:,:,:,kn,:),5),3),4)),'CDataMapping','scaled');
xlabel('\mu_m')
ylabel('\delta')
colorbar;
xlim([0 1.6])
ylim([0 800])

subplot(3,2,4),...
    image(Dwrge,mumrge.*86400,squeeze(sum(sum(sum(good(:,:,:,kn,:),1),4),5)),'CDataMapping','scaled');
xlabel('D_w')
ylabel('\mu_m')
colorbar;
xlim([0 1000])
ylim([0 1.6])

subplot(3,2,5),...
    image(Nrge,Dwrge,squeeze(sum(sum(sum(good(:,:,:,kn,:),1),4),2)),'CDataMapping','scaled');
xlabel('N_s')
ylabel('D_w')
colorbar;
xlim([11 21])
ylim([0 1000])

subplot(3,2,6),...
    image(Nrge,86400.*mumrge,squeeze(sum(sum(sum(good(:,:,:,kn,:),1),3),4)),'CDataMapping','scaled');
xlabel('N_s')
ylabel('\mu_n')
colorbar;
xlim([11 21])
ylim([0 1.6])
colormap(cmocean('deep',32))
set(gcf,'color','w')

kn=3:5
PROD=365./1000.*86400.*reshape(squeeze(OUT(:,:,:,kn,:,7)),[12 12 12 3 12]);
SINK=365./1000.*86400.*reshape(squeeze(OUT(:,:,:,kn,:,8)),[12 12 12 3 12]);
goodkn=good(:,:,:,kn,:);
    PROD=PROD.*goodkn;
SINK = SINK.*goodkn;

figure;
crge=[0.1 1.4]
subplot(3,2,1),...
    image(Dwrge,deltarge,squeeze(sum(sum(sum(PROD,2),4),5)./sum(sum(sum(goodkn,2),4),5)),'CDataMapping','scaled');
xlabel('D_w')
ylabel('\delta')
xlim([0 1000])
ylim([0 800])
h=colorbar;
ylabel(h,'mol N/m^2/yr')
caxis(crge)

subplot(3,2,2),...
    image(Nrge,deltarge,squeeze(sum(sum(sum(PROD,2),3),4)./sum(sum(sum(goodkn,2),3),4)),'CDataMapping','scaled');
xlabel('N_s')
ylabel('\delta')
h=colorbar;
ylabel(h,'mol N/m^2/yr')
xlim([11 21])
ylim([0 800])
caxis(crge)

subplot(3,2,3),...
    image(mumrge.*86400,deltarge,squeeze(sum(sum(sum(PROD,3),4),5)./sum(sum(sum(goodkn,3),4),5)),'CDataMapping','scaled');
xlabel('\mu_m')
ylabel('\delta')
h=colorbar;
ylabel(h,'mol N/m^2/yr')
xlim([0 1.6])
ylim([0 800])
caxis(crge)

subplot(3,2,4),...
    image(Dwrge,mumrge.*86400,squeeze(sum(sum(sum(PROD,1),4),5)./sum(sum(sum(goodkn,1),4),5)),'CDataMapping','scaled');
xlabel('D_w')
ylabel('\mu_m')
h=colorbar;
ylabel(h,'mol N/m^2/yr')
xlim([0 1000])
ylim([0 1.6])
caxis(crge)

subplot(3,2,5),...
    image(Nrge,Dwrge,squeeze(sum(sum(sum(PROD,2),4),1)./sum(sum(sum(goodkn,2),4),1)),'CDataMapping','scaled');
xlabel('N_s')
ylabel('D_w')
h=colorbar;
ylabel(h,'mol N/m^2/yr')
xlim([11 21])
ylim([0 1000])
caxis(crge)

subplot(3,2,6),...
    image(Nrge,86400.*mumrge,squeeze(sum(sum(sum(PROD,1),4),3)./sum(sum(sum(goodkn,1),4),3)),'CDataMapping','scaled');
xlabel('N_s')
ylabel('\mu_n')
h=colorbar;
ylabel(h,'mol N/m^2/yr')
xlim([11 21])
ylim([0 1.6])
caxis(crge)
colormap(cmocean('deep',32))
set(gcf,'color','w')


crge= [0 1];
figure;
subplot(3,2,1),...
    image(Dwrge,deltarge,squeeze(nansum(nansum(nansum(SINK./PROD,2),4),5)./nansum(nansum(nansum(goodkn,2),4),5)),'CDataMapping','scaled');
xlabel('D_w')
ylabel('\delta')
xlim([0 1000])
ylim([0 800])
colorbar;

caxis(crge)

subplot(3,2,2),...
    image(Nrge,deltarge,squeeze(nansum(nansum(nansum(SINK./PROD,2),3),4)./nansum(nansum(nansum(goodkn,2),3),4)),'CDataMapping','scaled');
xlabel('N_s')
ylabel('\delta')
colorbar;

xlim([11 21])
ylim([0 800])
caxis(crge)

subplot(3,2,3),...
    image(mumrge.*86400,deltarge,squeeze(nansum(nansum(nansum(SINK./PROD,3),4),5)./nansum(nansum(nansum(goodkn,3),4),5)),'CDataMapping','scaled');
xlabel('\mu_m')
ylabel('\delta')
colorbar;

xlim([0 1.6])
ylim([0 800])
caxis(crge)

subplot(3,2,4),...
    image(Dwrge,mumrge.*86400,squeeze(nansum(nansum(nansum(SINK./PROD,1),4),5)./nansum(nansum(nansum(goodkn,1),4),5)),'CDataMapping','scaled');
xlabel('D_w')
ylabel('\mu_m')
colorbar;

xlim([0 1000])
ylim([0 1.6])
caxis(crge)

subplot(3,2,5),...
    image(Nrge,Dwrge,squeeze(nansum(nansum(nansum(SINK./PROD,2),4),1)./nansum(nansum(nansum(goodkn,2),4),1)),'CDataMapping','scaled');
xlabel('N_s')
ylabel('D_w')
colorbar;

xlim([11 21])
ylim([0 1000])
caxis(crge)

subplot(3,2,6),...
    image(Nrge,86400.*mumrge,squeeze(nansum(nansum(nansum(SINK./PROD,1),4),3)./nansum(nansum(nansum(goodkn,1),4),3)),'CDataMapping','scaled');
xlabel('N_s')
ylabel('\mu_n')
colorbar;

xlim([11 21])
ylim([0 1.6])
caxis(crge)
colormap(cmocean('deep',32))
set(gcf,'color','w')

Dwodelta=zeros(size(SINK));
for ik = 1:12
    for ij = 1:12
Dwodelta(ij,:,ik,:,:)=Dwrge(ik)./deltarge(ij);
    end
end


Ng=zeros(size(SINK));
for ik = 1:12
    Ng(:,:,:,:,ik)=Nrge(ik);
end

mug=zeros(size(SINK));
for ik = 1:12
    mug(:,ik,:,:,:)=mumrge(ik);
end


dg=zeros(size(SINK));
for ik = 1:12
    dg(ik,:,:,:,:)=deltarge(ik);
end


%%
figure;
Dwodelta=Dwodelta.*goodkn;
scatter(Dwodelta(:),SINK(:)./PROD(:),28,Ng(:),'filled');
hold on
plot(linspace(.1,4,41),exp(-linspace(.1,4,41)),'k--','linewidth',1)
caxis([11 23])
colormap(cmocean('deep',32));
cbh=colorbar;
ylabel(cbh,'N_{sub}','fontsize',18)

A=1./Dwodelta(:);
B=SINK(:)./PROD(:);
[r,p]=corrcoef(A(~isnan(A+B)),B(~isnan(A+B)))

p=polyfit(A(~isnan(A+B)),B(~isnan(A+B)),1);
hold on
%plot(linspace(0,2,11),polyval(p,linspace(0,2,11)),'r-');
ylim([0 1])
xlim([0 4])
grid on;
xlabel('D_w/\delta')
ylabel('a-b');
set(gcf,'color','w')
set(gca,'fontsize',18)

figure;
scatter(Dwodelta(PROD~=0),Ng(PROD~=0),28,PROD(PROD~=0),'filled');
caxis([.7 1])
colormap(cmocean('deep',32));
cbh=colorbar;
ylabel(cbh,'mol N/m^2/yr','fontsize',18)

xlim([0 4])
ylim([11 23])
grid on;
xlabel('D_w/\delta')
ylabel('N_{sub}');
set(gcf,'color','w')
set(gca,'fontsize',18)


%%

figure;
subplot(2,2,1),...
scatter(1./Dwodelta(:),SINK(:)./PROD(:),15,Ng(:),'filled');
caxis([12 22])
colormap(cmocean('deep',32));
cbh=colorbar;
ylabel(cbh,'N_s','fontsize',12)

A=1./Dwodelta(:);
B=SINK(:)./PROD(:);
[r,p]=corrcoef(A(~isnan(A+B)),B(~isnan(A+B)))

p=polyfit(A(~isnan(A+B)),B(~isnan(A+B)),1);
hold on
%plot(linspace(0,1,11),polyval(p,linspace(0,1,11)),'r-');
ylim([0 1])
xlim([0 2])
grid on;
xlabel('\delta/D_w')
ylabel('a-b');

subplot(2,2,2),...
scatter(Ng(:),SINK(:)./PROD(:),15,1./Dwodelta(:),'filled');
caxis([0 2])
colormap(cmocean('deep',32))
cbh=colorbar;
ylabel(cbh,'\delta/D_w','fontsize',12)

A=Ng(:);
B=SINK(:)./PROD(:);
[r,p]=corrcoef(A(~isnan(A+B)),B(~isnan(A+B)))

p=polyfit(A(~isnan(A+B)),B(~isnan(A+B)),1);

hold on
plot(linspace(12,22,11),polyval(p,linspace(12,22,11)),'r-');
ylim([0 1])
xlim([12 22])
grid on;
xlabel('N_s')
ylabel('a-b');


subplot(2,2,3),...
scatter(Ng(:),1./Dwodelta(:),15,SINK(:)./PROD(:),'filled');
caxis([0 0.4])
colormap(cmocean('deep',32))
A1=Ng(~(isnan(A+B)));
B=1./Dwodelta(~(isnan(A+B)));
A=A1;
[r,p]=corrcoef(A(~isnan(A+B)),B(~isnan(A+B)))

p=polyfit(A(~isnan(A+B)),B(~isnan(A+B)),1);

hold on
plot(linspace(12,20,11),polyval(p,linspace(12,20,11)),'r-');
ylim([0 2])
xlim([12 22])
grid on;
xlabel('N_s')
ylabel('\delta/D_w');
cbh=colorbar;
ylabel(cbh,'export ratio','fontsize',12)
set(gcf,'color','w')

subplot(2,2,4),...
scatter(mug(:).*86400,PROD(:)./365.*1000,'k.');
A=mug(:).*86400;
B=PROD(:)./365.*1000;
B(B==0)=nan;
A(A==0)=nan;
[r,p]=corrcoef(A(~isnan(A+B)),B(~isnan(A+B)))

p=polyfit(A(~isnan(A+B)),B(~isnan(A+B)),1);

hold on
%plot(linspace(.4,1,11),polyval(p,linspace(.4,1,11)),'r-');
ylim([.1 1.4].*3)
xlim([0 1])
grid on;
xlabel('\mu_m mmol/m^3/d')
ylabel('PROD mmol/m^2/d');






