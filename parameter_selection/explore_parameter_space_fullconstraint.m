close all;
clear all;
% parameter selection analysis script 2 - with parameter constraints
load('sensitivity_v2.mat');
load('../obs_nitrate_woa13_NATL_climatology.mat','depthN','nitrate_natl_profile');
addpath ../1Dmodel/
addpath ../utility/
addpath ../utility/export_fig/
addpath ../utility/cmocean/


Dwin=D_from_b(500,.0205,.0005);


% CONSTRAINTS:


% mean N_d target: (from fitML - case 2)
target11 = 11.33; %11.33
weight11 = 1./2;

% maximum N_s to observed max surf. NO3 (fitML -case 3):
target4=9.66; %9.66
weight4=1./2;

% minimum N_s to observed August ML NO3 (fitML -case 2):
target5=1.32;
weight5=1./2;

% PROD target:
target7 = 0.7; % 0.7 mol N/m^2/yr
weight7=1./0.7; % 0.7
% scaling factor to convert model mmol N/m^2/s -> mol N/m^2/yr
fac7=86400*365./1000;

% delta target:
target12=500; %500
weight12=1./300;

% Nsub target:
target13=19.5; %19.5
weight13=1./2.5;

% kN target:
target14=3; %3
weight14= 1./3;


delrgeg=repmat(deltarge',[1 12 12 12 12]);
Nrgeg = permute(repmat(Nrge',[1 12 12 12 12]),[2 3 4 5 1]);
kNrgeg = permute(repmat(kNrge',[1 12 12 12 12]),[2 3 4 1 5]);

sqerr=(weight11*abs(OUT(:,:,:,:,:,11)-target11)).^2 ...
     +(weight4*abs(OUT(:,:,:,:,:,4)-target4)).^2 ...
     +(weight5*abs(OUT(:,:,:,:,:,5)-target5)).^2 ...
     +(weight7*abs(OUT(:,:,:,:,:,7).*fac7-target7)).^2 ... 
     +(weight12*(abs(delrgeg-target12))).^2 ... 
     +(weight14*(abs(kNrgeg-target14))).^2 ... 
     +(weight13*(abs(Nrgeg-target13))).^2;
 
%good=sqerr<0.12; % chosen to leave 1 solutions - delta = 300
good=sqerr<0.19; % .15 chosen to leave 1 solutions - delta = 500
%good=sqerr<0.26; % chosen to leave 1 solutions - delta = 700

%good=sqerr<0.502; % chosen to leave 100 solutions

%good=sqerr<0.085; % chosen to leave 1 solutions - Nsub = 22
%good=sqerr<0.15; % chosen to leave 1 solutions - delta = 500
%good=sqerr<0.285; % chosen to leave 1 solutions - Nsub = 17
%good=sqerr<0.44; % chosen to leave 1 solutions - Nsub = 15
%good=sqerr<0.51; % chosen to leave 1 solutions - Nsub = 13

%good=sqerr<0.15; % chosen to leave 1 solutions - kN = 1,3
%good=sqerr<0.2; % chosen to leave 1 solutions - kN=5

sum(good(:))


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
set(gcf,'color','w')
colormap(cmocean('deep',32))
crge=[0 16];
Dwrge = Dwin.*Dwrge;
goodkn=good(:,:,:,3:5,:);
sum(goodkn(:))
kn=1:12;

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

PROD=365./1000.*86400.*reshape(squeeze(OUT(:,:,:,kn,:,7)),[12 12 12 12 12]);
SINK=365./1000.*86400.*reshape(squeeze(OUT(:,:,:,kn,:,8)),[12 12 12 12 12]);
goodkn=good(:,:,:,kn,:);
    PROD=PROD.*goodkn;
SINK = SINK.*goodkn;

figure;
crge=[0.5 1.0]
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







