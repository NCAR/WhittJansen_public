function  [REAC,CONC,INT,Dnow,Dwinter,wen,a,b,prod] = construct_reac_NATL_box_budget(S,z,t,bgcparams)
N=length(z);
idxNs = (1:N)';
idxNd = N+(1:N)';
S=squeeze(S);
REAC = zeros(size(S));
Dwinter=Dwinterfn(bgcparams,t);
Dnow=Dfn(bgcparams.Db,(bgcparams.Dmldfrac).*Dwinter,mod(t./86400,360));
ven=venfn(bgcparams,t);
wen=dDdtfn_disc(bgcparams.Db,(bgcparams.Dmldfrac).*Dwinter,mod(t./86400,360))...
    -ven.*Dnow./bgcparams.Atllength;
I=Ifn(bgcparams.latlight,mod(t./86400,360));
a=aflxfn(Dnow,bgcparams.delta);
b=bflxfn(Dnow,(Dwinter-Dnow),bgcparams.delta);
% N consumption in the ML
prod=bgcparams.mum./(Dnow.*bgcparams.kz)...
    .*log((bgcparams.kI+I)./bgcparams.kI)...
    .*S(idxNs,:)./(bgcparams.kN+S(idxNs,:));

% restoring time-scale due to combined effect of AMOC advection and horizontal mixing
tau_lateral= (ven./(bgcparams.Atllength).*heaviside(ven) ... % heaviside isn't really necessary here since ven is by construction positive
    +bgcparams.kappah./(bgcparams.Atllength).^2 );
w_mix_vertical=  bgcparams.kappazML./Dnow;

% N surface     % N consumption that is lost to the layer below:
REAC(idxNs,:)=  - a.*prod ...
    ... % Gain/loss of N due to seasonal entrainment of deep water:
    - wen.*heaviside(wen)./Dnow.*(S(idxNs,:)-S(idxNd,:))...
    ... % Gain/loss of N from AMOC and lateral mixing:
    - tau_lateral.*(S(idxNs,:)-bgcparams.Nsub)...
    ... % Gain/loss of N from vert. mixing between boxes:
    - w_mix_vertical.*(S(idxNs,:)-S(idxNd,:))./Dnow;
% concentration budget:
dNsdt=zeros(size(tau_lateral));
dNsdt(2:end-1)=S(idxNs,3:end) - S(idxNs,1:end-2);
dNsdt(2:end-1)=dNsdt(2:end-1)./(t(3:end)-t(1:(end-2)));
dNsdt(1)=dNsdt(2);
dNsdt(end)=dNsdt(end-1);

CONC{1}.PROD=- a.*prod;
CONC{1}.ENT =- wen.*heaviside(wen)./Dnow.*(S(idxNs,:)-S(idxNd,:));
CONC{1}.LAT =- tau_lateral.*(S(idxNs,:)-bgcparams.Nsub);
CONC{1}.LATADV =-(ven./(bgcparams.Atllength).*heaviside(ven)).*(S(idxNs,:)-bgcparams.Nsub);
CONC{1}.LATMIX =-(bgcparams.kappah./(bgcparams.Atllength).^2).*(S(idxNs,:)-bgcparams.Nsub);
CONC{1}.MIX =- w_mix_vertical.*(S(idxNs,:)-S(idxNd,:))./Dnow;
CONC{1}.REAC= REAC(idxNs,:);
CONC{1}.RATE= dNsdt;

% integrated budget:
Nsint=Dnow.*S(idxNs,:);

dintNsdt=zeros(size(tau_lateral));
dintNsdt(2:end-1)=Nsint(3:end) - Nsint(1:end-2);
dintNsdt(2:end-1)=dintNsdt(2:end-1)./(t(3:end)-t(1:(end-2)));
dintNsdt(1)=dintNsdt(2);
dintNsdt(end)=dintNsdt(end-1);

NsdDdt=zeros(size(tau_lateral));
NsdDdt(2:end-1)=S(idxNs,2:end-1).*(Dnow(3:end)-Dnow(1:end-2))./(t(3:end)-t(1:end-2));
NsdDdt(1)=NsdDdt(2);
NsdDdt(end)=NsdDdt(end-1);

INT{1}.PROD=CONC{1}.PROD.*Dnow;
INT{1}.ENT=CONC{1}.ENT.*Dnow + NsdDdt;
INT{1}.LAT=CONC{1}.LAT.*Dnow;
INT{1}.LATADV=CONC{1}.LATADV.*Dnow;
INT{1}.LATMIX=CONC{1}.LATMIX.*Dnow;
INT{1}.MIX=CONC{1}.MIX.*Dnow;
INT{1}.REAC=CONC{1}.REAC.*Dnow + NsdDdt;
INT{1}.RATE=dintNsdt;

% N deep        % Gains of N from remineralization of the downward sinking flux (only a fraction b converges in this box)
REAC(idxNd,:)=  + b.*Dnow./(Dwinter-Dnow).*prod ...
    ...% Gains/Loss of N due to subduction of surface layer water:
    - wen.*heaviside(-wen)./(Dwinter-Dnow).*(S(idxNs,:)-S(idxNd,:))...
    ... % exchange of N with surface layer due to vertical mixing
    - tau_lateral.*(S(idxNd,:)-bgcparams.Nsub)...
    ...% vert. mix. abyss/pycnocline exchange:
    - bgcparams.kappaz./((Dwinter-Dnow).*(Dwinter)).*(S(idxNd,:)-bgcparams.Na)...
    ...% Gain/loss of N from vert. mixing between boxes:
    - w_mix_vertical.*(S(idxNd,:)-S(idxNs,:))./(Dwinter-Dnow);

% concentration budget:
dNddt=zeros(size(tau_lateral));
dNddt(2:end-1)=S(idxNd,3:end) - S(idxNd,1:end-2);
dNddt(2:end-1)=dNddt(2:end-1)./(t(3:end)-t(1:(end-2)));
dNddt(1)=dNddt(2);
dNddt(end)=dNddt(end-1);

CONC{2}.PROD=b.*Dnow./(Dwinter-Dnow).*prod;
CONC{2}.ENT =- wen.*heaviside(-wen)./(Dwinter-Dnow).*(S(idxNs,:)-S(idxNd,:));
CONC{2}.LAT =- tau_lateral.*(S(idxNd,:)-bgcparams.Nsub);
CONC{2}.LATADV =-(ven./(bgcparams.Atllength)).*(S(idxNd,:)-bgcparams.Nsub);
CONC{2}.LATMIX =-(bgcparams.kappah./(bgcparams.Atllength).^2).*(S(idxNd,:)-bgcparams.Nsub);
CONC{2}.MIX =- bgcparams.kappaz./((Dwinter-Dnow).*(Dwinter)).*(S(idxNd,:)-bgcparams.Na)...
    - w_mix_vertical.*(S(idxNd,:)-S(idxNs,:))./(Dwinter-Dnow);
CONC{2}.REAC= REAC(idxNd,:);
CONC{2}.RATE= dNddt;

% integrated budget:
Ndint=(Dwinter-Dnow).*S(idxNd,:);
dintNddt=zeros(size(tau_lateral));
dintNddt(2:end-1)=Ndint(3:end) - Ndint(1:end-2);
dintNddt(2:end-1)=dintNddt(2:end-1)./(t(3:end)-t(1:(end-2)));
dintNddt(1)=dintNddt(2);
dintNddt(end)=dintNddt(end-1);
NdDddt=zeros(size(tau_lateral));
NdDddt(2:end-1)=(Dwinter(3:end)-Dnow(3:end)-Dwinter(1:end-2)+Dnow(1:end-2))./(t(3:end)-t(1:end-2)).*S(idxNd,2:end-1);
NdDddt(1)=NdDddt(2);
NdDddt(end)=NdDddt(end-1);

INT{2}.PROD=CONC{2}.PROD.*(Dwinter-Dnow);
INT{2}.ENT=CONC{2}.ENT.*(Dwinter-Dnow) + NdDddt;
INT{2}.LAT=CONC{2}.LAT.*(Dwinter-Dnow);
INT{2}.LATADV=CONC{2}.LATADV.*(Dwinter-Dnow);
INT{2}.LATMIX=CONC{2}.LATMIX.*(Dwinter-Dnow);
INT{2}.MIX=CONC{2}.MIX.*(Dwinter-Dnow);
INT{2}.REAC=CONC{2}.REAC.*(Dwinter-Dnow)+ NdDddt;
INT{2}.RATE=dintNddt;

%% both boxes:

% integrated budget:
dintNdt=zeros(size(tau_lateral));
dintNdt(2:end-1)=Nsint(3:end)+ Ndint(3:end) - Nsint(1:end-2)- Ndint(1:end-2);
dintNdt(2:end-1)=dintNdt(2:end-1)./(t(3:end)-t(1:(end-2)));
dintNdt(1)=dintNdt(2);
dintNdt(end)=dintNdt(end-1);
NddDwdt=zeros(size(tau_lateral));

NddDwdt(2:end-1)=S(idxNd,2:end-1).*(Dwinter(3:end)-Dwinter(1:end-2))./(t(3:end)-t(1:end-2));
NddDwdt(end)=NddDwdt(end-1);
NddDwdt(1)=NddDwdt(2);

INT{3}.PROD=(b-a).*prod.*Dnow;
INT{3}.MIX=(-bgcparams.kappaz./(Dwinter).*(S(idxNd,:)-bgcparams.Na));
INT{3}.LAT=ven.*Dwinter./bgcparams.Atllength.*(bgcparams.Nsub-S(idxNd,:))...
    - (S(idxNd,:)-bgcparams.Nsub).*Dwinter.*bgcparams.kappah./(bgcparams.Atllength).^2 ...
    - (S(idxNs,:)-S(idxNd,:)).*Dnow.*bgcparams.kappah./(bgcparams.Atllength).^2;
INT{3}.LATADV=ven.*Dwinter./bgcparams.Atllength.*(bgcparams.Nsub-S(idxNd,:));
INT{3}.LATMIX=-(S(idxNd,:)-bgcparams.Nsub).*Dwinter.*bgcparams.kappah./(bgcparams.Atllength).^2 ...
    - (S(idxNs,:)-S(idxNd,:)).*Dnow.*bgcparams.kappah./(bgcparams.Atllength).^2;
INT{3}.SUBD=NddDwdt;
INT{3}.RATE=dintNdt;
INT{3}.REAC=REAC(1,:).*Dnow+REAC(2,:).*(Dwinter-Dnow)+NdDddt+NsdDdt;

% avg conc:
davgNdt=zeros(size(tau_lateral));
davgNdt(2:end-1)=((Nsint(3:end)+ Ndint(3:end))./Dwinter(3:end) - (Nsint(1:end-2)+ Ndint(1:end-2))./Dwinter(1:end-2));
davgNdt(2:end-1)=davgNdt(2:end-1)./(t(3:end)-t(1:(end-2)));
davgNdt(1)=davgNdt(2);
davgNdt(end)=davgNdt(end-1);

NIdonDwdt=zeros(size(tau_lateral));
NIdonDwdt(2:end-1)=(Nsint(2:end-1)+Ndint(2:end-1)).*(1./Dwinter(3:end)-1./Dwinter(1:end-2))./(t(3:end)-t(1:end-2));
NIdonDwdt(1)=NIdonDwdt(2);
NIdonDwdt(end)=NIdonDwdt(end-1);

CONC{3}.PROD=INT{3}.PROD./Dwinter;
CONC{3}.MIX=INT{3}.MIX./Dwinter;
CONC{3}.LAT=INT{3}.LAT./Dwinter;
CONC{3}.LATADV=INT{3}.LATADV./Dwinter;
CONC{3}.LATMIX=INT{3}.LATMIX./Dwinter;
CONC{3}.SUBD=INT{3}.SUBD./Dwinter + NIdonDwdt;
CONC{3}.RATE=davgNdt;

end
