% Example script for time evolution of ML depth and Psi through global waming

bgcparams=set_bgcparams_fn;
%100year daily time series for surface buoyancy in the north, assuming linear
%"warming"
b_n=linspace(0.0005,0.005,36500);
b_s=b_n+0.02;
D=0*b_n;
Psi=0*b_n;
for ii=1:length(b_n)
    dayofyear=mod(ii,365);
    Dwinter=D_from_b(bgcparams.h,b_s(ii),b_n(ii));
    D(ii)=Dfn(bgcparams.Db,Dwinter,dayofyear);
    Psi(ii)=Psimax(bgcparams.h,b_s(ii),b_n(ii));
end

figure(1);clf;
subplot(2,2,1)
plot(b_n,'r');hold on
subplot(2,2,2)
plot(Psi/1e6,'b');hold on
subplot(2,1,2)
plot(D,'b');hold on
