function dDdt=dDdtfn(Db,Dwinter,t)
flag = 0;
% sin^2 version, too much time with deep MLDs, unrealistic:
if flag == 1
yearday=mod(t,360);
dDdt=2.*(Dwinter-Db).*sin(yearday./360.*pi+pi/4).*cos(yearday./360.*pi+pi/4).*pi./360./86400;
else
% extreme value pdf version:
y = mod(t./360.*2.*pi+pi,2*pi);
%mu=(180+50)/180*pi;
mu=(180+65)/180*pi;
sigma=.6;
scalefact=evpdfdbw(mu,mu,sigma);
dDdt=(Dwinter-Db)./scalefact.*evpdfderiv(y,mu,sigma).*2.*pi./360./86400;
end
end
