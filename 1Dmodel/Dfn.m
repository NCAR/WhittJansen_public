function D=Dfn(Db,Dwinter,t)
flag = 0;
if flag == 1
% sin^2 version, too much time with deep MLDs, unrealistic:
yearday=mod(t,360);
D=Db+(Dwinter-Db).*sin(yearday./360.*pi+pi/4).^2;
else
% MLD seasonal cycle from Extreme Value Distribution function
y = mod(t./360.*2.*pi+pi,2*pi);
mu=(180+65)/180*pi;
sigma=.6;
scalefact=evpdfdbw(mu,mu,sigma);
D = (Dwinter-Db)./scalefact.*(evpdfdbw(y,mu,sigma))+Db;
end
end
