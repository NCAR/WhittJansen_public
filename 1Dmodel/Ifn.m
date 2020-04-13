function I = Ifn(lat,yearday)
Aatm = .5; % low fraction of incoming solar gets to the surface 
%- north atlantic is cloudy! 
PARfrac = .4;
S0=1367;
e = .0167;
eps1 = 2.*pi.*23./360; % obliquity
pibar = 2.*pi.*282.89./360;
% Rfac is always within 3% of 1 here
Rfac = (1+e.*cos(thetafn(yearday)-pibar)).^2;
% so don't use:
Rfac=1;
hrhs = -tan(2.*pi.*lat./360).*tan(eps1.*sin(thetafn(yearday)));
h0 = acos(hrhs);
mask1 = hrhs < -1;
h0(mask1) = pi;
mask2 = hrhs>1;
h0(mask2) = 0;
I = Aatm.*PARfrac.*S0./pi.*Rfac.*(h0.*sind(lat).*sin(eps1.*sin(thetafn(yearday)))...
    +cosd(lat).*cos(eps1.*sin(thetafn(yearday))).*sin(h0));

end

function theta = thetafn(yearday)
theta = 2.*pi.*(yearday-80)./360;
end