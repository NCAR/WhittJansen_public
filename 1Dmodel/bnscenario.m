function bn = bnscenario(bgcparams,t)
    bn0=bgcparams.bn0+bgcparams.bnrate.*t./86400;
    bn=bn0+(bgcparams.bn0./2).*sin(bgcparams.omega.*t);
end
