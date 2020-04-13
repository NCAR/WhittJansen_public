function Dwinter = Dwinterfn(bgcparams,t)
if bgcparams.bndmocflag == 1
    b_n=bnscenario(bgcparams,t);
else
    if isfield(bgcparams,'bndmocscale')
    b_n=bgcparams.bn0.*ones(size(t)).*bgcparams.bndmocscale;
    else
    b_n=bgcparams.bn0.*ones(size(t));
    end
end
    b_s=b_n+0.02;
    Dwinter=bgcparams.Dwinterfrac.*D_from_b(bgcparams.h,b_s,b_n);
end

