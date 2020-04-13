function ven = venfn(bgcparams,t)
if bgcparams.bnpsiflag==1
    b_n=bnscenario(bgcparams,t);
else
    b_n=bgcparams.bn0.*ones(size(t));
end
    b_s=b_n+0.02;
    Psi=Psimax(bgcparams.h,b_s,b_n);
    ven=Psi./bgcparams.Atlwidth./Dwinterfn(bgcparams,t);
end

