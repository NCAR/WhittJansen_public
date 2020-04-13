function Nsub = Nsubscenario(bgcparams,t)
    Nsub=bgcparams.Nsub+bgcparams.Nsubrate.*t./86400;
end
