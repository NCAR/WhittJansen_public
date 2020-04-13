function  REAC = construct_reac_NATL_box(S,z,t,bgcparams)
N=length(z);
REAC = zeros(size(S));
    idxNs = (1:N)';
    idxNd = N+(1:N)';
   
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
            .*S(idxNs)./(bgcparams.kN+S(idxNs));
        
    Nsub = Nsubscenario(bgcparams,t);
    
    % restoring time-scale due to combined effect of AMOC advection and horizontal mixing
    tau_lateral= (ven./(bgcparams.Atllength).*heaviside(ven) ... % heaviside isn't really necessary here since ven is by construction positive 
                  +bgcparams.kappah./(bgcparams.Atllength).^2 );
    w_mix_vertical=  bgcparams.kappazML./Dnow;
        
    % N surface     % N consumption that is lost to the layer below:
    REAC(idxNs)=  - a.*prod ...
                ... % Gain/loss of N due to seasonal entrainment of deep water:
                  - wen.*heaviside(wen)./Dnow.*(S(idxNs)-S(idxNd))...
                ... % Gain/loss of N from AMOC and lateral mixing:
                  - tau_lateral.*(S(idxNs)-Nsub)...
                ... % Gain/loss of N from vert. mixing between boxes:
                  - w_mix_vertical.*(S(idxNs)-S(idxNd))./Dnow;

    % N deep        % Gains of N from remineralization of the downward sinking flux (only a fraction b converges in this box)      
    REAC(idxNd)=  + b.*Dnow./(Dwinter-Dnow).*prod ...
                 ...% Gains/Loss of N due to subduction of surface layer water:
                  - wen.*heaviside(-wen)./(Dwinter-Dnow).*(S(idxNs)-S(idxNd))...
                 ... % Gain/loss of N from AMOC and lateral mixing:
                  - tau_lateral.*(S(idxNd)-Nsub)...
                 ...% vert. mix. abyss/pycnocline exchange:
                  - bgcparams.kappaz./((Dwinter-Dnow).*(Dwinter)).*(S(idxNd)-bgcparams.Na)...
                 ...% Gain/loss of N from vert. mixing between boxes:
                  - w_mix_vertical.*(S(idxNd)-S(idxNs))./(Dwinter-Dnow);

                 
end
