function bgcparams=set_bgcparams_fn()
% script to set bgcparams

%% N loss parameters
bgcparams.mum = 0.48./86400; % maximum nutrient uptake rate mmol N /m^3 per day (assuming nitrogen units). Not well constrained. McKinley et al. (2004) use 1.06 
% in the subpolar N. Atlantic. probably could reasonably be somewhere between .2 and 2
bgcparams.kN = 3.2; % mmol N /m^3 nutrient half-saturation; nutrient uptake rate declines for N<~ 1; Dutkiewicz et al. 2001 DSR take it to be 1; not well constrained
bgcparams.kz = 1./10; % 1/m attenuation of light with depth; Jerlov type II/III water (Simonot and Le Treut 1986 JGR Climatology, Paulson and Simpson 1977) -constrained within a factor of 2
bgcparams.kI = 30; % W/m^2 light half-saturation; nutrient uptake rate declines for irradiances near 30 W/m^2 - coincidentally the same as McKinley et al. (2004) - not that well constrained; perhaps could be 5-100?
% & Dutkiewicz et al. 2001 DSR
bgcparams.latlight = 60;  % target latitude for the irradiance/light seasonal cycle


%% Physical parameters
bgcparams.bn0=0.0005; % (initial) surface buoyancy in the north
bgcparams.omega=0;
bgcparams.bnrate=0;
bgcparams.bnpsiflag =1; % set to zero to hold psi at psi(0.0005)
bgcparams.bndmocflag = 1; % set to zero to hold dmoc at dmoc(.0005)
bgcparams.bndmocscale = 1; % scale up or down bn to calculate Dwinter depth
bgcparams.h=500; % depth scale of the pycnocline in the tropics

bgcparams.Db = 40; % shallowest MLD/biologically active layer in summer; assumed deep relative to 1/kz  
% so that no nutrient consumption occurs below the surface layer
% looks consistent with obs, but some nutrient consumption may occur below the MLD during summer,
% perhaps to a depth between 50 and 75 m

% vertical diffusivity in the pycnocline:
bgcparams.kappaz=2e-5;%2.5e-3; 
%bgcparams.kappaz=0;
% the depth scale for vertical mixing is defined to be the winter mixed layer depth ~10^3 m 
% for a 400m winter MLD 1e-4 yields a restoring timescale of ~50 years
bgcparams.kappazML=2e-5; % the diapycnal mixing rate at the bottom of the ML
%bgcparams.kappazML=0;
% the corresponding depth scale is the instanteneous ML depth

bgcparams.kappah=2.5e3; % a horizontal diffusivity; the corresponding length scale is bgcparams.Atllength below
%bgcparams.kappah=0;
% 2.5e3 yields a restoring timescale of 80 years
% for comparison, present day AMOC yields a restoring timescale ~ 10 years

% the maximum fraction of Dwinter in the surface layer
bgcparams.Dmldfrac=.99; 
% active upper-pycnocline layers extend to only a fraction of the depth scale of the AMOC
% upper layer 
% this is a fudge factor that enables us to separate the seasonal cycle of mixing from the 
% AMOC depth, which does not exactly represent the average MLD across the
% entire N. Atlantic  (Dwinter/Dmoc)
bgcparams.Dwinterfrac = 0.17;



%% Sinking flux e-folding depth scale
bgcparams.delta=500; % e.g. Armstrong et al. 2002; Dutkiewicz et al. 2001 DSR have the same value: 400 m 
% this is not that well constrained - and really need to represent losses to
% density classes that don't outcrop, which may be shallower than Dwinter
% elsewhere in the basin?
% abyssal and subtropical nutrient concentration
bgcparams.Na = 18;
bgcparams.Nsub = 19;
% this may decline at a rate mmol/m3 per day
bgcparams.Nsubrate=0;

% basin geometry
bgcparams.Atlwidth=4000E3;
bgcparams.Atllength=2750E3;
end
