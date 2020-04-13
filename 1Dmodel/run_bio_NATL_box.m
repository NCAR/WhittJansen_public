function [Soutmat, toutmat] = run_bio_NATL_box(bgcparams,Sinit,yrs,noutperyr,dovis)
%% Main model timestepping file
% execute this script to run the model. It will visualize as running
%


%% Step 1: start timer
%clear all
tic
if dovis==1, figure('Position',[10 10 1200 700]); end

%% Step 2: select run duration and number of output points, set up timestepping
% total simulation duration
total_sim_time = 86400*360*yrs; % seconds
% number of ouput points
nout = yrs*noutperyr;
nyrsvis=10;% visualize every 10 years, while running
nvismod=noutperyr*nyrsvis; 
% time step
dt = 86400./4;
% number of steps in the simulation
nt = ceil(total_sim_time/dt);
dtout = total_sim_time/nout;
nstepsout = floor(dtout/dt);

%% Step 3: set up grid

% number of grid points
N = 1;
x=1;
dx=1;


%% Step 4: initialize output variables
ntracers = 2;

% output array
Soutmat = zeros(N,ntracers,nout);
toutmat = zeros(1,nout);

%% Step 5: specify the initial condition
% model state matrix

%S = zeros(N,ntracers);
%S(:,1) = 13; % Ns at t=0
%S(:,2) = 13; % Nd at t=0
S=Sinit;
%S=zeros(size(Sequib));
%% Step 6: prep for timestepping
S0 = zeros(N*ntracers,1);
toutmat(1)=0;
Soutmat(:,:,1)=S(:,:);
RHSm1 = zeros(N*ntracers,1); %Adams bashforth rhs
% stack the state vector
for ij = 1:ntracers
    S0(N*(ij-1)+(1:N))=S(:,ij);
end

%% Step 7: specify the biogeochemical model parameters

%% Step 8 specify the physical variables for offline bgc


%% Step 9: run the timestepping loop
display('done initialization, starting timestepping...')

ct = 0;
nt
for it = 1:nt
    
    
    %% reaction vector
    REAC = construct_reac_NATL_box(S0,x,it.*dt,bgcparams);
    
    % construct RHS
    if it > 1
        RHSA = (S0+3.*dt./2.*(REAC)-dt./2.*RHSm1); % adams bashforth
    else
        RHSA = S0+dt.*(REAC); % forward euler
    end
    RHSm1 = REAC;
    
    S0 = RHSA;
    
    if mod(it,nstepsout)==0
        ct=ct+1;
        % save
        for ij = 1:ntracers
            Soutmat(:,ij,ct)=S0(N*(ij-1)+(1:N));
        end
        toutmat(ct)=dt*it;
        
        % visualize
        if (dovis==1 && mod(ct,nvismod) == 0)
            visualize_running_bgc_script_NATL_box_v2;
            disp(it*dt/86400/360)
            pause(.001)
        end
        % report flux
        %pause(.001)
    end
    % reset small values to finite small numbers (don't let N crash)
    mask1  = S0<1E-6;
    S0(mask1) = 1e-6;
    
end
toc
end