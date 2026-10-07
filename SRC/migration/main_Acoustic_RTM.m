%% ========================================================================
%  Acoustic Synthetic Data Generation and Reverse Time Migration (RTM) Using Modified Nearly Analytical
%  Discritization (MNAD) Numerical Solver
% ========================================================================
%  Purpose:
%    Generates synthetic acoustic seismic data, background data, and
%    appying the RTM on the scattered data
%    Loads a velocity
%    model, builds a PML-padded domain, sets up source/receiver
%    acquisition geometry, builds the MNAD coefficient
%    matrices, and forward-models shot gathers in parallel.
%
%  Compatible with both synthetic acquisitions OBS / streamer geometries
%
%  Outputs:
%    - Illuminations             : output directories icluding illumination 
%    - Wavefields                : output directories icluding nt/nf/2 time
%      slice of wavefield
%    - data_Acoustic.mat         : forward-modeled shot gather data
%    - profile                   : key parameters into a profile struct 
% ========================================================================
%  Author: Toktam Zand
%  Contact: tzand@igf.edu.pl

clear ; clc ; %close all ;

IW = true ; % Create output folders for illumination maps and wavefield
            % snapshots (used by later steps in the RTM workflow, 
            % not by this script directly). If you no need change IW =
            % false ;
if IW == true
    if ~exist('Illuminations','dir'), mkdir Illuminations; end
    if ~exist('Wavefields','dir'), mkdir Wavefields; end
end

% % Set the number inside parpool based on your system
poo = gcp('nocreate');
if isempty(poo), parpool(2);clear poo; end

% Log start time for run-time tracking
date = string(datetime('now')) ;
formatSpec = 'Start -> Date: %s \n';
fprintf(formatSpec,date)

% % Thickness of PML (Perfectly Matched Layer) absorbing boundary, in grid cells
d_pml = 40 ;

% % ______________________________________________________________________% %
% %  ---------------------- Synthetic Data Generation ------------------

%% --- Load velocity model ---
% newname = struct2cell(load('Velocity.mat')) ;
% V = newname{1} ;
V = 1500.*ones(100,200) ;
V(70:end,:) = 2500 ; 
[nz,nx] = size(V) ;            % Number of grid points in x and z directions
V_min = min(V(:)) ;
V_max = max(V(:)) ;

%% --- Frequency content ---
fc = 15 ;                       % Dominant frequency of the Ricker wavelet [Hz]
f_max  = 3.*fc ;               % Approximate maximum frequency content (for e.g. stability/dispersion checks

profile.fc = fc ;    
profile.f_max = f_max ; 

%% --- Domain Discritization ---
% Compute stable time step [s] and spatial steps [m]

dx = V_min/(2.5*f_max) ;
dz = dx ;
% % dx and dz may change regarding your need but must be smaller than what
% calculated above corresponding maximum velocity and frequency 

dt = (0.65*dx)/V_max ;


%% --- Time axis ---
% The recording time
T_max = 3 ;                               % Set the total recording duration in seconds
dt = 0.8.*dt ;                            % Time step [s] 
t = 0:dt:T_max-dt ;                       % Time vector
nt = length(t) ;                          % Number of time samples
fs = 1/dt ;                               % Sampling rate [Hz]

profile.dt = dt ;        
profile.T_max = T_max ;        

%% --- Spatial grid (x-direction) ---
dx = floor(dx) ;                         % Grid spacing in x [m]
x = (0:nx-1)*dx ;

profile.dx = dx ;         
profile.nx = nx ;         

%% --- Spatial grid (z-direction) ---
dz = floor(dz) ;                          % Number of grid points and spacing [m] in z
z = (0:nz-1)*dz ;

profile.dz = dz ;         
profile.nz = nz ;        

%% --- Acquisition geometry (synthetic) ---
zsc = 2*dz ;                                   % Set sources depth [m]
zrc = 4*dz ;                                   % Set receivers depth [m]
h_offset = 1*dx ;                              % Set the horizontal offset of receiver array [m] if needed
ds = 10*dx ;                                    % Source spacing in lateral direction [m]
dr = 1*dx ;                                   % Receiver spacing in lateral direction [m]
slag = 1*dx ;                               % Source lag from the sides of model in lateral direction [m]
rlag = dx ;                                 % Source lag from the sides of model in lateral direction [m]

% Build source/receiver coordinates (in meters) for all shots
[Source, Receiver] = acquisi(V,dx,dz,zsc,zrc,h_offset,ds,dr,slag,rlag) ;
ns = length(Source) ;                                   % The number of source shots

% Preallocate per-shot coordinate cells
xs = cell(1,ns) ;
zs = cell(1,ns) ;
xr = cell(1,ns) ;
zr = cell(1,ns) ;

% Convert source/receiver positions from meters to grid indices,
% and shift by d_pml to account for the PML padding added to the model
for ss = 1:ns
    xs{ss} = d_pml + Source(1,ss)./dx ;        % Source x grid index
    zs{ss} = d_pml + Source(2,ss)./dz ;        % Source z grid index

    xr{ss} = d_pml + Receiver{1,ss}./dx ;      % Receiver x grid indices
    zr{ss} = d_pml + Receiver{2,ss}./dz ;      % Receiver z grid indices
end

profile.ns = ns ;
profile.Source = Source ;
profile.Receiver = Receiver ;

%% --- Source wavelet ---
% Generate Ricker wavelet
t0 = 1.1/fc ;        % [s]   delay wavelet
tau = pi*(t-t0)*fc ;         % [-]   Ricker wavelet argument
w  = ((1 - 2*tau.^2) .* exp( -tau.^2 ))' ; % Analytical Ricker Wavelet solution

profile.w = w ;

%% --- Pad velocity model with PML boundary layers ---
% Extend the model edges outward by d_pml cells on each side (x then z),
% by repeating the edge values — this is a common way to build the PML
% "buffer" region  around the physical domain.
v = [repmat(V(:,1),1,d_pml) , V , repmat(V(:,end),1,d_pml)] ;
v = [repmat(v(1,:),d_pml,1) ; v ; repmat(v(end,:),d_pml,1)] ;

%% --- Build MNAD coefficient matrices ---
% Each *_gen function builds one element of the MNAD kernal tensor A used to
% propagate acoustic wavefield "u" and gradient fields vs x and z axces "v" and "w", respectively. 
% The cofficients provided in *_gen can be find in Appendix A of the
% manuscript

Auu = Auu_gen(dx,dz,dt,v) ;
Auv = Auv_gen(dx,dz,dt,v) ;
Auw = Auw_gen(dx,dz,dt,v) ;
Avu = Avu_gen(dx,dz,dt,v) ;
Avv = Avv_gen(dx,dz,dt,v) ;
Avw = Avw_gen(dx,dz,dt,v) ;
Awu = Awu_gen(dx,dz,dt,v) ;
Awv = Awv_gen(dx,dz,dt,v) ;
Aww = Aww_gen(dx,dz,dt,v) ;
save As Auu Auv Auw Avu Avv Avw Awu Awv Aww

%% --- Forward model shot gathers in parallel ---
% Each worker forward-propagates one shot and returns its recorded data
% Function "fwd_PDE_exact_MNAD" solve the Equation 9 in the manuscript.

data = cell(1,ns);
parfor ss = 1:ns
    data{ss} = fwd_PDE_exact_MNAD(dx,xs{ss},xr{ss},dz,zs{ss},zr{ss},v,dt,nt,w,d_pml,ss,Auu,Auv,Auw,Avu,Avv,Avw,Awu,Awv,Aww) ;
end

% Log completion time
date = string(datetime('now')) ;
formatSpec = 'Data generated -> Date: %s \n';
fprintf(formatSpec,date)


save data_MNAD_Acoustic data profile


%% --- Generating the background velocity model --- 
% Smooth background velocity model is generating with applying the gausiian filter on the original velocity model 
% Pad velocity model with PML boundary layers by extending the model edges outward by d_pml cells on each side (x then z),
% by repeating the edge values — this is a common way to build the PML
% "buffer" region  around the physical domain.
sindx = 3 ;
V_s = v_smooth(V,sindx,sindx,sindx) ;

vs = [repmat(V_s(:,1),1,d_pml) , V_s , repmat(V_s(:,end),1,d_pml)] ;
vs = [repmat(vs(1,:),d_pml,1) ; vs ; repmat(vs(end,:),d_pml,1)] ;

%% --- Build MNAD coefficient matrices for background velocity model ---
% Each *_gen function builds one element of the MNAD kernal tensor A used to
% propagate acoustic wavefield "u" and gradient fields vs x and z axces "v" and "w", respectively. 
% The cofficients provided in *_gen can be find in Appendix A of the
% manuscript

Auu = Auu_gen(dx,dz,dt,vs) ;
Auv = Auv_gen(dx,dz,dt,vs) ;
Auw = Auw_gen(dx,dz,dt,vs) ;
Avu = Avu_gen(dx,dz,dt,vs) ;
Avv = Avv_gen(dx,dz,dt,vs) ;
Avw = Avw_gen(dx,dz,dt,vs) ;
Awu = Awu_gen(dx,dz,dt,vs) ;
Awv = Awv_gen(dx,dz,dt,vs) ;
Aww = Aww_gen(dx,dz,dt,vs) ;
save As Auu Auv Auw Avu Avv Avw Awu Awv Aww

%% --- Forward model shot gathers in parallel ---
% Each worker forward-propagates one shot and returns its recorded data
% Function "fwd_PDE_exact_MNAD" solve the Equation 9 in the manuscript.
% This data is generationg by wave propagation in the background model
% generally lack of reflection data and will deducted from the recorded
% data before migration to provide the scattered data

datas = cell(1,ns);     % Data correspond to the background velocity model
data_sct = cell(1,ns);  % Scattered data which is the input of migration step
parfor ss = 1:ns
    datas{ss} = fwd_PDE_exact_MNAD(dx,xs{ss},xr{ss},dz,zs{ss},zr{ss},vs,dt,nt,w,d_pml,ss,Auu,Auv,Auw,Avu,Avv,Avw,Awu,Awv,Aww,1) ;
    data_sct{ss} = data{ss}{1} - datas{ss}{1} ;
end

% Log completion time
date = string(datetime('now')) ;
formatSpec = 'Scattered data generated -> Date: %s \n';
fprintf(formatSpec,date)


save data_scattered_MNAD_Acoustic data_sct profile

%% Generating the approximation of Hessian for preconditioning
illum_tot = 0 ;
for ss = 1:ns
    load(['./Illuminations/illum_' num2str(ss)],'illum')
    illum_tot = illum_tot + illum ;
    delete(['./Illuminations/illum_' num2str(ss) '.mat'])
end
illum = illum_tot./ns ;
save(['./Illuminations/illum_' num2str(0)],'illum')

%% -------------------- RTM --------------------
save acs nx nz 
Gt = @(d,ss) Born_adj_MNAD(dx,dz,nt,dt,vs,xr{ss},zr{ss},d,d_pml,ss,Auu,Auv,Auw,Avu,Avv,Avw,Awu,Awv,Aww) ;

tic
parfor ss = 1:ns
    grad{ss} = feval(Gt,(data_sct{ss}),ss) ;
end
Elapsed = toc

mhh = zeros(nz,nx,ns) ;
for ss = 1:ns
    mhh(:,:,ss) = reshape(grad{ss},nz,nx) ;
end

load(['./Illuminations/illum_' num2str(0)]) ;
img_RTM = mean(mhh,3)./illum ;


save img_RTM_MNAD img_RTM 

