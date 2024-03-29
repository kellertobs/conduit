clear; close all;

% set run parameters
runID    =  'testheatbidir';           % run identifier
restart  =  0;                   % restart from file (0: new run; <1: restart from last; >1: restart from specified frame)
nop      =  50;                  % output frame plotted/saved every 'nop' time steps
plot_op  =  1;                   % switch on (1) to live plot results
save_op  =  1;                   % switch on (1) to save output to file
plot_cv  =  1;                   % switch on (1) to live plot iterative convergence

% set model domain parameters
L        =  40;                  % conduit length
R        =  4;                   % conduit radius
N        =  400 + 2;             % number of grid points in z-direction (incl. 2 ghosts)
h        =  L/(N-2);             % grid spacing
 
% set model timing parameters
M        =  1e5;                 % number of time steps to take
tend     =  3600*24*50;          % end time for simulation [s]
dt       =  1e-2;                % initial time step [s]

% set model vesicularity parameters
SlugNo   =  1;                   % slug control number (>1 for shorter more concentrated bubbly slugs)
f1       =  0.00;                % amplitude of random noise
f2       =  0.01;                % amplitude of bubbly core or slug
kf       =  1e-7;                % volume diffusivity [m^2/s]
smth     =  (N/100)^2;           % regularisation of initial bubble distribution

% set model rheology parameters
eta0     =  1e3;                 % background melt viscosity [Pas]
A        = -2.0;                 % bubble weakening exponent
B        =  3.0;                 % crystal stiffening exponent

% set model buoyancy parameters
rhom     =  2400;                % melt phase density [kg/m3]
rhoc     =  2700;                % crystal phase density [kg/m3] 
rhof     =  200;                 % bubble phase density [kg/m3]
g0       =  9.81;                % gravity [m/s2]

% set model temperature/crystallinity parameters
Tw       =  300;
T0       =  1100;                % initial magma temperature [degC]
Tsol     =  800;                 % solidus temperature [degC]
Tliq     =  1200;                % liquidus temperature [degC]
kTm      =  4;                   % thermal diffusivity [m2/s]
kTc      =  1;                   % thermal diffusivity [m2/s]
kTf      =  0.02;                % thermal diffusivity [m2/s]
Cm       =  1400;                % heat capacity [J/kg/K]
Cc       =  1000;                % heat capacity [J/kg/K]
Cf       =  2000;                % heat capacity [J/kg/K]
LH       =  400e3;               % latent heat [J/kg]
tau_c    =  (R/4)^2./kT*rhom*C;  % conduit wall cooling time [s]

% set numerical model parameters
nup      =  100;                  % nonlinear coefficients, residual norms updated every 'nup' iterations
CFL      =  1.0;                 % (physical) time stepping courant number (multiplies stable step) [0,1]
ADVN     =  'FRM';               % advection scheme ('UPW2', 'UPW3', or 'FRM')
theta    =  0.50;                % time-stepping scheme selector (1=BE, 1/2=CN, 0=FE)
rtol     =  1e-4;                % outer its relative tolerance
atol     =  1e-6;                % outer its absolute tolerance
maxit    =  5e3;                 % maximum outer its
alpha    =  0.95;                % inner its step size (multiple of stable step) [0,1]
beta     =  0.20;                % iterative damping parameter [0,1]
zeta     =  2/3;
delta    =  1.0;                 % regularisation of viscosity
etactr   =  1e3;                 % minimum viscosity for regularisation

% create output directory
if ~isfolder(['../out/',runID])
    mkdir(['../out/',runID]);
end

% save input parameters and runtime options (unless restarting)
if restart == 0 
    parfile = ['../out/',runID,'/',runID,'_par'];
    save(parfile);
end

% run code
run('../src/main')
