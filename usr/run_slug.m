clear;

% set run parameters
runID    =  'slug_upw3';         % run identifier
restart  =  0;                   % restart from file (0: new run; <1: restart from last; >1: restart from specified frame)
nop      =  25;                  % output frame plotted/saved every 'nop' time steps
plot_op  =  1;                   % switch on (1) to live plot results
save_op  =  1;                   % switch on (1) to save output to file
plot_cv  =  1;                   % switch on (1) to live plot iterative convergence

% set model domain parameters
L        =  60;                  % conduit length
R        =  4;                   % conduit radius
N        =  600 + 2;             % number of grid points in z-direction (incl. 2 ghosts)
h        =  L/(N-2);             % grid spacing
 
% set model timing parameters
M        =  1e5;                 % number of time steps to take
tend     =  3600*24;             % end time for simulation [s]
dt       =  1e-2;                % initial time step [s]

% set model vesicularity parameters
SlugNo   =  3;                   % slug control number (>1 for shorter more concentrated bubbly slugs)
f1       =  0.02;                % amplitude of random noise
f2       =  0.20;                % amplitude of bubbly core or slug
kf       =  1e-6;                % volume diffusivity [m^2/s]
smth     =  1/4/h^2;             % regularisation of initial bubble distribution

% set model rheology parameters
eta0     =  1e4;                 % background melt viscosity [Pas]
A        = -1.0;                 % bubble weakening exponent
B        =  2.0;                 % crystal stiffening exponent

% set model buoyancy parameters
rhom     =  2400;                % melt phase density [kg/m3]
rhoc     =  2700;                % crystal phase density [kg/m3] 
rhof     =  200;                 % bubble phase density [kg/m3]
g0       =  9.81;                % gravity [m/s2]

% set model temperature/crystallinity parameters
T0       =  1100;                % initial magma temperature [degC]
Tsol     =  800;                 % solidus temperature [degC]
Tliq     =  1200;                % liquidus temperature [degC]
kT       =  8;                   % thermal diffusivity [m2/s]
C        =  1200;                % heat capacity [J/kg/K]
LH       =  400e3;               % latent heat [J/kg]
tau_c    =  (R/4)^2./kT*rhom*C;  % conduit wall cooling time [s]

% set numerical model parameters
nup      =  50;                  % nonlinear coefficients, residual norms updated every 'nup' iterations
CFL      =  1.0;                 % (physical) time stepping courant number (multiplies stable step) [0,1]
ADVN     =  'UPW3';              % advection scheme ('UPW2', 'UPW3', or 'FRM')
theta    =  0.50;                % time-stepping scheme selector (1=BE, 1/2=CN, 0=FE)
rtol     =  1e-5;                % outer its relative tolerance
atol     =  1e-8;                % outer its absolute tolerance
maxit    =  1e4;                 % maximum outer its
alpha    =  0.90;                % inner its step size (multiple of stable step) [0,1]
beta     =  0.30;                % iterative damping parameter [0,1]
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
