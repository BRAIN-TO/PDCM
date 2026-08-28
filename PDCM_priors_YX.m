function [M] = PDCM_priors_YX(R,A,B,C)
% PDCM_priors defines structure of parameters for BOLD response. It applies
%                 parameter values as describe in Table 1 of Havlicek et
%                 al., 2015 Table 1A
%
% INPUT:    
%           R - Number of regions
%           A - endogenous fixed connectivity zeros(N,N)
%           B - modulatory input (changes) in connectivity zeros(N,N)
%           C - direct input
%
% OUTPUT:   
%           M - Model including P0 (structure with all default parameters)
%
% AUTHOR:       Yuexin Xi, 28 August, 2026
%
% REFERENCE: Havlicek M, Roebroeck A, Friston K, Gardumi A, Ivanov D, 
%               Uludag K. Physiologically informed dynamic causal modeling 
%               of fMRI data. Neuroimage. 2015 Nov 15;122:355-72. 
%               doi: 10.1016/j.neuroimage.2015.07.078.
%
% EXAMPLE: 
%--------------------------------------------------------------------------
M.R  = R;    % Number of regions

% Neuronal parameter:
%--------------------------------------------------------------------------
P0.sigma = 0.5;   % excitotary self-connection
P0.mu    = 0.8; % inhibitory-excitatory connection following the old script
P0.lam   = 0.2; % inhibitory gain factor
P0.A     = A;   % intrinsic connections between neuron populations
P0.B     = B;   % modulatory effects on A
P0.C     = C;   % 
P0.Bmu   = 0;   % modulatory effects on mu
P0.Blam  = 0;   % modulatory effects on lambda

% NVC parameters:
%--------------------------------------------------------------------------
P0.c1      = 0.6;   % decay of vasoactive signal
P0.c2      = 1.5;   % gain of vasoactive signal
P0.c3      = 0.6;   % decay of blood inflow signal

% LAMINAR HEMODYNAMIC MODEL:
%--------------------------------------------------------------------------
% Baseline physiological parameters:

P0.V0t   = 2.5;   % Total (regional) amount of CBV0 in the gray matter (in mL) [1-6]
P0.w_v = 0.5;   % CBV0 fraction of microvasculature (i.e. venules here )with respect to the total amount 

P0.t0 = 2;     % Mean transit time through microvasculature(in second)
P0.E0   = 0.4;   % Baseline oxygen extraction fraction

% Parameters describing relative relationship between physiological variable:
% CBF-CBV coupling (steady-state)
P0.alpha = 0.35;  

% CBF-CMRO2 coupling (steady-state)
P0.nr = 3;         % n-ratio   (Ref. Buxton et al. (2004) NeuroImage) for CBF-CMRO2 coupling

% CBF-CBV dynamic uncoupling 
P0.tau_in = 3; %  - inflation 
P0.tau_de = 6; %  - deflation

% BOLD SIGNAL MODEL:
%--------------------------------------------------------------------------
P0.Hct = 0.38; % Hematocrit fraction
P0.gyro   = 2*pi*42.6*10^6;    % Gyromagnetic constant for Hydrogen
P0.suscep = 0.264*10^-6;       % Susceptibility difference between fully oxygenated and deoxygenated blood

% Water proton density:
P0.rho_t  = 0.89;                   % For gray matter tissue 
P0.rho_b  = 0.95 - P0.Hct*0.22;   % For blood (venules) Ref. Lu et al. (2002) NeuroImage

% Relaxation rates for 7 T (in sec-1)
P0.R2s_t  = 34; % For gray matter tissue
P0.R2s_b  = 85; % For blood (venules)

% Slope of change in R2* of blood with change in extraction fration during activation 
P0.r0v    = 228;    % For 7T
P0.M0     = 100;    % What is M0?

M.P0 = P0;
M.x  = zeros(R,6);
%M.xn = zeros(R,4);          
%M.xk = zeros(R,4);