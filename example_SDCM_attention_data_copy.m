clear all; close all;

% This script changes the input order to show it does not affect the result


% run S-DCM on the SPM Attention dataset
%--------------------------------------------------------------------------

load('SPM.mat');  % load SPM data structure file

% Prepare data:
load('VOI_V1_1.mat');     % load data for V1 regions
DCM.Y.y(:,1) = xY.u;
%c
DCM.xY(:,1) = xY;
load('VOI_V5_1.mat');     % load data for V5 regions
DCM.Y.y(:,2) = xY.u;
%c
DCM.xY(:,2) = xY;
load('VOI_SPC_1.mat');    % load data for SPC regions
DCM.Y.y(:,3) = xY.u;
%c
DCM.xY(:,3) = xY;

DCM.Y.dt   = SPM.xY.RT;     % TR
DCM.Y.X0   = xY.X0;  % low freqency fluctuations % SDCM

% small rescalling
scale   = max(max((DCM.Y.y))) - min(min((DCM.Y.y)));
scale   = 4/max(scale,4);
DCM.Y.y     = DCM.Y.y*scale;
DCM.Y.scale = scale;



% input index 
u_idx = [1 2 3]; % SDCM
% Specify inputs for PDCM model
%--------------------------------------------------------------------------
Sess   = SPM.Sess(1);
U.name = {};
U.u    = [];
for i = 1:length(u_idx)
    u = u_idx(i);
    for  j = 1:length(Sess.U(u).name)
        U.u             = [U.u Sess.U(u).u(33:end,j)];
        U.name{end + 1} = Sess.U(u).name{j};
    end
end
U.dt   = Sess.U(1).dt;
% U.u(U.u>1) = 1; % SDCM
DCM.U  = U;

% specify model parameters (or scanling constants)
% B0      = 3; % field strength % SDCM
TE      = 0.04;     % echo time (secs)
nr      = size(DCM.Y.y,2);
M.delays = ones(nr,1)*DCM.Y.dt/2; % SDCM
M.TE    = TE;
% M.B0    = B0; % SDCM
M.m     = nr;
% M.n     = 6;   % does not change results between script and spm     
M.l     = nr;
M.N     = 64; %???
M.dt    = DCM.U.dt; % does not change results between script and spm
M.ns    = size(DCM.Y.y,1);
M.TE    = TE;
M.x     = zeros(M.m,5); % SDCM
M.n     = size(M.x(:),1); % SDCM
M.IS    = 'spm_int_sdcm'; % SDCM

M.f   = @spm_fx_fmri;     % physiological model function SDCM
M.g   = @spm_gx_fmri;     % BOLD model function SDCM
M.Tn  = [];                    %    
M.Tc  = [];
M.Tv  = [];
M.Tm  = [];




n           = nr;
% Connectivity parameters
pE.A        = zeros(n);  % endogenous
pE.B        = zeros(n,n,3); % modulatory % SDCM

pE.D        = zeros(n,n,0);   % nonlinear modulation % SDCM
pE.C        = [[1 0 0]',zeros(n,2)]; % encoding of driving inputs
% neuronal parameters (scaling constants)
% pE.mu       = zeros(n,1); % SDCM
% pE.lambda   = zeros(n,1); % SDCM
% pE.sigma    = zeros(1); % SDCM
% pE.Bmu      = []; % SDCM
% pE.Blambda  = []; % SDCM
% NVC parameters (scaling constants)
pE.decay    = zeros(1,1); % SDCM
% pE.decay2   = zeros(n,1); % SDCM
% pE.ga       = zeros(n,1); % SDCM
% Hemodynamic parameters (scaling constants)
pE.transit  = zeros(n,1); 
% pE.alpha     = zeros(n,1); % SDCM
% pE.visco_de  = zeros(n,1); % SDCM
% pE.visco_in  = zeros(n,1); % SDCM
% pE.nratio    = zeros(n,1); % SDCM
% pE.V0        = zeros(1,1); % SDCM
pE.epsilon   = zeros(1,1); % SDCM

% specify which parameters will be estimated (by specifying prior variance)
spC          = spm_unvec(spm_vec(pE)*0,pE);
spC.C        = [pE.C]*exp(0);
% specify connectivity structure % SDCM
A0           = [1 1 0; % V5 -> V1
                1 1 1; % V1 -> V5 and SPC -> V5
                0 1 1]; % V5-> SPC
spC.A        = A0*exp(0);
spC.B(:,:,2) = [0 0 0;     % related to the first input u(:,1) 
                1 0 0;     % modulataion by motion (V1->V5)
                0 0 0];        
spC.B(:,:,3) = [0 0 0;     % related to the second input u(:,2) 
                0 0 1;     % modulataion by attention (SPC->V5)
                0 0 0];                   

% spC.mu       = ones(n,1)*exp(-2); % SDCM
% spC.sigma    = ones(1)*exp(-1); % SDCM
% spC.lambda   = ones(n,1)*exp(-2); % SDCM
% spC.decay2   = ones(n,1)*exp(-2)*0; % SDCM
spC.decay    = ones(1,1)*exp(-2)*0; % SDCM

spC.transit  = ones(n,1)*exp(-4);
spC.epsilon = ones(1,1)*exp(-4)*0; % SDCM

% spC.visco_in = ones(n,1)*exp(-1); % SDCM
% spC.visco_de = ones(n,1)*exp(-1); % SDCM

% spC.V0       = ones(1,1)*exp(-4); % SDCM

pC           = diag(spm_vec(spC));

M.pE         = pE;
M.pC         = pC;
DCM.M        = M;

% for SDCM 
DCM.v = size(DCM.Y.y,1);
% SDCM
load("sdcm_spm_parameter.mat","M", "Y");
DCM.M.pC = M.pC;
DCM.M.pE = M.pE;
DCM.M.dt = 0.5;
% for SDCM
DCM.M.hE  = M.hE;
DCM.M.hC  = M.hC;
%DCM.Y = Y;

% Run the model inversion:
[Ep,Cp,Eh,F] = spm_nlsi_GN(DCM.M,DCM.U,DCM.Y);
% Ep - estimated parameters (same structure as pE above)


% get the time-courses with estimated paramteres
[y X]        = spm_int_sdcm(Ep,DCM.M,DCM.U);

% y - BOLD time-courses (time x region)
% X - physiological time-courses (time x (physiological variable per region))

% save variables in DCM_pdcm.mat
DCM.F = F;
DCM.Ep = Ep;
DCM.Cp = Cp;
DCM.y = y;
DCM.X = X;
save("DCM_sdcm.mat","DCM","F","Ep","Cp");

