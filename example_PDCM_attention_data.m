clear all; close all;


% run P-DCM on the SPM Attention dataset
%--------------------------------------------------------------------------

load('SPM.mat');  % load SPM data structure file

% Prepare data:
load('VOI_V1_1.mat');     % load data for V1 regions
DCM.Y.y(:,1) = xY.u;
load('VOI_V5_1.mat');
DCM.Y.y(:,2) = xY.u;
load('VOI_SPC_1.mat');
DCM.Y.y(:,3) = xY.u;

DCM.Y.dt   = SPM.xY.RT;     % TR
DCM.Y.X0   = [ones(size(xY.X0,1),1),xY.X0(:,2:6)];  % low freqency fluctuations

% small rescalling
scale   = max(max((DCM.Y.y))) - min(min((DCM.Y.y)));
scale   = 4/max(scale,4);
DCM.Y.y     = DCM.Y.y*scale;
DCM.Y.scale = scale;



% input index 
u_idx = [2 3 1]; % first two (motiona and attention) will be modulatory 
                 % last one (the third) will be driving
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
U.u(U.u>1) = 1;
DCM.U  = U;


% specify model parameters (or scanling constants)
B0      = 3;
TE      = 0.04;     % echo time (secs)
nr      = size(DCM.Y.y,2);
M.delays = ones(1,nr)*DCM.Y.dt/2; 
M.TE    = TE;
M.B0    = B0;
M.m     = nr;
M.n     = 6;         
M.l     = nr;
M.N     = 64;
M.dt    = DCM.U.dt;
M.ns    = size(DCM.Y.y,1);
M.TE    = TE;
M.B0    = B0;
M.x     = zeros(M.m,6); 
M.IS    = 'spm_int_IT';

M.f   = @spm_fx_fmri_pdcm;     % physiological model function
M.g   = @spm_gx_fmri_pdcm;     % BOLD model function
M.Tn  = [];                    %    
M.Tc  = [];
M.Tv  = [];
M.Tm  = [];

n           = nr;
% Connectivity parameters
pE.A        = zeros(n);  % endogenous
pE.B        = zeros(n,n,2); % modulatory
pE.D        = zeros(n);    % nonlinear modulation 
pE.C        = [zeros(n,2),[1 0 0]']; % encoding of driving inputs
% neuronal parameters (scaling constants)
pE.mu       = zeros(n,1);
pE.lambda   = zeros(n,1);
pE.sigma    = zeros(1);
pE.Bmu      = [];
pE.Blambda  = [];
% NVC parameters (scaling constants)
pE.decay2   = zeros(n,1);
pE.ga       = zeros(n,1);
% Hemodynamic parameters (scaling constants)
pE.transit  = zeros(n,1); 
pE.alpha     = zeros(n,1); 
pE.visco_de  = zeros(n,1); 
pE.visco_in  = zeros(n,1); 
pE.nratio    = zeros(n,1);
pE.V0        = zeros(1,1);

% specify which parameters will be estimated (by specifying prior variance)
spC          = spm_unvec(spm_vec(pE)*0,pE);
spC.C        = [pE.C]*exp(0);
% specify connectivity structure 
A0           = [0 1 0; % V5 -> V1
                1 0 1; % V1 -> V5 and SPC -> V5
                0 1 0]; % V5-> SPC
spC.A        = A0*exp(0);
spC.B(:,:,1) = [0 0 0;     % related to the first input u(:,1) 
                1 0 0;     % modulataion by motion (V1->V5)
                0 0 0];    
    
spC.B(:,:,2) = [0 0 0;     % related to the second input u(:,2) 
                0 0 1;     % modulataion by attention (SPC->V5)
                0 0 0];                   
            
spC.mu       = ones(n,1)*exp(-2);
spC.sigma    = ones(1)*exp(-1);
spC.lambda   = ones(n,1)*exp(-2);
spC.decay2   = ones(n,1)*exp(-2)*0;

spC.transit  = ones(n,1)*exp(-4);

spC.visco_in = ones(n,1)*exp(-1);
spC.visco_de = ones(n,1)*exp(-1);

spC.V0       = ones(1,1)*exp(-4);

pC           = diag(spm_vec(spC));

M.pE         = pE;
M.pC         = pC;
DCM.M        = M;


% Run the model inversion:
[Ep,Cp,Eh,F] = spm_nlsi_GN(M,DCM.U,DCM.Y);
% Ep - estimated parameters (same structure as pE above)


% get the time-courses with estimated paramteres
[y X]        = spm_int_IT(Ep,DCM.M,DCM.U);

% y - BOLD time-courses (time x region)
% X - physiological time-courses (time x (physiological variable per region))






