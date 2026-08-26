clear all; close all;


% run P-DCM on the SPM Attention dataset
%--------------------------------------------------------------------------

load('/home/yuexin/Documents/DCM_CAMH_singlesubject/U.mat');  % load DCM_info, will not be used directly as DCM
DCM_info=DCM;
clear DCM;

fprintf("TR = %d \n", DCM_info.Y.dt)     % TR
if ~DCM_info.Y.X0
    error("DCM.Y.X0 missing.");
end

load('/home/yuexin/Documents/DCM_CAMH_singlesubject/sub1_r1_ACC-R.csv')
DCM.Y.y(:,1) = sub1_r1_ACC_R;
load('/home/yuexin/Documents/DCM_CAMH_singlesubject/sub1_r2_AI-L.csv')
DCM.Y.y(:,2) = sub1_r2_AI_L;
load('/home/yuexin/Documents/DCM_CAMH_singlesubject/sub1_r3_AI-R.csv')
DCM.Y.y(:,3) = sub1_r3_AI_R;
load('/home/yuexin/Documents/DCM_CAMH_singlesubject/sub1_r4_dlPFC-R.csv')
DCM.Y.y(:,4) = sub1_r4_dlPFC_R;
load('/home/yuexin/Documents/DCM_CAMH_singlesubject/sub1_r5_mPFC-R.csv')
DCM.Y.y(:,5) = sub1_r5_mPFC_R;
load('/home/yuexin/Documents/DCM_CAMH_singlesubject/sub1_r6_IOG-L.csv')
DCM.Y.y(:,6) = sub1_r6_IOG_L;
load('/home/yuexin/Documents/DCM_CAMH_singlesubject/sub1_r7_IOG-R.csv')
DCM.Y.y(:,7) = sub1_r7_IOG_R;
load('/home/yuexin/Documents/DCM_CAMH_singlesubject/sub1_r8_thal-R.csv')
DCM.Y.y(:,8) = sub1_r8_thal_R;
load('/home/yuexin/Documents/DCM_CAMH_singlesubject/sub1_r9_VTA-R.csv')
DCM.Y.y(:,9) = sub1_r9_VTA_R;


% small rescalling
%scale   = max(max((DCM.Y.y))) - min(min((DCM.Y.y)));
%scale   = 4/max(scale,4);
%DCM.Y.y     = DCM.Y.y*scale;
%DCM.Y.scale = scale;


% specify model parameters (or scanling constants)
B0      = 3; % field strength
TE      = DCM_info.TE;     % echo time (secs)
nr      = size(DCM.Y.y,2);
M.delays = ones(1,nr)*DCM_info.Y.dt/2; 
M.TE    = TE;
M.B0    = B0;
M.m     = nr;
%M.n     = 6;         
M.l     = nr;
%M.N     = 64;
M.dt    = DCM_info.U.dt;
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
pE.A        = DCM_info.a.*exp(-2);  % endogenous
pE.B        = zeros(n,n); % modulatory

pE.D        = zeros(n);    % nonlinear modulation 
pE.C        = [0	0;
                0	0;
                0	0;
                0	0;
                0	0;
                1	0;
                1	0;
                0	0;
                0   0]*exp(0); % encoding of driving inputs
% neuronal parameters (scaling constants)
pE.mu       = zeros(1,1);
pE.lambda   = zeros(1,1);
pE.sigma    = zeros(1);
pE.Bmu      = [];
pE.Blambda  = [];
% NVC parameters (scaling constants)
pE.decay2   = zeros(1,1);
pE.ga       = zeros(1,1);
% Hemodynamic parameters (scaling constants)
pE.transit  = zeros(1,1); 
pE.alpha     = zeros(1,1); 
pE.visco_de  = zeros(1,1); 
pE.visco_in  = zeros(1,1); 
pE.nratio    = zeros(1,1);
pE.V0        = zeros(1,1);

% specify which parameters will be estimated (by specifying prior variance)
spC          = spm_unvec(spm_vec(pE)*0,pE);
spC.C        = [0	0;
                0	0;
                0	0;
                0	0;
                0	0;
                1	0;
                1	0;
                0	0;
                0   0]*exp(0);
% specify connectivity structure 
spC.A        = DCM_info.a*exp(0);
spC.B(:,:,1) = DCM_info.b*exp(0);        
                

spC.mu       = ones(1,1)*exp(-2);
spC.sigma    = ones(1)*exp(-1);
spC.lambda   = ones(1,1)*exp(-2);

spC.transit  = ones(1,1)*exp(-4);

spC.visco_in = ones(1,1)*exp(-1);
spC.visco_de = ones(1,1)*exp(-1);

spC.V0       = ones(1,1)*exp(-4);

pC           = diag(spm_vec(spC));

M.pE         = pE;
M.pC         = pC;
DCM.M        = M;

load('/home/yuexin/Documents/DCM_CAMH_singlesubject/stimulus_u.csv');
DCM.U.dt = DCM_info.U.dt;
DCM.U.name = DCM_info.U.name;
DCM.U.u = zeros(size(DCM_info.U.u));
for i = 1:size(stimulus_u,2)
    for j =1:size(stimulus_u,1)
        if stimulus_u(j,i)>0
            duration=stimulus_u(j,i);
            DCM.U.u(j:j+duration,i) = 1;
        end
    end
end
% Run the model inversion:
[Ep,Cp,Eh,F] = spm_nlsi_GN(M,DCM.U,DCM.Y);
% Ep - estimated parameters (same structure as pE above)


% get the time-courses with estimated paramteres
[y X]        = spm_int_IT(Ep,DCM.M,DCM.U);

% y - BOLD time-courses (time x region)
% X - physiological time-courses (time x (physiological variable per region))

% save variables in DCM_pdcm.mat
DCM.F = F;
DCM.Ep = Ep;
DCM.Cp = Cp;
DCM.y = y;
DCM.X = X;
DCM.v = length(DCM.y);
save("DCM_pdcm_CAMH.mat","DCM","F","Ep","Cp");

