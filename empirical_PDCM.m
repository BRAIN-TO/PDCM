clear all; close all;


% Specify P-DCM 
%--------------------------------------------------------------------------
%% Specify model parameters
B0      = 3; % field strength
TE      = 0.04;     % echo time (secs)
TR      = 2.0;  % Repetition time (secs)

%% Times series
% load timeseries files
timeseries = readtable('/home/yuexin/Documents/[External] Re_ Task Based fMRI DCM Data/DCM_CAMH_singlesubject_sub-0006/DCM_CAMH_singlesubject_sub-0006/sub-0006_ROI-timeseries_9ROI.csv'); 
% get the names of ROIs
roi_names = timeseries.Properties.VariableNames;
% get timeseries
timeseries = table2array(timeseries);

% DCM.Y.y should be a ns x nr array.
% ns: length of timeseries. 
% nr: number of ROIs
Y.y     = timeseries;
ns      = size(Y.y,1);
nr      = size(Y.y,2);
M.l     = nr;
R       = nr;

%if ~DCM_info.Y.X0
%    error("DCM.Y.X0 missing.");
%end

%% External Inputs
taskinputs  = readtable('/home/yuexin/Documents/[External] Re_ Task Based fMRI DCM Data/DCM_CAMH_singlesubject_sub-0006/DCM_CAMH_singlesubject_sub-0006/sub-0006_U_events_decoded.csv');
% cnames contain the names of all inputs, driving input at the end
cnames      = {'modulatory','driving'};
% stimulus_ u should be a cell array of size n_events x (2 n_inputs)
% e.g., cell 1 is driving input onset, cell 2 is driving input duration.
% cell 3 is modulatory input, cell 4 is modulatory input duration.
idx = taskinputs.OutcomexVolatility_amp > 0;
task_modulatory = taskinputs(idx,:);
stimulus_u  = {table2array(taskinputs(:,3))' ones(size(taskinputs(:,3)))';table2array(task_modulatory(:,3))' ones(size(task_modulatory(:,3)))'};
% onsets should be a cell array of size 1 x nu
% nu: number of inputs
% Each cell contains a double array of onsets
onset{1}      = [stimulus_u(2,1)',stimulus_u(1,1)];
% duration should be a cell array of size 1 x nu
% nu: number of inputs
% Each cell contains a double array of duration, corresponding to onsets
duration{1}    = [{ones(1,80)*(TR/600)}, {ones(1,140)*(TR/600)}];

dt = TR / 60;


%% create DCM.U and DCM.Y
cutoff = Inf;   % if some low pass filtering has to be done ... otherwise Inf for none

DCM = create_SPM_file_for_DCM(Y,ns,TR,onset,duration,cnames,cutoff,roi_names,TR/dt);

%% Prepare M
M   = PDCM_priors_YX(R,zeros(R,R),zeros(R,R),zeros(R,2));
%% Other Model specification
M.delays = ones(1,nr)*(TR/2);
M.TE    = TE;
M.B0    = B0;
M.m     = nr;
%M.n     = 6;         
%M.N     = 64;
M.dt    = DCM.U.dt;

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


%% Specify Connectivity parameters
A_matrix    = readtable('/home/yuexin/Documents/[External] Re_ Task Based fMRI DCM Data/DCM_CAMH_singlesubject_sub-0006/DCM_CAMH_singlesubject_sub-0006/matrix_A_endogenous_9ROI.csv');
pE.A        = table2array(A_matrix);
B_matrix1    = readtable('/home/yuexin/Documents/[External] Re_ Task Based fMRI DCM Data/DCM_CAMH_singlesubject_sub-0006/DCM_CAMH_singlesubject_sub-0006/matrix_B_model1_9ROI.csv');
pE.B        = table2array(B_matrix);
C_matrix    = readtable('/home/yuexin/Documents/[External] Re_ Task Based fMRI DCM Data/DCM_CAMH_singlesubject_sub-0006/DCM_CAMH_singlesubject_sub-0006/matrix_C_driving_9ROI.csv');
pE.C        = [zeros(R,size(DCM.U.u,2)) table2array(C_matrix)];

pE.C
pE.A        = pE.A.*exp(-2);  % adjust prior as needed
pE.B        = zeros(R,R); % modulatory

pE.D        = zeros(R);    % nonlinear modulation 
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

%% specify which parameters will be estimated (by specifying prior variance)
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

%% Create SPM file for DCM



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

