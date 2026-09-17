clear all; close all;


%% Specify P-DCM, Users need to fill out this section.
%--------------------------------------------------------------------------
% Specify model parameters
B0          = 3;        % field strength
TE          = 0.04;     % echo time (secs)
TR          = 2.0;      % Repetition time (secs)
dt          = TR / 60;  % sampling rate for DCM
n_events    = 2;        % number of regressors (events) in total
R           = 9;        % number of ROIs

% Specify matrices
%--------------------------------------------------------------------------
% pE.A should be a R x R matrix where the entry specifies intrinsic
% connections.
A_table    = readtable('/home/yuexin/Documents/[External] Re_ Task Based fMRI DCM Data/DCM_CAMH_singlesubject_sub-0006/DCM_CAMH_singlesubject_sub-0006/matrix_A_endogenous_9ROI.csv','ReadRowNames',true);
A_matrix        = table2array(A_table);
% pE.B should be a R x R matrix where the entry specifies whether the 
% intrinsic connection (corresponding to A) experience connectivity changes.
B_table    = readtable('/home/yuexin/Documents/[External] Re_ Task Based fMRI DCM Data/DCM_CAMH_singlesubject_sub-0006/DCM_CAMH_singlesubject_sub-0006/matrix_B_model2_9ROI.csv','ReadRowNames',true);
B_matrix        = table2array(B_table);
% pE.C should be a R x n_events matrix, the entries in the last column
% represent which regions receives driving inputs
C_table    = readtable('/home/yuexin/Documents/[External] Re_ Task Based fMRI DCM Data/DCM_CAMH_singlesubject_sub-0006/DCM_CAMH_singlesubject_sub-0006/matrix_C_driving_9ROI.csv','ReadRowNames',true);
C_matrix        = [zeros(R,n_events-1) table2array(C_table)];
% pE.Bmu should be a 1 x n_events array, if Bmu or Blambda has
% event-related changes.
pE.Bmu      = [];
pE.Blambda  = [];

% External Inputs
%--------------------------------------------------------------------------
taskinputs  = readtable('/home/yuexin/Documents/[External] Re_ Task Based fMRI DCM Data/DCM_CAMH_singlesubject_sub-0006/DCM_CAMH_singlesubject_sub-0006/sub-0006_U_events_decoded.csv');
% cnames contain the names of all inputs, driving input at the end
cnames      = {'modulatory','driving'};
% stimulus_ u should be a cell array of size n_events x (2 n_inputs)
% e.g., cell 1 is driving input onset, cell 2 is driving input duration.
% cell 3 is modulatory input, cell 4 is modulatory input duration.
idx = taskinputs.OutcomexVolatility_amp > 0;
task_modulatory = taskinputs(idx,:);
stimulus_u  = {table2array(taskinputs(:,3))' ones(size(taskinputs(:,3)))';table2array(task_modulatory(:,3))' ones(size(task_modulatory(:,3)))'};
% onsets should be a cell array of size 1 x n_events
% Each cell contains a double array of onsets
onset{1}      = [stimulus_u(2,1),stimulus_u(1,1)];
% duration should be a cell array of size 1 x n_events
% Each cell contains a double array of duration, corresponding to onsets
duration{1}    = [{zeros(1,80)}, {zeros(1,140)}];

% load timeseries files
%--------------------------------------------------------------------------
timeseries = readtable('/home/yuexin/Documents/[External] Re_ Task Based fMRI DCM Data/DCM_CAMH_singlesubject_sub-0006/DCM_CAMH_singlesubject_sub-0006/sub-0006_ROI-timeseries_9ROI.csv'); 
% get the names of ROIs
roi_names = timeseries.Properties.VariableNames;
% time series should be a n_TR x R matrix
timeseries = table2array(timeseries);


%% Times series
% DCM.Y.y should be a n_TR x R matrix.
% ns: length of timeseries. 
% nr: number of ROIs
Y.y     = timeseries;
ns      = size(Y.y,1); % number of measurements (how many TRs in this run)
nr      = size(Y.y,2); % number of ROIs


%% create DCM.U and DCM.Y
cutoff = Inf;   % if some low pass filtering has to be done ... otherwise Inf for none

DCM = create_SPM_file_for_DCM(Y,ns,TR,onset,duration,cnames,cutoff,roi_names,TR/dt);


%% Prepare DCM.M
M   = PDCM_priors_YX(R,zeros(R,R),zeros(R,R),zeros(R,2));
% Other Model specification
M.delays = ones(1,nr)*(TR/2);
M.TE    = TE;
M.B0    = B0;
M.m     = nr;
M.dt    = DCM.U.dt;
M.x     = zeros(nr,6);
M.IS    = 'spm_int_IT';
M.l     = nr;
M.ns    = ns;
M.nr    = nr;
M.f   = @spm_fx_fmri_pdcm;     % physiological model function
M.g   = @spm_gx_fmri_pdcm;     % BOLD model function
M.Tn  = [];                    
M.Tc  = [];
M.Tv  = [];
M.Tm  = [];


%% Specify Connectivity parameters priors
pE.A        = A_matrix.*exp(-2); % adjust intrinsic connectivity prior as needed
pE.B        = B_matrix.*0;       % modulatory, prior = 0, adjust as needed
pE.D        = zeros(R);         % nonlinear modulation. not used for now 
pE.C        = C_matrix.*exp(-2);    % adjust driving inputs as needed
% neuronal parameters (scaling constants)
pE.mu       = zeros(1,1);
pE.lambda   = zeros(1,1);
pE.sigma    = zeros(1);
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
% specify connectivity structure 
spC.A        = A_matrix.*exp(0); % adjust variance as needed
spC.B        = B_matrix.*exp(0); % adjust variance as needed     
spC.C        = C_matrix*exp(0);                
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


%% Run the model inversion:
[Ep,Cp,Eh,F] = spm_nlsi_GN(DCM.M,DCM.U,DCM.Y);
% Ep - estimated parameters (same structure as pE above)


% get the time-courses with estimated paramteres
[yp X]        = spm_int_IT(Ep,DCM.M,DCM.U);

% y - BOLD time-courses (time x region)
% X - physiological time-courses (time x (physiological variable per region))

% save variables
DCM.Ep = Ep; % posterior means and covariances of estimated parameters
DCM.F = F;      % evidence
DCM.Cp = Cp; 
DCM.Eh = Eh;
DCM.Yp  = yp; % predicted responses
DCM.Xp  = X;
DCM.X = X;
DCM.v = length(DCM.y);

% Fields below are added just to fit in spm dcm QC spm_dcm_fmri_check
DCM.options.two_state = 0;
DCM.v = length(DCM.Y.y);
DCM.y = yp;
DCM.R = DCM.Y.y - DCM.y;

save("DCM_pdcm_CAMH_model2.mat","DCM","F","Ep","Cp");

