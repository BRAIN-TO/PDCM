function [g,dgdx] = spm_gx_fmri_pdcm(x,u,P,M)
% Simulated BOLD response to input
% FORMAT [g,dgdx] = spm_gx_fmri_pdcm(x,u,P,M)
% g          - BOLD response (%)
% x          - state vector     (see spm_fx_fmri)
% P          - Parameter vector (see spm_fx_fmri)
% M          - model specification structure (see spm_nlsi)
%__________________________________________________________________________
%
% This function implements the BOLD signal model described in: 
% 
% Havlicek M., Roebroeck, A., Friston, K., Gardumi, A., Ivanov, D., Uludag, 
% K. Physiologically informed dynamic causal modeling of fMRI data 
% (2015) NeuroImage 122: 355-372.
%
% as updated verison of 
%
% Stephan KE, Weiskopf N, Drysdale PM, Robinson PA, Friston KJ (2007)
% Comparing hemodynamic models with DCM. NeuroImage 38: 387-401.
%__________________________________________________________________________
% Martin Havlicek 
 
% Biophysical constants 
%==========================================================================
 
% time to echo (TE) 
%--------------------------------------------------------------------------
n   = M.m;

TE  = M.TE;
 
% resting venous volume (%)
%--------------------------------------------------------------------------
if length(P.V0)<n
    V0  = 0.03*ones(n,1)*exp(P.V0);
else
    V0  = 0.03*exp(P.V0);
end
% slope r0 of intravascular relaxation rate R_iv as a function of oxygen 
% saturation S:  R_iv = r0*[(1 - S)-(1 - S0)] (Hz)
%--------------------------------------------------------------------------

 
% resting oxygen extraction fraction
%--------------------------------------------------------------------------
E0  = 0.4;


Hct  = 0.38;       % Hematocrit fraction

B0     = M.B0;              % Field strenght        
gyro   = 2*pi*42.6*10^6;    % Gyromagnetic constant 
suscep = 0.264*10^-6;       % Susceptibility difference

nu0   = suscep*gyro*Hct*B0;


% Water proton density 
rho_t  = 0.89;  % In GM tissue
rho_b  = 0.95 - Hct*0.22;  % In blood  Ref. Lu et al. (2002) NeuroImage


% Relaxation rates (in sec-1):

if B0 == 7
    R2s_t  = 34;         % For tissue (for 7T)
    R2s_b  = 85;           % For venous blood (for 7T)
    r0     = 228;          % Slope of change in R2* of blood with change in extraction fration during activation 
elseif B0 == 3
    R2s_t  = 25;         
    R2s_b  = 40;         
    r0     = 108; 
 
elseif B0 == 1.5
    R2s_t  = 16;         
    R2s_b  = 12;         
    r0     = 15; 
else
    R2s_t = interp1([3 7]',[25 34],B0,'linear','extrap');
    R2s_b = interp1([3 7]',[40 85],B0,'linear','extrap');
    r0    = interp1([3 7]',[108 228],B0,'linear','extrap');
end

% (Baseline) Intra-to-extra-vascular signal ratio
ep   = rho_b./rho_t.*exp(-TE*R2s_b)./exp(-TE*R2s_t);        % For venules


%-Coefficients in BOLD signal model
%==========================================================================

k1     = 4.3.*nu0.*E0.*TE;
k2     = ep.*r0.*E0.*TE;
k3     = 1 - ep; 


%-Output equation of BOLD signal model
%==========================================================================
f   = exp(x(:,3));
v   = exp(x(:,4));
q   = exp(x(:,5));
g  = V0.*(k1.*(1 - q) + k2.*(1 - q./v) + k3.*(1-v))*100;

