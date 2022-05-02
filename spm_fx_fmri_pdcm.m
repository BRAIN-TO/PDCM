function [f dfdx] = spm_fx_fmri_pdcm(x,u,P,M)
% state equation for a dynamic [bilinear/nonlinear/Balloon] model of fMRI
% responses
% FORMAT [y] = spm_fx_fmri(x,u,P,M)

% x      - state vector
%   x(:,1) - excitatory neuronal activity             ue
%   x(:,2) - vascular signal                          s
%   x(:,3) - rCBF                                  ln(f)
%   x(:,4) - venous volume                         ln(v)
%   x(:,5) - deoyxHb                               ln(q)
%  [x(:,6) - inhibitory neuronal activity             ui]
%
% y      - dx/dt
%
%___________________________________________________________________________
%
% References for hemodynamic & neuronal state equations:
%
% 1. Havlicek M, Roebroeck A, Friston KJ, Gardumi, A, Ivanov D, Uludag K,
%    Physiologically informed dynamic causal modeling of fMRI data 
%    NeuroImage 122: 355-372, 2015.
% 2. Buxton RB, Wong EC & Frank LR. Dynamics of blood flow and oxygenation
%    changes during brain activation: The Balloon model. MRM 39:855-864,
%    1998.
% 3. Friston KJ, Mechelli A, Turner R, Price CJ. Nonlinear responses in
%    fMRI: the Balloon model, Volterra kernels, and other hemodynamics.
%    Neuroimage 12:466-477, 2000.
% 4. Stephan KE, Kasper L, Harrison LM, Daunizeau J, den Ouden HE,
%    Breakspear M, Friston KJ. Nonlinear dynamic causal models for fMRI.
%    Neuroimage 42:649-662, 2008.
% 5. Marreiros AC, Kiebel SJ, Friston KJ. Dynamic causal modelling for
%    fMRI: a two-state model.
%    Neuroimage. 2008 Jan 1;39(1):269-78.

%
%__________________________________________________________________________

% Martin Havlicek
% $Id: spm_fx_fmri_pdcm.m 2019 $

%----------------------------------------------------------------------
A   = full(P.A);                       % linear parameters
B   = full(P.B);                       % bi-linear parameters
C   = P.C/16;                          % exogenous parameters
D   = full(P.D);                       % nonlinear parameters

n   = size(A,2);
if isempty(M.Tn) && length(P.mu) == 1
   Tn = ones(n,1);
else
   Tn = M.Tn; 
end

if isempty(M.Tv) && length(P.visco_in) == 1
   Tv = ones(n,1);
else
   Tv = M.Tv; 
end

if isempty(M.Tm) && length(P.nratio) == 1
   Tm = ones(n,1);
else
   Tm = M.Tm; 
end

if isempty(M.Tc) && length(P.decay2) == 1
   Tc = ones(n,1);
else
   Tc = M.Tc; 
end


% Local neuronal parameters:
%==========================================================================
%   N(1) - inhibitory-excitatory connection (IE)            mu     (Hz)
%   N(2) - inhibitory gain  (EI and II)                     lambda (Hz)
%--------------------------------------------------------------------------
N     = [0.8 0.2];
sigma = 0.5*exp(P.sigma);
A     = A - diag(diag(A)) - diag(sigma.*exp(diag(A)));

nb = size(B,3);
for i = 1:nb
    A = A + u(i)*B(:,:,i);
end

nd = size(D,3);
for i = 1:nd
    A = A + x(i,1)*D(:,:,i);
end


EE    = A;
mu    = P.mu;
lam   = P.lambda;
if length(mu)<n
   mu = Tn*mu; 
end
if length(lam)<n
   lam = Tn*lam; 
end

if ~isempty(P.Bmu)
    for i = 1:size(P.Bmu,2);
        mu = mu + P.Bmu(:,i)*u(nb+i);
    end
end
if ~isempty(P.Blambda)
    for i = 1:size(P.Blambda,2); % different rows are associated with different inputs
        lam = lam + P.Blambda(:,i)*u(nb+i); % same inputs for mu and lambda
    end
end

IE    = diag(sigma.*N(1).*exp(mu)); % global scaling by sigma
EI    = diag(N(2).*exp(lam));
II    = EI;

% Neurovascular parameters:
%==========================================================================
%   V(1) -     decay1 (Hz)
%   V(2) -     gain   (Hz)
%   V(3) -     decay2 (Hz)

V     = [0.6 1.5 0.6];
%--------------------------------------------------------------------------
de1   = V(1).*ones(n,1);

if length(P.ga)<n
    ga = V(2).*sum(Tc*exp(P.ga),2); 
else
    ga = V(2).*exp(P.ga);
end

if length(P.decay2)<n
    de2 = V(3).*sum(Tc*exp(P.decay2),2); 
else
    de2 = V(3).*exp(P.decay2);
end
%
% Hemodynamic parameters:
%==========================================================================
%   H(1) - mean transit time                              transit (sec)
%   H(2) - Grubb's exponents                              alpha   (-)
%   H(3) - n-ratio  (f-1)/(m-1)                           nr      (-)
%   H(4) - viscoelastic time      (inflation)             visco   (sec)
%   H(5) - viscoelastic time      (deflatiob)             visco   (sec)
%--------------------------------------------------------------------------
H     = [2 0.35 3 3 6];
% transit time
%--------------------------------------------------------------------------
if length(P.transit)<n
    tt = H(1).*sum(Tv*exp(P.transit),2); 
else
    tt = H(1).*exp(P.transit);
end
% alpha
%------------------------------------------------------------------------
if length(P.alpha)<n
    al = H(2).*sum(Tv*exp(P.alpha),2); 
else
    al = H(2).*exp(P.alpha);
end
% n-ration
% %------------------------------------------------------------------------
if length(P.nratio)<n
    nr = H(3).*sum(Tm*exp(P.nratio),2);
else
    nr      = H(3).*exp(P.nratio);
end
% viscoelastic constant for inflation and deflation phase
if length(P.visco_in)<n
    ve_in = H(4).*sum(Tv*exp(P.visco_in),2); 
else
    ve_in   = H(4).*exp(P.visco_in);
end
if length(P.visco_de)<n
   ve_de = H(5).*sum(Tv*exp(P.visco_de),2);
else
   ve_de = H(5).*exp(P.visco_de);
end

ve       = ve_in;
x(:,3:5) = exp(x(:,3:5));

%d        = max([size(B,3),size(P.Bmu,2),size(P.Blambda,2)]);
u        = u(:);


% model equations:
f(:,1)   = EE*x(:,1)  - IE*x(:,6) + C*u;
f(:,6)   = -II*x(:,6) + EI*x(:,1);


% % m = (f-1 + nr)/nr - oxygen metabolism
% %--------------------------------------------------------------------------
m        = (x(:,3)-1+nr)./nr;
% implement differential state equation y = dx/dt (hemodynamic)
%--------------------------------------------------------------------------
f(:,2)   =  x(:,1) - de1.*(x(:,2));
dfin     = (ga.*x(:,2) - de2.*(x(:,3)-1));
f(:,3)   =  dfin./x(:,3);
%
% simple test for inflation and deflation
fv_de       = (tt.*x(:,4).^(1./al) + ve_de.*x(:,3))./(tt+ve_de);

dv_de = (x(:,3) - fv_de)./tt;

ve(dv_de<0)  = ve_de(dv_de<0);
% % Fout = f(v) - outflow
% %--------------------------------------------------------------------------
fv       = (tt.*x(:,4).^(1./al) + ve.*x(:,3))./(tt+ve);
f(:,4)   = (x(:,3) - fv)./(tt.*x(:,4));
f(:,5)   = (m - fv.*x(:,5)./x(:,4))./(tt.*x(:,5));
f        = f(:);

if nargout>1
% state independent part of Jacobian matrix:


dfdx{1,1} = EE;
for i = 1:size(D,3)
    Di  = D(:,:,i) + diag((diag(EE) - 1).*diag(D(:,:,i)));
    dfdx{1,1}(:,i) = dfdx{1,1}(:,i) + Di*x(:,1);
end
dfdx{1,6} = -IE;
dfdx{2,1} = speye(n,n);
dfdx{2,2} = diag(-de1);
dfdx{3,2} = diag(ga./x(:,3));
dfdx{3,3} = diag(-(de2+ga.*x(:,2))./(x(:,3).^2));
dfdx{4,3} = diag(1./(x(:,4).*(tt+ve)));
dfdx{4,4} = diag(-(al.*x(:,3)+(x(:,4).^(1./al)).*(1-al))./(al.*(x(:,4).^2).*(tt + ve)));
dfdx{5,3} = diag((1./nr - (ve.*x(:,5))./(x(:,4).*(tt + ve)))./(tt.*x(:,5)));
dfdx{5,4} = diag(((al.*ve.*x(:,3) - tt.*x(:,4).^(1./al) + al.*tt.*x(:,4).^(1./al)))./(al.*tt.*(x(:,4).^2).*(tt + ve)));
dfdx{5,5} = diag(-(nr + x(:,3) -1)./(nr.*tt.*x(:,5).^2));
dfdx{6,1} = EI;
dfdx{6,6} = -II;

dfdx      = spm_cat(dfdx);

end

