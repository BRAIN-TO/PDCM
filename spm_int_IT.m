function [y X] = spm_int_IT(P,M,U)
% integrates a nonlinear system dx/dt = f(x,u);
% FORMAT [y] = spm_int(P,M,U)
% P   - model parameters
% M   - model structure
%   M.delays - sampling delays (s); a vector with a delay for each output
%
% U   - input structure or matrix
%
% y   - response y = g(x,u,P)
%__________________________________________________________________________

%
% at v = M.ns is the number of samples [default v = size(U.u,1)]
%
% It will also handle timing delays if specified in M.delays
%
%--------------------------------------------------------------------------
%
 
% Martin Havlicek
% $Id: spm_int_IT.m (modified verison of spm_int by Karl Friston)
 
 
% convert U to U.u if necessary
%--------------------------------------------------------------------------
if ~isstruct(U), u.u = U; U = u; end
try, dt = U.dt; catch, U.dt = 1; end
 
% number of times to sample (v) and number of microtime bins (u)
%--------------------------------------------------------------------------
u       = size(U.u,1);
try,  v = M.ns;  catch, v = u;   end


% get expansion point
%--------------------------------------------------------------------------
x = [spm_vec(M.x)];
 
% add [0] states if not specified
%--------------------------------------------------------------------------
try
    M.f = spm_funcheck(M.f);
catch
    M.f = @(x,u,P,M) sparse(0,1);
    M.x = sparse(0,0);
end

 
% output nonlinearity, if specified
%--------------------------------------------------------------------------
try
    g   = spm_funcheck(M.g);
catch
    g   = @(x,u,P,M) x;
    M.g = g;
end
 
 
% delays
%--------------------------------------------------------------------------
try
    D  = max(round(M.delays/U.dt),1);
catch
    D  = ones(M.l,1)*round(u/v);
end


% Evaluation times (t) and indicator array for inputs (su) and output (sy)
%==========================================================================
 
% assume that the input can change anytime
%--------------------------------------------------------------------------
su    = ones(1,u);
 
% get times that the response is sampled
%--------------------------------------------------------------------------
s     = ceil((0:v - 1)*u/v);
for j = 1:M.l
    i       = s + D(j);
    sy(j,:) = sparse(1,i,1:v,1,u);
end
 
% time in seconds
%--------------------------------------------------------------------------
t     = find(su | any(sy));
su    = full(su(:,t));
sy    = full(sy(:,t));
dt    = [diff(t) 0]*U.dt;
 
 
% Integrate
%--------------------------------------------------------------------------
y     = zeros(M.l,v);
X     = zeros(length(x),v);
U.u   = full(U.u);

for i = 1:length(t)
 
      q          = spm_unvec(x,M.x);
     [fx, dfdx]   = M.f(q,U.u(t(i),:),P,M);

    % output sampled
    %----------------------------------------------------------------------
    if any(sy(:,i))
        j      = find(sy(:,i));
        s      = sy(j(1),i);
        q      = spm_vec(M.g(q,[],P,M));
        y(j,s) = q(j);
    end
 
    % compute updated states (with Ito-Taylor approximation)
    %----------------------------------------------------------------------

    x         = x + dt(i)*fx(:) + 0.5*dt(i)^2*(full(dfdx)*fx(:));

    
    % check for convergence
    %----------------------------------------------------------------------
    if norm(x,1) > 1e6, break, end
    X(:,i) = x(:);
end

y      = real(y');
X      = real(X');
