function ci = racah_ci(n,l,alpha,beta,gamma,states)
% Calculates the configuration interaction operator matrix, after Racah IV.

if isvector(alpha) && length(alpha)==3
  beta = alpha(2);
  gamma = alpha(3);
  alpha = alpha(1);
  if exist('beta') && iscell(beta)
    states = beta;
  end
elseif isscalar(alpha) & exist('beta') && isscalar(beta) && exist('gamma') && isscalar(gamma)
else
  error('alpha must be a scalar (and beta and gamma must be scalars) or length-3 vector');
end  

if ~exist('states')
  states = racah_states(n,l);
end

ci = zeros(length(states));
for i = 1:length(states)
  S = states{i}{1}; L = racah_lconv(states{i}{2}); 
  v = states{i}{3}; U = states{i}{4}; W = racah_vtow(v,S,l);
  ci(i,i) = alpha * L*(L+1) + beta  * racah_g(U) + gamma * racah_g(W);
end

%crosswhite_f5ci
%ci = alpha.*alphamat + beta.*betamat + gamma.*gammamat;
