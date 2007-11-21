function susceptibility = ic_susc(Hic,mu,T)

if iscell(Hic)
  V = Hic{1}; E = Hic{2};
  if ~isvector(E)
    E = diag(E);
  end
  if (length(E) ~= size(V,1))
    error('V and E must have same size or length');
  end
  E = real(E) - min(real(E));
else
  [V,E] = eig(Hic);
  E = diag(E);
  E = E - min(E);
end

% Physical constants. Taken from NIST Reference on Constants, Units, and 
% Uncertainty, http://physics.nist.gov/cuu/Constants/
mu_B  = 927.400949e-26;    % J/Tesla - Bohr magneton
%mu_B  = 0.46686437;        % cm^{-1} / Tesla - Bohr magneton
k_B   = 1.3806505e-23;     % J/K - Boltzmann constant
Q_e   = 1.60217653e-19;    % C - Charge of electron
N_A   = 6.0221415e23;      % Avogadro's number
h     = 6.62606896e-34;    % Js - Planck's constant
c     = 299792458;         % m/s - speed of light in vacuum

% Defining beta = 1/kT here saves computation later on.
beta = 1 ./ (k_B*T);

% Converts energy levels from meV to J
%E = E .* (Q_e/1000);
E = E .* (h*c*100);

% Defines smallness criteria
small = 1e-5 * (Q_e/1000);

% Defines the limit of energy levels to calculate for
upElimit = 10000 * (h*c*100); % 15000cm-1

iVs = find(E<(upElimit));     % Limits calculations to lower levels (<1eV) to save time
for id_i = 1:length(iVs)
  tic
  ind_i = iVs(id_i);
  % Calculates the degenerate terms
  i_degen = find(abs(E-E(ind_i))<small);
  me_degen = zeros(1,length(i_degen));
  for id_j = 1:length(i_degen)
    ind_j = i_degen(id_j);
    me_degen(ind_j) = V(:,ind_i)' * mu * V(:,ind_j);
  end
  degen(ind_i) = sum(-me_degen.^2);

  % Calculates the nondegenerate terms
  i_ndegen = find( (abs(E-E(ind_i))>=small) & E<(upElimit) );  % Limits calcs to levels <1eV to save time 
  me_ndegen = zeros(1,length(i_ndegen)); ndegen(:,ind_i) = zeros(length(T),1);
  for id_j = 1:length(i_ndegen)
    ind_j = i_ndegen(id_j);
    me_ndegen = V(:,ind_i)' * mu * V(:,ind_j);
    temp_fact = ( exp(-beta.*E(ind_i)) - exp(-beta.*E(ind_j)) ) ./ ( E(ind_j) - E(ind_i) );
    % Sum sigma'' over gamma'
    ndegen(:,ind_i) = ndegen(:,ind_i) + ( (me_ndegen^2) .* temp_fact' );
  end

  % Calculates the elements of the partition function: exp(-Ei(H)/kT)
  Z(:,ind_i) = exp(-beta .* E(ind_i));
  display(sprintf('Loop %1g of %1g completed in %0.2g min',id_i,length(iVs),toc/60));
end

for ind_T = 1:size(T,2)
  %sum_gamma = sum( (degen.*Z(ind_T,:)) + ndegen(ind_T,:) );
  %susceptibility(ind_T) = ( mu_B^2*beta(ind_T)/sum(Z(ind_T,:)) ) * sum_gamma;
  %susceptibility(ind_T) = sum( ((mu_B^2*beta(ind_T)/sum(Z(ind_T,:))).*degen.*Z(ind_T,:)) + ndegen(ind_T,:) );
  sum_gamma = sum( (beta(ind_T).*degen.*Z(ind_T,:)) ) + sum( ndegen(ind_T,:) );
  susceptibility(ind_T) = ( mu_B^2/sum(Z(ind_T,:)) ) * sum_gamma;
end

%susceptibility = susceptibility./mu_B;  % in u_B/T/atom


% Calculates the magnetisation M(ind_H,ind_T) per unit volume per atom;
%M = g * exp_JH;  % in u_B/atom
%M = g * mu_B * exp_JH .* N_A; % J/T/m^3/mol = A/m/mol
%To get magnetisation in emu/g (cgs units): Get magnetisation in J/T/m^3/kg by
%  M = g*exp_JH * mu_B * N_A/(molar mass). Then multiply by 4pi*10^{-7}.
%To get susceptibility in emu/mol (cgs units): Take M in J/T/m^3/mol, divide by
%  field in Tesla, and divide result by 4pi
%NB. H = B/mu_0 = B/(4pi*10^{-7}) where B is in Tesla. 
%  and chi = M/H not M/B => chi = M/(field in Tesla)*(4pi*10^{-7}) in m^3/mol.
%NB. 1T = 1Netwon/A/m 

