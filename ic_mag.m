function mag = ic_mag(Hic,H,Jdir,T,mu)

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
k_B   = 1.3806505e-23;     % J/K - Boltzmann constant
Q_e   = 1.60217653e-19;    % C - Charge of electron
N_A   = 6.0221415e23;      % Avogadro's number

% Normalises the direction vector
Jdir = Jdir ./ sqrt(Jdir*Jdir');

% Calculates the magnetic operator matrix
Jmat = mu{1}.*Jdir(1) + mu{2}.*Jdir(2) + mu{3}.*Jdir(3);

% Defining beta = 1/kT here saves computation later on.
beta = 1 ./ (k_B*T);

for ind_H = 1:size(H,2)
% Calculates the total Hamiltonian as a function of field (last index)
  Hmltn = Hcf + (-H(ind_H)*mu_B).*Jmat;

% Calculates the eigenvectors V and eigenvalues (enegies) E
% Where:            ---
%        | V  >  =  >    a  |j, j    >
%           i       ---i  i      z,i
%
  [V, Edummy] = eig(Hmltn);

% Reduce the energy levels to a vector.
  Edummy = sort( Edummy(logical(eye(size(Edummy,1)))) );
% Sets energy levels relative to lowest level.
  E = Edummy - min(Edummy);

% Converts energy levels from meV to J
  E = E .* (Q_e/1000);

  mj = -J:J;
  for ind_j = 1:size(V,2)
% Calculates the matrix elements <Vi|J.H|Vi>
    me = V(:,ind_j)' * Jmat * V(:,ind_j);

% Calculates the elements of the expectation of J.H: <i|J.H|i> exp(-Ei(H)/kT)
    JH(:,ind_j) = me * exp(-beta .* E(ind_j));

% Calculates the elements of the partition function: exp(-Ei(H)/kT)
    Z(:,ind_j) = exp(-beta .* E(ind_j));
  end

  for ind_T = 1:size(T,2)
% Calculates the expectation <<J.H>> = sum(<i|J.H|i>exp(-Ei/kT)) / sum(exp(-Ei/kT))
    exp_JH(ind_H,ind_T) = sum(JH(ind_T,:)) / sum(Z(ind_T,:));
  end

end

% Calculates the magnetisation M(ind_H,ind_T) per unit volume per atom;
M = exp_JH;  % in u_B/atom
%M = g * mu_B * exp_JH .* N_A; % J/T/m^3/mol = A/m/mol
%To get magnetisation in emu/g (cgs units): Get magnetisation in J/T/m^3/kg by
%  M = g*exp_JH * mu_B * N_A/(molar mass). Then multiply by 4pi*10^{-7}.
%To get susceptibility in emu/mol (cgs units): Take M in J/T/m^3/mol, divide by
%  field in Tesla, and divide result by 4pi
%NB. H = B/mu_0 = B/(4pi*10^{-7}) where B is in Tesla. 
%  and chi = M/H not M/B => chi = M/(field in Tesla)*(4pi*10^{-7}) in m^3/mol.
%NB. 1T = 1Netwon/A/m 
