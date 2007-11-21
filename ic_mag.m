function mag = ic_mag(Hic,H,Jdir,T,mu)

%if iscell(Hic)
%  V = Hic{1}; E = Hic{2};
%  if ~isvector(E)
%    E = diag(E);
%  end
%  if (length(E) ~= size(V,1))
%    error('V and E must have same size or length');
%  end
%  E = real(E) - min(real(E));
%else
%  [V,E] = eig(Hic);
%  E = diag(E);
%  E = E - min(E);
%end

% Physical constants. Taken from NIST Reference on Constants, Units, and 
% Uncertainty, http://physics.nist.gov/cuu/Constants/
%mu_B  = 927.400949e-26;    % J/Tesla - Bohr magneton
%Hic = Hic .* 0.1239842436; % Converts to meV. 
%mu_B = 5.78838263e-2;    % meV/T - Bohr magneton
mu_B  = 0.46686437;        % cm^{-1} / Tesla - Bohr magneton
k_B   = 1.3806505e-23;     % J/K - Boltzmann constant
Q_e   = 1.60217653e-19;    % C - Charge of electron
N_A   = 6.0221415e23;      % Avogadro's number
h     = 6.62606896e-34;    % Js - Planck's constant
c     = 299792458;         % m/s - speed of light in vacuum

% Normalises the direction vector
Jdir = Jdir ./ sqrt(Jdir*Jdir');

% Calculates the magnetic operator matrix
Jmat = mu{1}.*Jdir(1) + mu{2}.*Jdir(2) + mu{3}.*Jdir(3);

% Defining beta = 1/kT here saves computation later on.
beta = 1 ./ (k_B*T);

for ind_H = 1:size(H,2)
% Calculates the total Hamiltonian as a function of field (last index)
  Hmltn = Hic + (-H(ind_H)*mu_B).*Jmat;

% Calculates the eigenvectors V and eigenvalues (enegies) E
% Where:            ---
%        | V  >  =  >    a  |j, j    >
%           i       ---i  i      z,i
%
  [V, E] = eig(Hmltn);

% Reduce the energy levels to a vector.
  E = real(diag(E)); 
% Sets energy levels relative to lowest level.
  E = E - min(E);

% Converts energy levels from meV to J
%  E = E .* (Q_e/1000);

% Converts energy levels from cm^{-1} to J - NB. E = hf = hc/lambda in metres!
  E = E .* (h*c*100);

  for ind_j = 1:length(E)
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
mag = exp_JH;  % in u_B/atom
%M = g * mu_B * exp_JH .* N_A; % J/T/m^3/mol = A/m/mol
%To get magnetisation in emu/g (cgs units): Get magnetisation in J/T/m^3/kg by
%  M = g*exp_JH * mu_B * N_A/(molar mass). Then multiply by 4pi*10^{-7}.
%To get susceptibility in emu/mol (cgs units): Take M in J/T/m^3/mol, divide by
%  field in Tesla, and divide result by 4pi
%NB. H = B/mu_0 = B/(4pi*10^{-7}) where B is in Tesla. 
%  and chi = M/H not M/B => chi = M/(field in Tesla)*(4pi*10^{-7}) in m^3/mol.
%NB. 1T = 1Netwon/A/m 
