function [H_cf,st_CF] = ic_hmltn(conf,F,xi,B,alpha,T,M,P,r_fl)
% Calculates the intermediate-coupling Hamiltonian after the method of Chan and Lam 1970
% At present, this function only works for the f-electron configurations.
%
% Syntax:
%                                                                                                 n
% Inputs:  conf - either: scalar - conf = n, where n is the number of electrons in configuration f
%                     or: vector - conf = [n l] giving the configuration l^n. l=3. n=1 to (2l+1)
%                     or: string - conf = 'f5' giving the l^n configuration in words. l=f
%          F    -         vector - F = [F0 F2 F4 F6] are the Slater integrals (in meV) which
%                                  determine the electrostatic interactions.
%          xi   -         scalar - the spin-orbit parameter (in meV)
%          B    -     cell array - B = {B2 B4 B6} are the Crystal Field parameters in Stevens
%                                  normalisation in meV.
%          r_fl -        logical - an optional flag to show that F = [E0 E1 E2 E3] are the Racah
%                                  parameters (linear combinations of Slater integrals) instead.

% By Duc Le 2007 - duc.le@ucl.ac.uk

icfact = [-sqrt(15/7)/2       sqrt(11/14)        -sqrt(429/7)/10  ];

if isnumeric(conf)
  if isscalar(conf)
    n = conf;
  elseif isvector(conf) && length(conf)==2 && conf(2)==3
    n = conf(1);
  else
    error('conf must be either a scalar or vector of length 2');
  end
elseif ischar(conf) && length(conf)==2 && conf(1)=='f'
  n = double(conf(2))-48;
else
  error('conf must be either a numeric scalar < 7, vector [n<7 l=3], or string "fn"');
end

% Checks n is integer and less than 7
if mod(n,1)~=0 | n>7 | n<0
  error('number of electrons n must be integer between 0 and 7');
end
% Checks F is correct format
if ~isnumeric(F) | length(F)~=4
  error('F must be a numeric vector of length 4');
end
% Checks xi
if ~isnumeric(xi) | ~isscalar(xi)
  error('xi must be a numeric scalar');
end
% Checks B
if iscell(B) & length(B)==3
  B2 = B{1}; B4 = B{2}; B6 = B{3};
  if ~isvector(B2) | length(B2)~=5 | ~isvector(B4) | length(B4)~=9 | ~isvector(B6) | length(B6)~=13
    error('B2=B{1},B4=B{2},B6=B{3} must be vectors of length 5, 9, and 13 respectively');
  end
else
  error('B must be a cell array of length 3')
end
% Checks alpha
if exist('alpha','var')
  if isscalar(alpha) || islogical(alpha)
    r_fl = alpha;
    clear alpha;
  elseif isvector(alpha) && length(alpha)==3
    cimat = racah_ci(n,3,alpha);
  else
    error('alpha must be a length-3 vector');
  end
end
% Checks T
if exist('T','var')
  if isscalar(T) || islogical(T)
    r_fl = T;
    clear T;
  elseif isvector(T) && (length(T)==6 || length(T)>=8)
    timat = judd_tichain(n,3,T);
  elseif isvector(T) && length(T)==3 && exist('M','var') && length(M)==4
    P = M; M = T; clear T;
  else
    error('T must be a length 6 (or >=8) vector');
  end
end
% Checks M,P
if exist('M','var')
  if isscalar(M) || islogical(M)
    r_fl = T;
    clear T;
  elseif isvector(M) && length(M)==3 && exist('P','var') && length(P)==4
    HssHsoo = judd_HssHsoo(n,3,M,P);
  elseif isvector(M) && length(M)==12
    HssHsoo = judd_HssHsoo(n,3,M);
    if exist('P','var') && (isscalar(P) || islogical(P))
      r_fl = P;
      clear P;
    end
  else
    error('You must specify both M (length 3) and P (length 4), or a single length 12 vector a_i');
  end
end
% Checks r_fl
if exist('r_fl')
  E = F;
else            % Converts Slater integrals into Racah parameters, after Racah IV
  E = [( F(1) - 10*F(2) -  33*F(3) -  286*F(4) );   ...
       (        70*F(2) + 231*F(3) + 2002*F(4) )/9; ...
       (           F(2) -   3*F(3) +    7*F(4) )/9; ...
       (         5*F(2) +   6*F(3) -   91*F(4) )/3];
end

% Calculates the electrostatic Hamiltonian
if exist('emat','file')==2
  load('emat');
  Emat = E(1).*s_emat{n,1} + E(2).*s_emat{n,2} + E(3).*s_emat{n,3} + E(4).*s_emat{n,4};
else
  Emat = racah_emat(n,3,E);
end

% Calculates the Configuration Interaction Hamiltonian
if exist('cimat')
  Emat = Emat + cimat;         % The two-electron operators: alpha*L^2 + beta*G(R7) + gamma*G(G2)
end
if exist('timat')
  Emat = Emat + timat;         % The three electron operators: sum_i(T^i*t_i) i=2,3,4,6,7,8 
end

Emat = sparse(Emat);

% Because the Spin-Orbit and Crystal Field Hamiltonian matrices take a very long time to compute, they
% have been pre-computed and stored in a file: UVmatrices.mat
% If however this file does not exist, it must be generated again using the genmatrices script.
%if ~exist('UVmatrics.mat','file')
%  genmatrices
%end

% Calcalates the spin-orbit Hamiltonian
if n<5
  matfile = ['UVmat' sprintf('%1g',n) '.mat'];
else
  matfile = ['UVmat' sprintf('%1g',n) 'SO.mat'];
end
if exist(matfile,'file')==2
  load(matfile);
  H_SO = H_SO.*xi;
else
  display(sprintf('Calculating Spin-Orbit Matrix for n=%1g',n)); tic
  [H_SO,st_SO] = racah_so(n,3,1);
  if n<5
    display(sprintf('Time elapsed = %0.5g min',toc/60));
    display(sprintf('Calculating CF Matrix for k=2')); tic; U2 = fast_ukq(n,3,2); 
    display(sprintf('Time elapsed = %0.5g min',toc/60));
    display(sprintf('Calculating CF Matrix for k=4')); tic; U4 = fast_ukq(n,3,4);
    display(sprintf('Time elapsed = %0.5g min',toc/60));
    display(sprintf('Calculating CF Matrix for k=6')); tic; U6 = fast_ukq(n,3,6);
    display(sprintf('Time elapsed = %0.5g min',toc/60));
    save(matfile,'U2','U4','U6','H_SO','st_SO');
  else
    save(matfile,'H_SO','st_SO');
  end
  H_SO = H_SO.*xi;
end
% Re-expresses the electrostatic energy matrix from the |vUSL> basis to the |vUSLJ> basis
lstSO = length(st_SO); H_El = sparse(lstSO,lstSO);
for i = 1:lstSO
  for j = 1:lstSO
    if (st_SO{i}{5}==st_SO{j}{5}) % J==J'
      H_El(i,j) = Emat(st_SO{i}{6},st_SO{j}{6});
    end
  end
end
% Calculates the Russell-Saunders (LS) coupling basis states from the Electrostatic+SO Hamiltonians.
H_LS = H_El + H_SO;
%[Vls,Els] = eig(H_LS);

% Adds in the spin-spin, spin-other-orbit, and spin-orbit configuration interactions if applicable
if exist('HssHsoo')
  H_LS = H_LS + HssHsoo;
end

%H_cf = H_LS; st_CF = st_SO; if 0
% Calculates the projection operator matrix, after Chan and Lam for the first 6 J-multiplets
%Els = diag(Els) - min(min(Els)); Els_sort = sort(Els); P = zeros(size(H_LS));
%for i = 1:6
%  st_trunc(i) = find(Els_sort(i)==Els);
%  phi_r = Vls(:,st_trunc(i));
%  P = P + (phi_r * phi_r');
%end

% Re-expresses the LS matrix from the |vUSLJ> basis to the |vUSLJm> basis, same as the CF matrix.
index = 0;
for i = 1:length(st_SO)
  J = st_SO{i}{5};
  for mJ = -J:J
    index = index + 1;
    i_ic(index) = i;
    st_CF{index} = {st_SO{i}{1:5} mJ};
  end
end

% Calculates the crystal field Hamiltonian
% NB. Times to calculate the CF matrices are: n=2,B2: 27.305s; n=3,B2: 483.62s; n=4,B2: 2957s; n=5,B2: 9802.8s
%     Size of CF matrices (NxN): n=2,B2: 91; n=3,B2: 364; n=4,B2: 1001; n=5,B2: 2002;
H_cf = sparse(length(st_CF),length(st_CF));
if n<5
  matfile = ['UVmat' sprintf('%1g',n) '.mat'];
  load(matfile);
  U = {U2 U4 U6};
  for k = 1:3
    for q = 1:(4*k+1)
      H_cf = H_cf + B{k}(q) .* U{k}{q};
    end
  end
else
%  % Loading and calculating all values of k and q. Elasped time for n=5: 438.13s
%  matfile = ['UVmat' sprintf('%1g',n) 'U2.mat'];
%  if exist(matfile,'file')==2; load(matfile); else; U2 = racah_Ukq(n,3,2); end
%  for q = 1:(4*1+1); H_cf = H_cf + B{1}(q) .* U2{q}; end; clear U2;
%
%  matfile = ['UVmat' sprintf('%1g',n) 'U4.mat'];
%  if exist(matfile,'file')==2; load(matfile); else; U4 = racah_Ukq(n,3,4); end
%  for q = 1:(4*2+1); H_cf = H_cf + B{2}(q) .* U4{q}; end; clear U4;
%
%  matfile = ['UVmat' sprintf('%1g',n) 'U6.mat'];
%  if exist(matfile,'file')==2; load(matfile); else; U6 = racah_Ukq(n,3,6); end
%  for q = 1:(4*3+1); H_cf = H_cf + B{3}(q) .* U6{q}; end; clear U6;

  for k = 1:3
    for q = 1:(4*k+1)
      if B{k}(q) ~= 0 
        if q<2*k; matfile = sprintf('UVmat%1gU%1gm%1g.mat',n,2*k,-(q-1-(2*k)));
        else;     matfile = sprintf('UVmat%1gU%1gk%1g.mat',n,2*k, (q-1-(2*k)));
	end
        if exist(matfile,'file')==2
          load(matfile); 
	  U = U./icfact(k);
        else
	  display(sprintf('Calculating CF Matrix k=%1g,q=%1g',2*k,q-1-(2*k))); tic;
          U = fast_ukq(n,3,2*k,q-1-(2*k));
          display(sprintf('Time elapsed = %0.5g min',toc/60));
          save(matfile,'U');
        end
        H_cf = H_cf + B{k}(q) .* U;
%        H_cf = H_cf + B{k}(q) .* U{q};
	% Time to calculate matrix: n=5, t=100s
        clear U;
      end
    end
  end

end

% Re-expresses the LS matrix from the |vUSLJ> basis to the |vUSLJm> basis, same as the CF matrix.
for i = 1:length(st_CF)
  for j = 1:length(st_CF)
    if (st_CF{i}{6}==st_CF{j}{6})   % mJ==mJ'
      H_cf(i,j) = H_LS(i_ic(i),i_ic(j)) + H_cf(i,j);
    end
  end
end

%end
% Calculates the truncated LS-coupling basis states.
%[V,E] = eig(P * (H_LS + H_cf) * P');

% NB. Times for MATLAB to diagonalised NxN matrices: N=500, t=1.94888s; N=1000, t=17.714922s; N=5000, t=3221.718853s

% NB. for n=5, time to set up and diagonalise H_ic from loading Ukq, is 726.75s
