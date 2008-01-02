function I = ic_pes(conf,F,xi,alpha,r_fl)
% Calculates the photo-emmission spectrum under intermediate-coupling, after Gerken, 1983
%
% Syntax:  I = ic_pes(conf,F,xi,alpha,r_fl)
%                                                                                                  n
% Inputs:  conf  - either: scalar - conf = n, where n is the number of electrons in configuration f
%                      or: vector - conf = [n l] giving the configuration l^n. l=3. n=1 to (2l+1)
%                      or: string - conf = 'f5' giving the l^n configuration in words. l=f
%          F     -         vector - F = [F0 F2 F4 F6] are the Slater integrals (in meV) which
%                                   determine the electrostatic interactions.
%          xi    -         scalar - the spin-orbit parameter (in meV)
%          B     -     cell array - B = {B2 B4 B6} are the Crystal Field parameters in Stevens
%                                   normalisation in meV.
%          alpha -         vector - [alpha beta gamma] of the parameters for the two-electron
%                                   configuration interaction
%          r_fl  -        logical - an optional flag to show that F = [E0 E1 E2 E3] are the Racah
%                                   parameters (linear combinations of Slater integrals) instead.

% By Duc Le 2007 - duc.le@ucl.ac.uk

if isnumeric(conf)
  if isscalar(conf)
    n = conf;
  elseif isvector(conf) && length(conf)==2 && conf(2)==3
    n = conf(1);
  else
    error('conf must be either a scalar or vector of length 2');
  end
elseif ischar(conf) 
  if length(conf)==2 && conf(1)=='f'
    n = double(conf(2))-48;
  elseif nargin==1
    [E,xi,alpha,n] = parlookup(conf);
    [Ep,xip,alphap] = parlookup(conf,1);
  end
else
  error('conf must be either a numeric scalar < 7, vector [n<7 l=3], or string "fn"');
end

% Checks n is integer and less than 7
if mod(n,1)~=0 | n>7 | n<2
  error('number of electrons n must be integer between 2 and 7');
end
% Checks F is correct format
if exist('F','var') && (~isnumeric(F) | length(F)~=4)
  error('F must be a numeric vector of length 4');
end
% Checks xi
if ~isnumeric(xi) | ~isscalar(xi)
  error('xi must be a numeric scalar');
end
% Checks B
%if iscell(B) & length(B)==3
%  B2 = B{1}; B4 = B{2}; B6 = B{3};
%  if ~isvector(B2) | length(B2)~=5 | ~isvector(B4) | length(B4)~=9 | ~isvector(B6) | length(B6)~=13
%    error('B2=B{1},B4=B{2},B6=B{3} must be vectors of length 5, 9, and 13 respectively');
%  end
%else
%  error('B must be a cell array of length 3')
%end
cimat = 0;
% Checks alpha
if exist('alpha','var')
  if isscalar(alpha) || islogical(alpha)
    r_fl = alpha;
    clear alpha;
  elseif isvector(alpha) && length(alpha)==3
    cimatn   = racah_ci(n,3,alpha);
    cimatnm1 = racah_ci(n-1,3,alphap);
  else
    error('alpha must be a length-3 vector');
  end
end
% Checks r_fl
if ~exist('E','var')
  if exist('r_fl','var')
    E = F;
  else            % Converts Slater integrals into Racah parameters, after Racah IV
    E = [( F(1) - 10*F(2) -  33*F(3) -  286*F(4) );   ...
         (        70*F(2) + 231*F(3) + 2002*F(4) )/9; ...
         (           F(2) -   3*F(3) +    7*F(4) )/9; ...
         (         5*F(2) +   6*F(3) -   91*F(4) )/3];
  end
end

% Calculates the IC Hamiltonian for fn and f(n-1) configurations
%Hic_n   = ic_hmltn(n,E,xi,B,alpha,1);
%Hic_nm1 = ic_hmltn(n-1,E,xi,B,alpha,1);
%[Vn,En]     = eig(Hic_n);
%[Vnm1,Enm1] = eig(Hic_nm1);

[Vn,En,st_n]       = diagicmat(n,E,xi,cimatn);       En   = diag(En)   - min(diag(En));
[Vnm1,Enm1,st_nm1] = diagicmat(n-1,Ep,xip,cimatnm1); Enm1 = diag(Enm1) - min(diag(Enm1));

% Transitions only between states of same J
iE = find(En==0); Cn = Vn(:,iE); J = st_n{find(abs(Vn(:,iE))==max(abs(Vn(:,iE))))}{5};
% Calculates the intensity for transitions from the ground state - min(En) to any Enm1
I = zeros(length(st_nm1),2);
for iE = 1:length(st_nm1)
  I(iE,1) = Enm1(iE)/8066; I(iE,2) = 0; Jp = st_nm1{find(abs(Vnm1(:,iE))==max(abs(Vnm1(:,iE))))}{5};
  for j = (5/2):(7/2)
    Isum = 0;
    for in = 1:length(st_n)
      if st_n{in}{5} == J
        Cjtsl = Cn(in); 
        S = st_n{in}{1}; Ls = st_n{in}{2}; L = racah_lconv(Ls); v = st_n{in}{3}; U = st_n{in}{4};
        for im = 1:length(st_nm1)
	  if st_nm1{im}{5} == Jp
            Cjtslp = Vnm1(im,iE); Sp = st_nm1{im}{1}; Lp = racah_lconv(st_nm1{im}{2});
            Isum = Isum + ( Cjtsl*Cjtslp * sqrt((2*S+1)*(2*L+1)) * ninej([1/2 3 j; S L J; Sp Lp Jp]) ...
                          * racah_cfp(n,U,v,S,Ls,st_nm1{im}{4},st_nm1{im}{3},Sp,st_nm1{im}{2}) );
          end
        end
      end  % If st_n{in}{5}==J
    end
    I(iE,2) = I(iE,2) + ( (2*j+1) * Isum^2 );
  end  % for j
  I(iE,2) = n*(2*Jp+1) * I(iE,2);
end
 
%-------------------------------------------------------------------------------------------%
%                   Function to calculate and diagonalise the IC matrix                     %
%-------------------------------------------------------------------------------------------%
function [Vls,Els,st_SO] = diagicmat(n,E,xi,cimat)
% Calculates the electrostatic Hamiltonian
if exist('emat','file')==2
  load('emat');
  Emat = E(1).*s_emat{n,1} + E(2).*s_emat{n,2} + E(3).*s_emat{n,3} + E(4).*s_emat{n,4};
else
  Emat = racah_emat(n,3,E);
end
% Calculates the configuration interaction Hamiltonian
if ~isscalar(cimat);
  Emat = Emat + cimat;
end
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
  display(sprintf('Time elapsed = %0.5g min',toc/60));
  save(matfile,'H_SO','st_SO');
  H_SO = H_SO.*xi;
end
% Re-expresses the electrostatic energy matrix from the |vUSL> basis to the |vUSLJ> basis
for i = 1:length(st_SO)
  for j = 1:length(st_SO)
    if (st_SO{i}{5}==st_SO{j}{5}) % J==J'
      H_El(i,j) = Emat(st_SO{i}{6},st_SO{j}{6});
    end
  end
end
% Calculates the Russell-Saunders (LS) coupling basis states from the Electrostatic+SO Hamiltonians.
H_LS = H_El + H_SO;
[Vls,Els] = eig(H_LS);

%-------------------------------------------------------------------------------------------%
%                 Function to lookup the E, spin-orbit and CI parameters                    %
%-------------------------------------------------------------------------------------------%
function [E,xi,alpha,n] = parlookup(conf,p)
switch conf
  case 'Pr3+'; fourf = 2;
  case 'Nd3+'; fourf = 3;
  case 'Pm3+'; fourf = 4;
  case 'Sm3+'; fourf = 5;
  case 'Eu3+'; fourf = 6;
  case 'Gd3+'; fourf = 7;
  case 'Tb3+'; fourf = 8;
  case 'Dy3+'; fourf = 9;
  case 'Ho3+'; fourf = 10;
  case 'Er3+'; fourf = 11;
  case 'Tm3+'; fourf = 12;
  case 'Pa3+'; fivef = 2;
  case 'U3+' ; fivef = 3;
  case 'Np3+'; fivef = 4;
  case 'Pu3+'; fivef = 5;
  case 'Am3+'; fivef = 6;
  case 'Cm3+'; fivef = 7;
  case 'Bk3+'; fivef = 8;
  case 'Cf3+'; fivef = 9;
  case 'Es3+'; fivef = 10;
  case 'Fm3+'; fivef = 11;
  case 'Md3+'; fivef = 12;
end
%          E^1    E^2    E^2      xi     alpha  beta    gamma
p4f = {[0 4548.2 21.937 466.73]  740.75 [21.255   -799.94 1342.9  ];   % Pr3+
       [0 4739.3 23.999 485.96]  884.58 [ 0.5611  -117.15 1321.3  ];   % Nd3+
       [0 4921.6 24.522 525.53] 1000.8  [10.991   -244.88  789.74 ];   % Pm3+
       [0 5496.9 25.809 556.40] 1157.3  [22.250   -742.55  796.64 ];   % Sm3+
       [0 5573.0 26.708 557.39] 1326.0  [25.363   -580.25 1155.7  ];   % Eu3+
       [0 5761.0 28.02  582.0 ] 1450.0  [22.55    -103.7   997.0  ];   % Gd3+
       [0 6021.5 29.03  608.54] 1709.5  [20.131   -370.7  1255.9  ];   % Tb3+
       [0 6119.6 30.012 610.14] 1932.0  [37.062  -1139.1  2395.3  ];   % Dy3+
       [0 6440.6 30.22  624.39] 2141.3  [23.635   -807.2  1278.4  ];   % Ho3+
       [0 6769.9 32.388 646.62] 2380.7  [18.347   -509.28  649.71 ];   % Er3+
       [0 7142.4 33.795 674.27] 2628.7  [14.677   -631.79    0    ]};  % Tm3+
%          E^1    E^2    E^2      xi     alpha  beta    gamma
p5f = {[0 2650   11.3   245.4 ] 1350    [28       -730    1100    ];   % Pa3+ (Extrapolated)
       [0 2878.5 11.848 258.25] 1623.0  [27.6     -722    1000    ];   % U3+  (Crosswhite JChemPhys 72 5103 1980)
       [0 3168.9 13.679 297.75] 1932.6  [32.713   -756.12 1084.2  ];   % Np3+ (Carnall JChemPhys 61 4993 1974)
       [0 3634.5 15.356 342.15] 2272.2  [31.064   -675.0    29.743];   % Pu3+ (Carnall JChemPhys 53 2922 1970)
       [0 3582.8 17.276 334.30] 2593.3  [21.634   -158.48 1240.4  ];   % Am3+ (Pappalardo JChemPhys 51 1182 1969)
       [0 3955.4 17.789 369.79] 2876.1  [27.895   -925.21 1119.6  ];   % Cm3+ (Carnall JChemPhys 63 3510 1975)
       [0 4127.2 19.537 387.83] 3252.8  [31.123  -1283.7  1247.5  ];   % Bk3+ (Carnall JChemPhys 58 3614 1973c)
       [0 4173.2 19.84  391.40] 3601.7  [35.371   -748.02 1200.0  ];   % Cf3+ (Carnall JChemPhys 58 1938 1973b)
       [0 4445.9 21.446 418.64] 4014.7  [22.505   -722.53 1000.0  ];   % Es3+ (Carnall JChemPhys 59 1785 1973a)
       [0 4800   22.9   453   ] 4450    [23       -730    1100    ];   % Fm3+ (Extrapolated)
       [0 5000   23.8   472   ] 4900    [23       -730    1100    ]};  % Md3+ (Extrapolated)
% Parses parameters for output
if exist('fourf')
  n = fourf;
  if nargin==1
    E     = p4f{n,1};
    xi    = p4f{n,2};
    alpha = p4f{n,3};
  else
    E     = p4f{n-1,1};
    xi    = p4f{n-1,2};
    alpha = p4f{n-1,3};
  end
elseif exist('fivef')
  n = fivef;
  if nargin==1
    E     = p5f{n,1};
    xi    = p5f{n,2};
    alpha = p5f{n,3};
  else
    E     = p5f{n-1,1};
    xi    = p5f{n-1,2};
    alpha = p5f{n-1,3};
  end
else
  error('Ion string not recognised');
end
