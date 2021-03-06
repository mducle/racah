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
if mod(n,1)~=0 | n>14 | n<0
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
elseif iscell(B) && length(B)==4 && ischar(B{4})
  if strncmp(B{4},'s',1)
    B = B(1:3); B = wy2stev(B,'s');
  else
    B = B(1:3);
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
    if n<3; warning('3-particle CI interaction requires n>=3. n=%1g, so ignoring T^i parameters',n); else
    timat = judd_tichain(n,3,T); end
  elseif isvector(T) && length(T)==3 && exist('M','var') && length(M)==4
    P = M; M = T; clear T;
  else
    error('T must be a length 6 (or >=8) vector');
  end
end
% Checks M,P
if exist('M','var')
  if isscalar(M) || islogical(M)
    r_fl = M;
    clear M;
  %elseif isvector(M) && length(M)==3 && exist('P','var') && length(P)==4
  %  HssHsoo = judd_HssHsoo(n,3,M,P);
  elseif isvector(M) 
    if length(M)==3; 
      if exist('P','var') && length(P)==4
        %HssHsoo = judd_HssHsoo(n,3,M,P);
	HssHsoo = crosswhite_HssHsoo(n,3,M,P(2:4));
      elseif exist('P','var') && length(P)==3
        %HssHsoo = judd_HssHsoo(n,3,M,[0 P]); 
	HssHsoo = crosswhite_HssHsoo(n,3,M,P);
      else
        %HssHsoo = judd_HssHsoo(n,3,M,[0 0 0 0]); 
	HssHsoo = crosswhite_HssHsoo(n,3,M,[0 0 0]);
      end
    elseif length(M)==4; 
      P = M; clear M; 
    else
      error('M must be a length 3 or 4 vector'); 
    end
  %elseif isvector(M) && length(M)==12
  %  HssHsoo = judd_HssHsoo(n,3,M);
  %  if exist('P','var') && (isscalar(P) || islogical(P))
  %    r_fl = P;
  %    clear P;
  %  end
  %else
  %  error('You must specify both M (length 3) and P (length 4), or a single length 12 vector a_i');
  %end
  else
    error('M must be a vector');
  end
end
if exist('P','var')
  if isscalar(P) || islogical(P)
    r_fl = P;
    clear P;
  elseif isvector(P)
    if length(P)==4
      %HssHsoo = judd_HssHsoo(n,3,[0 0 0],P);
      Pmat = crosswhite_HssHsoo(n,3,M,P(2:4));
    elseif length(P)==3
      Pmat = crosswhite_HssHsoo(n,3,M,P);
      %HssHsoo = judd_HssHsoo(n,3,[0 0 0],[0 P]);
    end
    if exist('HssHsoo','var'); HssHsoo = HssHsoo + Pmat; else; HssHsoo = Pmat; end
  else
    error('P must be a length 4 vector');
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
%lstSO = length(st_SO); H_El = sparse(lstSO,lstSO);
%for i = 1:lstSO
%  for j = 1:lstSO
%    if (st_SO{i}{5}==st_SO{j}{5}) % J==J'
%      H_El(i,j) = Emat(st_SO{i}{6},st_SO{j}{6});
%    end
%  end
%end
  st_SI = racah_states(n,3); id = 0; icv = 0; icv1 = 0; icv2 = 0;
  for i = 1:length(st_SI)
    L = racah_lconv(st_SI{i}{2});
    S = st_SI{i}{1};
    Jmin = abs(L-S); Jmax = L+S;
    cvSI2SO(i,1) = icv1+1; icv1=icv1+(Jmax-Jmin+1); cvSI2SO(i,2) = icv1;      % Keeps a running index for the convH2H subroutine
    cvSI2CF(i,1) = icv2+1; icv2=icv2+sum(2.*[Jmin:Jmax]+1);                   %   first to convert from |vUSL> to |vUSLJ>
    cvSI2CF(i,2) = icv2;                                                      % This time for conversion from |vUSL> to |vUSLJmJ>
    for J = Jmin:Jmax
      id = id + 1;
      cvSO2CF(id,1) = icv+1; icv=icv+(2*J+1); cvSO2CF(id,2) = icv;            % Keeps a running index for the convH2H subroutine
    end
  end
% Calculates the Russell-Saunders (LS) coupling basis states from the Electrostatic+SO Hamiltonians.
H_LS = convH2H(Emat,length(H_SO),cvSI2SO) + H_SO;

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
%index = 0;
%for i = 1:length(st_SO)
%  J = st_SO{i}{5};
%  for mJ = -J:J
%    index = index + 1;
%    i_ic(index) = i;
%    st_CF{index} = {st_SO{i}{1:5} mJ};
%  end
%end

% Calculates the crystal field Hamiltonian
%
%         ---   k    [  k        q  k    ]     ---  k  k     ---   k [  k       q  k  ]
% V   = i >    B     | O   - (-1)  O     |  +  >   B  O   +  >    B  | O  + (-1)  O   |
%  cf     ---   -|q| [  |q|         -|q| ]     ---  0  0     ---   q [  q          -q ]
%        k,q<0                                  k           k,q>0  
%
%H_cf = sparse(length(st_CF),length(st_CF));
H_cf = sparse(icv2,icv2);
if n<5
  matfile = ['UVmat' sprintf('%1g',n) '.mat'];
  load(matfile);
  U = {U2 U4 U6};
  for k = 1:3
    %for q = 1:(4*k+1)
    %  H_cf = H_cf + B{k}(q) .* (U{k}{q}./icfact(k));
    for iq = 1:(4*k+1)
      q = iq-1-(2*k); qp = 1+(2*k)+abs(q); qm = 1+(2*k)-abs(q);
      H_cf = H_cf + ( (B{k}(iq)/icfact(k)) .* (U{k}{qp} + (sign(q)*(-1)^q).*U{k}{qm}) ); 
    end
  end
else

  for k = 1:3
    %for q = 1:(4*k+1)
    for iq = 1:(4*k+1)
      q = iq-1-(2*k); qp = 1+(2*k)+abs(q); qm = 1+(2*k)-abs(q);
      %if B{k}(q) ~= 0 
      %  if q<2*k; matfile = sprintf('UVmat%1gU%1gm%1g.mat',n,2*k,-(q-1-(2*k)));
      %  else;     matfile = sprintf('UVmat%1gU%1gk%1g.mat',n,2*k, (q-1-(2*k))); end
      %  if exist(matfile,'file')==2
      %    load(matfile); 
      %  else
      %    display(sprintf('Calculating CF Matrix k=%1g,q=%1g',2*k,q-1-(2*k))); tic;
      %    U = fast_ukq(n,3,2*k,q-1-(2*k));
      %    display(sprintf('Time elapsed = %0.5g min',toc/60));
      %    save(matfile,'U');
      %  end
      %  H_cf = H_cf + B{k}(q) .* (U./icfact(k));
      %  clear U;
      %end
      if B{k}(iq) ~= 0
        matfileP = sprintf('UVmat%1gU%1gk%1g.mat',n,2*k,q);
        matfileM = sprintf('UVmat%1gU%1gm%1g.mat',n,2*k,q);
        if exist(matfileP,'file')==2
          load(matfileP); 
        else
          display(sprintf('Calculating CF Matrix k=%1g,q=%1g',2*k,abs(q))); tic;
          [U,st_CF] = fast_ukq(n,3,2*k,abs(q));
          display(sprintf('Time elapsed = %0.5g min',toc/60));
          save(matfileP,'U','st_CF');
	end
	H_cf = H_cf + B{k}(iq) .* (U./icfact(k));
	clear U;
	if q~=0
          if exist(matfileM,'file')==2
            load(matfileM); 
          else
            display(sprintf('Calculating CF Matrix k=%1g,q=%1g',2*k,-abs(q))); tic;
            [U,st_CF] = fast_ukq(n,3,2*k,-abs(q));
            display(sprintf('Time elapsed = %0.5g min',toc/60));
            save(matfileM,'U','st_CF');
          end
          H_cf = H_cf + B{k}(iq) .* (sign(q)*(-1)^q).*(U./icfact(k));
	  clear U;
	end
      end
    end
  end

end

% Re-expresses the LS matrix from the |vUSLJ> basis to the |vUSLJm> basis, same as the CF matrix.
%for i = 1:length(st_CF)
%  for j = 1:length(st_CF)
%    if (st_CF{i}{6}==st_CF{j}{6})   % mJ==mJ'
%      H_cf(i,j) = H_LS(i_ic(i),i_ic(j)) + H_cf(i,j);
%    end
%  end
%end
H_cf = convH2H(H_LS,length(H_cf),cvSO2CF) + H_cf;

% ------------------------------------------------------------------------------------------------------------------------------- %
% Routine to convert a matrix from one basis to another, e.g from |vUSL> to |vUSLJmJ>
% ------------------------------------------------------------------------------------------------------------------------------- %
function Hout = convH2H(Hin,lnOut,cv)
  Hout = sparse(lnOut,lnOut);
  [inz,jnz] = find(Hin);
  for id = 1:length(inz)
    i = inz(id); j = jnz(id);
    m = cv(i,2)-cv(i,1)+1; n = cv(j,2)-cv(j,1)+1;
    Hout(cv(i,1):cv(i,2),cv(j,1):cv(j,2)) = spdiags(ones(min(n,m),1).*Hin(i,j),0,n,m);   % Elements only on the diagonals
  end
