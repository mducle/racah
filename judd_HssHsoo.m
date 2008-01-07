function H = judd_HssHsoo(n,l,M,P,cellout)
% Calculates the spin-spin and spin-other-orbit Hamiltonian from M^k, P^k parameters
%
% Syntax:  H = judd_HssHsoo(n,l,M,P,c)
%
% Inputs:  n - scalar - number of equivalent electrons
%          l - scalar - orbital angular momentum quantum number of the electrons
%          M - vector - either a 1x3 vector [M0 M2 M4] of Marvin's integrals 
%                     - or a 1x12 vector of the parameters a_i of Crosswhite and Judd
%          P - vector - a 1x4 vector [P0 P2 P4 P6] of the spin-other-orbit/CI parameters
%                       (if M the 1x12 a_i parameters, this input is ignored)
%          c - any    - if a fifth input is given, the function will assume you want
%                       a cell array of the z_i operators in the |vUSLJ> basis
%
% Outputs: H - matrix - the Hamiltonian in the |vUSLJ> basis, labeled as per Nielson and Koster
%              cell   - a cell array {H z1 z2 ... z11 zc} of both Hamiltonian and zi operators

% By Duc Le 2007 - duc.le@ucl.ac.uk

%zcnt = 14;

if ~isnumeric(n) || ~isnumeric(l) || ~isscalar(n) || ~isscalar(l)
  error('n and l must be numerical scalars');
%elseif isnumeric(M) && isvector(M) 
%  if length(M)==zcnt
%    a = M;
%  elseif isnumeric(P) && isvector(P) && length(M)==3 && length(P)==4
%    a = judd_MPtoa(M,P);
%  else
%    error('M and P must be numerical vectors of length 3 and 4 respectively');
%  end
elseif ~isnumeric(M) || length(M)~=3
  error('M must be a numerical vector of length 3');
elseif ~isnumeric(P) || length(P)~=4
  error('P must be a numerical vector of length 4');
end

%zmat = {};
%if exist('HssHsoo.mat','file')==2
%  load 'HssHsoo.mat'
%  if ~isempty(zmat)
%    zmatstruct = zmat;
%    zmat = {};
%  end
%else
%  save('HssHsoo.mat','zmat');
%end

mk = {}; pk = {};
if exist('HssHsoo.mat','file')==2
  load 'HssHsoo.mat'
  if ~isempty(mk); mkstruct = mk; mk = {}; end
  if ~isempty(pk); pkstruct = pk; pk = {}; end
else
  save('HssHsoo.mat','mk','pk');
end

confstring = sprintf('%c%1g',lower(racah_lconv(l)),n);

%if isfield(mk,confstring)
%  H = sparse(length(zmat.(confstring)),length(zmat.(confstring)));
%  for z = 1:zcnt
%    H = H + a(z).*zmat.(confstring){z};
%  end
%else
if isfield(mk,confstring) && isfield(pk,confstring); lmk = length(mk.(confstring)); H = sparse(lmk,lmk); 
  for k = 1:3; H = H + M(k).*(mk.(confstring){k}+mk.(confstring){k+3}); end
  for k = 1:4; H = H + P(k).*pk.(confstring){k}; end
else

%-----------------------------------------------------------------------------------------------------%
%if ~exist('T22') 
%  display('Calculating z1-z4'); T22 = judd_T22(n,l);
%  save('HssHsoo.mat','T22','-append');
%elseif ~isfield(T22,confstring);
%  confcalc = fieldnames(T22); lastconf = confcalc{length(confcalc)}; nold = sscanf(lastconf,'%c%1g');
%  display('Calculating z1-z4'); T22 = judd_T22(n,l,nold(2),T22);
%  save('HssHsoo.mat','T22','-append');
%end
if ~exist('rpk','var') 
  display('Calculating reduced matrix elements of the coefficients of P^k'); rpk = judd_pk(n,l);
  save('HssHsoo.mat','rpk','-append');
elseif ~isfield(rpk,confstring);
  confcalc = fieldnames(rpk); lastconf = confcalc{length(confcalc)}; nold = sscanf(lastconf,'%c%1g');
  display('Calculating reduced matrix elements of the coefficients of P^k'); rpk = judd_pk(n,l,nold(2),rpk);
  save('HssHsoo.mat','rpk','-append');
end
%if ~exist('T11') 
%  display('Calculating z5-z11,zc'); T11 = judd_T11(n,l);
%  save('HssHsoo.mat','T11','-append');
%elseif ~isfield(T11,confstring);
%  confcalc = fieldnames(T11); lastconf = confcalc{length(confcalc)}; nold = sscanf(lastconf,'%c%1g');
%  display('Calculating z5-z11,zc'); T11 = judd_T11(n,l,nold(2),T11);
%  save('HssHsoo.mat','T11','-append');
%end
if n<4
  if ~exist('rmk','var') 
    display('Calculating reduced matrix elements of the coefficients of M^k'); rmk = judd_mk(n,l);
    save('HssHsoo.mat','rmk','-append');
  elseif ~isfield(rmk,confstring);
    confcalc = fieldnames(rmk); lastconf = confcalc{length(confcalc)}; nold = sscanf(lastconf,'%c%1g');
    display('Calculating reduced matrix elements of the coefficients of M^k'); rmk = judd_mk(n,l,nold(2),rmk);
    save('HssHsoo.mat','rmk','-append');
  end
else
  if ~exist('rmk','var') || ~isfield(rmk,confstring)
    display('Calculating reduced matrix elements of the coefficients of M^k');
    [tmp_mkmat,tmpHcell] = horie_mk(n,l); rmk.(confstring) = tmpHcell([2 5]);
    save('HssHsoo.mat','rmk','-append');
  end
end

%redmat = {T22.(confstring){:} T11.(confstring){:}};
redmat = {rmk.(confstring){1}{:} rmk.(confstring){2}{:} rpk.(confstring){:}};

statesLS = racah_states(n,l); index = 0;
for i = 1:length(statesLS)
  L = racah_lconv(statesLS{i}{2});
  S = statesLS{i}{1};
  Jmin = abs(L-S); Jmax = L+S;
  for iJ = Jmin:Jmax
    index = index + 1;
    states{index} = {S statesLS{i}{2:4} iJ i};
  end
end

num_states = length(states); H = sparse(num_states,num_states);
%for z = 1:zcnt; zmat{z} = H; end
for k = 1:6; mk{k} = H; end; for k = 1:4; pk{k} = H; end;

% Equation for the matrix elements is:
%                                      S'+L+J { S' L' J }         (tt)
% <vUSLJ|H|v'U'S'L'J'> = d(J,J') x (-1)       { L  S  t } <vUSL||T    ||v'U'S'L'>
%
% So, to satisfy the six-j symbols we have to have the triads (S'St) and (LL't) non-zero
% i.e. |S-S'|<= t <= S+S' and |L-L'|<= t <= L+L' with t being either 1 or 2.
%
% Ref: Judd, Crosswhite and Crosswhite, Phys. Rev. 169(1968)130

for i = 1:num_states
  S = states{i}{1}; Ls = states{i}{2}; L = racah_lconv(Ls); J = states{i}{5}; irm = states{i}{6};
  for j = 1:num_states
    Sp = states{j}{1}; Lsp = states{j}{2}; Lp = racah_lconv(Lsp); Jp = states{j}{5}; jrm = states{j}{6};
    if J==Jp && (S-Sp)<=2 && (L-Lp)<=2
      for z = 1:3
        %zmat{z}(i,j) = (-1)^(Sp+L+J) * sixj([Sp Lp J; L S 2]) * redmat{z}(irm,jrm);
        mk{z}(i,j) = (-1)^(Sp+L+J) * sixj([Sp Lp J; L S 2]) * redmat{z}(irm,jrm);
      end
    if (S-Sp)<=1 && (L-Lp)<=1
      %for z = 5:12
      for z = 4:6
        %zmat{z}(i,j) = (-1)^(Sp+L+J) * sixj([Sp Lp J; L S 1]) * redmat{z}(irm,jrm);
        mk{z}(i,j) = (-1)^(Sp+L+J) * sixj([Sp Lp J; L S 1]) * redmat{z}(irm,jrm);
      end
      for z = 1:4
        pk{z}(i,j) = (-1)^(Sp+L+J) * sixj([Sp Lp J; L S 1]) * redmat{z+6}(irm,jrm);
      end
    end 
    end % if J==Jp
  end
end

%for z = 1:zcnt
%  H = H + a(z).*zmat{z};
%end
for k = 1:3; H = H + M(k).*(mk{k}+mk{k+3}); end
for k = 1:4; H = H + P(k).*pk{k}; end

if nargin==5
  H = {H mk{:} pk{:}};
end

%if exist('zmatstruct')
%  zmatstruct.(confstring) = zmat;
%  zmat = zmatstruct;
%else
%  zmat.(confstring) = zmat;
%end
if exist('mkstruct'); mkstruct.(confstring) = mk; mk = mkstruct; else; mk.(confstring) = mk; end
if exist('pkstruct'); pkstruct.(confstring) = pk; pk = pkstruct; else; pk.(confstring) = pk; end

save('HssHsoo.mat','mk','pk','-append');

end % if isfield(mk,confstring)
