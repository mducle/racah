function theta_k = theta_k(conf,state,statep)
% Calculates the operator equivalent factors theta_k (k=2,4,6) for a f-electron configuration or state
%
% Syntax:  theta_k = theta_k(conf)
%
% Inputs:  conf    - either: a cell array {n l} of the configuration (n=2-7)
%                        or: a string 'fn' (n=2-7) e.g. conf='f5'
%          state   -         a cell array {S L v U J} of the state to calculate for.
%          statep  -         a cell array {S' L' v' U' J} of the state to calculate for.
%
% Outputs: theta_k - vector: a vector [theta_2 theta_4 theta_6] of the operator equivalent factors.

% By Duc Le 2007 - duc.le@ucl.ac.uk

% Tue Mar 13 01:36:21 GMT 2007 - initial release 

theta_k = zeros(1,3);

if iscell(conf) & length(conf)==2
  n = conf{1};
  l = conf{2};
elseif ischar(conf) & length(conf)==2
  n = double(conf(2))-48;
  l = racah_lconv(upper(conf(1)));
else
  error('Input must be a cell array or string representing the f-electron configuration');
end

if ~exist('state')
  LSJ = hundsrule(n,l);
  statesLS = racah_states(n,l);
  for i = 1:length(statesLS)
    if [racah_lconv(statesLS{i}{2}) statesLS{i}{1}]==LSJ(1:2)
      break;
    end;
  end
  state = statesLS{i};
  J = LSJ(3); Jp = J;
elseif iscell(state) & length(state)>4
  J = state{5};
else
  error('Input "state" must be a cell array {S L v U J ...}, length 5 or more');
end

if ~exist('statep')
  statep = state;
  Jp = J;
elseif iscell(statep) & length(statep)>4
  Jp = statep{5};
  if J~=Jp
    return
  end
else
  error('Input "statep" must be a cell array {Sp Lp vp Up J ...}, length 5 or more');
end

statesparent = racah_states(n-1,l);
  
S  = state{1};  Ls  = state{2};  L  = racah_lconv(Ls);  v  = state{3};  U  = state{4};
Sp = statep{1}; Lsp = statep{2}; Lp = racah_lconv(Lsp); vp = statep{3}; Up = statep{4};

for ik = 1:3
  k = ik*2;

  if S==Sp
    [coefp,parentp,ind_parentp]= racah_parents(n,l,vp,Up,Sp,Lsp);
    [coefs,parents,ind_parents]= racah_parents(n,l,v,U,S,Ls);
    ind_p = []; 
    for z = 1:length(ind_parents)                                    % Finds index of states with
      ind_p = [ind_p ind_parentp(find(ind_parents(z)==ind_parentp))];% cfp(psi) and cfp(psi') 
    end                                                              % both non-zero.
    sumcfp = 0;
    for z = ind_p
      Lparent = racah_lconv(statesparent{z}{2});
      cfp   = coefs(find(z==ind_parents));                           % Looks up the coef. frac.
      cfp_p = coefp(find(z==ind_parentp));                           % parentage for psi and psi'
      W = racahW([3 L 3 Lp Lparent k]); 
      sumcfp = sumcfp + ( cfp * cfp_p * (-1)^(Lparent+k-3-Lp) * sqrt((2*L+1)*(2*Lp+1)) * W );
    end
    Umat = n * sumcfp;

    % Calculates the reduced matrix element (psi_J||U^k||psi'_J)
    redmat(ik) = (-1)^(S+k-L-Jp) * sqrt((2*J+1)*(2*Jp+1)) * racahW([L J Lp Jp S k]) * Umat;

  else
    redmat(ik) = 0;
  end

end

if 2*J>2; theta_k(1) =  -8*sqrt( (7*factorial(2*J-2)) /(15*factorial(2*J+3)) ) * redmat(1); end
if 2*J>4; theta_k(2) =  16*sqrt((14*factorial(2*J-4)) /(11*factorial(2*J+5)) ) * redmat(2); end
if 2*J>6; theta_k(3) =-640*sqrt( (7*factorial(2*J-6))/(429*factorial(2*J+7)) ) * redmat(3); end
