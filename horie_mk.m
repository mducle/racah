function [mkmat,Hcell] = horie_mk(n,l,mk)
% Calculates the coefficients of the Marvin integrals for the spin-spin and spin-other-orbit interactions
%
% Syntax:  [mkmat,Hcell] = horie_mk(n,l,mk)
%
% Inputs:  n     - scalar - number of equivalent electrons
%                - cell   - {n l}
%                - string - 'ln' 
%          l     - scalar - orbital angular momentum number of equivalent electrons
%          l     - char   - orbital angular momentum of equivalent electrons as leter, e.g. 's','p','d','f'
%          mk    - vector - [M0 M2 M4] - Marvin's integrals parameters
%
% Outputs: mkmat - matrix - the matrix of coeffients*Mk values - i.e. Hss+Hsoo
%          Hcell - cell   - { {m0 m2 m4} {Hss_m0 Hss_m2 Hss_m4} {Hsoo1_m0 _m2 _m4} {Hsoo2_m0 _m2 _m4} }
%                           cell array of the coefficients of the spin-spin Hss_, spin-other-orbit Hsoo1_
%                           and Hsoo2_ and their sum {m0 m2 m4}

% By Duc Le 2007 - duc.le@ucl.ac.uk

% Parses input arguments
if iscell(n) && length(n)==2
  l = n{2}; n = n{1};
elseif ischar(n) && length(n)==2
  mk = l;
  l = racah_lconv(n(1));
  n = double(n(2))-48;
elseif isnumeric(n) && isscalar(n)
elseif isscalar(l)
  if ischar(l); l = racah_lconv(l); end
else
  error('n must be either lenght 2 cell or string, or a numeric scalar. l must be scalar numeral or char');
end
if l~=3
  error('Sorry only f^n configurations implemented so far...');
end
if n<1 || n>(2*(2*l+1))
  error ('n must be in the interval [1,2(2l+1)]');
end

% Calculates the matrix elements of the V^{1k} and U^k tensors
%for ik = 1:7; V1k{ik} = racah_V1x(n,l,ik-1); Uk{ik} = racah_Umat(n,l,ik-1); end
% Calculates ourselves - faster as we only calculate the cfp's once for all k=0-6 - cuts time in half
states       = racah_states(n,l);
statesparent = racah_states(n-1,l);
s = 1/2; lnst = length(states);
kend = 2*l+1;    % Calcalates from V^{10} to V^{1(kend-1)}

for ik = 1:kend; V1k{ik} = sparse(lnst,lnst); Uk{ik} = sparse(lnst,lnst); end

% Calculates the factors to multiply the lower triangle by: (-1)^((L-S)-(L'-S'))
for i = 1:lnst; 
  Ls{i} = states{i}{2}; L(i) = racah_lconv(Ls{i}); S(i) = states{i}{1}; 
  v(i) = states{i}{3}; U{i} = states{i}{4}; LSdiff(i) = L(i)-S(i); 
end
for i = 1:lnst; lTriFact(:,i) = LSdiff(i)-LSdiff; end; lTriFact = tril((-1).^lTriFact,-1);

for i = 1:lnst
  [coefs,parents,ind_parents]= racah_parents(n,l,v(i),U{i},S(i),Ls{i});
  for j = i:lnst              % Lower triangle is related by: <x|H|x'> = (-1)^((L-S)-(L'-S'))<x'|H|x>
    if abs(S(j)-S(i))>1; continue; end                                    % Triangular conditions on the
    if abs(L(j)-L(i))>(kend-1); continue; end                             % 6j / W symbols.
    [coefp,parentp,ind_parentp]= racah_parents(n,l,v(j),U{j},S(j),Ls{j});
    ind_p = []; 
    for p = 1:length(ind_parents)                                         % Finds index of states with
      ind_p = [ind_p ind_parentp(find(ind_parents(p)==ind_parentp))];     % cfp(psi) and cfp(psi') 
    end                                                                   % both non-zero.
    for ik = 1:kend; sumcfpV(ik) = 0; sumcfpU(ik) = 0; end
    for p = ind_p
      Sparent = statesparent{p}{1}; Lparent = racah_lconv(statesparent{p}{2});
      vparent = statesparent{p}{3}; Uparent = statesparent{p}{4};
      if abs(S(i)-Sparent)>s || abs(S(j)-Sparent)>s; continue; end
      if abs(L(i)-Lparent)>l || abs(L(j)-Lparent)>l; continue; end
      cfp   = coefs(find(p==ind_parents));                                % Looks up the coef. frac.
      cfp_p = coefp(find(p==ind_parentp));                                % parentage for psi and psi'
      WS = racahW([Sparent S(i) s 1 s S(j)]);
      nonkprodV = cfp * cfp_p * WS;                                       % Non k-dependent product
      if S(i)==S(j); nonkprodU = cfp * cfp_p * (-1)^Lparent; end;
      for ik = 1:kend
        WL = racahW([Lparent L(i) l (ik-1) l L(j)]);
	sumcfpV(ik) = sumcfpV(ik) + ( nonkprodV * WL );                   % For <psi||V^(1k)||psi'>
	if S(i)==S(j)
          W = racahW([l L(i) l L(j) Lparent (ik-1)]);
          sumcfpU(ik) = sumcfpU(ik) + ( (-1)^(ik-1) * nonkprodU * W );    % For <psi||U^(k)||psi'>
	end
      end
    end
    nonkprodV = n * sqrt( 3*(2*L(i)+1)*(2*L(j)+1)*(2*S(i)+1)*(2*S(j)+1)/2 );
    if S(i)==S(j); nonkprodU = n * (-1)^(-l-L(j)) * sqrt((2*L(i)+1)*(2*L(j)+1)); end;
    for ik = 1:kend
      V1k{ik}(i,j) = nonkprodV * sumcfpV(ik);
      if S(i)==S(j); Uk{ik}(i,j)  = nonkprodU * sumcfpU(ik); end
    end
  end
end

for ik = 1:kend
  V1k{ik} = V1k{ik} + V1k{ik}'.*lTriFact;
  Uk{ik}  = Uk{ik}  + Uk{ik}'.*lTriFact;
end

%mkmat = V1k; Hcell = Uk; if 0

% Calculates the reduced matrix elements of the tensor <l||C^{k}||l>
for ik = 1:kend
  C(ik) = (-1)^l * sqrt((2*l+1)*(2*l+1)) * threej([l (ik-1) l; 0 0 0]);
end

% Calculates the reduced matrix elements of the double tensor <l||U^{Kk}||l>
for ik = 2:2:kend
  k = ik-1; % Calculates for odd k's only
  U{ik-1}(ik) = (-1)^((k-1)/2) * (2*l+1) * sqrt(prod([k+1 2*k-1 2*k+1])/3) ...
                  * sqrt(prod(factorial([k-1 k+1 2*l-k]))/factorial(2*l+k+1)) ...
                  * factorial(l+(k+1)/2) / prod(factorial([(k-1)/2 (k+1)/2 l-(k+1)/2]));
  U{ik+1}(ik) = -sqrt( k*(2*k+3)/((k+1)*(2*k-1)) ) * U{ik-1}(ik);
end

% We follow the expressions of Horie 1953 (Prog. Theoret. Phys v10, p296)

% Calculates the coefficients for the spin-spin interactions after Horie:
%
%                                          2l
%                      S'+L+J  { S S' 2 }  ---     k-1 ---         1,k-1                1,k+1
% (psi|Hss|psi') = (-1)        { L L' J }  >    Z M    >    (psi||V     ||psi'') (psi||V     ||psi'')
%                                          ---   k     ---
%                                        k(odd)=1      psi''
%                        S+L+S'+L'
%                  * (-1)          { S S' 2   } { L   L'  2   }
%                                  { 1 1  S'' } { k-1 k+1 L'' }
%
%                                            1/2     k-1          k+1
% where Z  = -4 [ 5k(k+1)(2k-1)(2k+1)(2k+3) ]    (l|C   ||l) (l||C   ||l)
%        k
%
for ik = 2:2:kend; k = ik-1; Zk(ik) = (-1)^k * 4*sqrt(prod([5*k k+1 2*k+[-1 1 3]])) * C(ik-1) * C(ik+1); end
for i = 1:lnst
  for j = 1:lnst
    for ik = 2:2:kend
      k = ik-1; sump = 0; phasefact = (-1)^(S(i)+L(i)+S(j)+L(j));
      for pp = find( V1k{ik-1}(i,:).*V1k{ik+1}(:,j)' )
        sump = sump + ( V1k{ik-1}(i,pp)*V1k{ik+1}(pp,j) * phasefact ...
                        * sixj([S(i) S(j) 2; 1 1 S(pp)]) * sixj([L(i) L(j) 2; k+1 k-1 L(pp)]) );
      end
      Hcell{2}{ik/2}(i,j) = sump * Zk(ik);
    end
  end
end

% Calculates the coefficients of the spin-other-orbit interactions Hsoo = HsooI + HsooII:
%
% 
%       I           [     0              3       1/2 ----     K+k      (K)          (Kk)
% (psi|H   |psi') = [ -2nM  +  6 { ------------ }    >    (-1)    (l||C   ||l) (l||U    ||l)
%       soo         [              l(l+1)(2l+1)      ----
%                                                     kK
%                                      ]          11                      1/2
%                    x  { K l l }  k-1 ] x (psi||V  ||psi') [l(l+1)(2l+1)]
%                       { l k 1 } M    ]
%
% Where the sum runs over odd k up to (2l-1), K = k +/- 1, and U^(Kk) is defined as:
%
%      k-1,k           (k-1)/2                              1/2
% (l||U     ||l) = (-1)         (2l+1) [(k+1)(2k-1)(2k+1)/3]
%
%                                         1/2
%                  [ (k-1)!(k+1)!(2l-k)! ]     ((2l+k+1)/2)!  [k+1]  [2l-k+1]
%                x [ ------------------- ]     -------------  [---]! [------]!
%                  [      (2l+k+1)!      ]      ((k-1)/2)!    [ 2 ]  [  2   ]
%
%      k+1,k        [   k(2k+3)   ]1/2      k-1,k
% (l||U     ||l) = -[ ----------- ]    (l||U     ||l)
%                   [ (k+1)(2k-1) ]
%
%Hcell{3}{1} = (-2*n) .* V1k{2} .* sqrt(l*(l+1)*(2*l+1));
%for ik = 4:2:kend
for ik = 2:2:kend
  k = ik-1; sumK = 0;
  for iK = (ik-1):2:(ik+1)
    K = iK-1; sumK = sumK + ( (-1)^(K+k) * C(iK) * U{iK}(ik) * sixj([K l l; l k 1]) );
  end
  Hcell{3}{ik/2} = 6*sqrt(3)*sumK .* V1k{2};
end
% And:
%
%
%       II               1/2 ----     k+K+L+L'      (K)          (K,k)      k-1
% (psi|H   |psi') = -2(3)    >    (-1)         (l||C   ||l) (l||U     ||l) M
%       soo                  ----
%                             kK
%
%                     [  ---        (K)                  (1k)        { L K  L'' }
%                   x [  >   (psi||U   ||psi'') (psi''||V    ||psi') { k L' 1   }
%                     [  ---                                         
%                       psi''
%
%                        ---        (k)                  (1K)        { L k  L'' } ]
%                    + 2 >   (psi||U   ||psi'') (psi''||V    ||psi') { K L' 1   } ]
%                        ---                                                      ]
%                       psi''
for i = 1:lnst
  for j = 1:lnst
    for ik = 2:2:kend
      k = ik-1; sumK = 0; 
      for iK = (ik-1):2:(ik+1)
        K = iK-1; sump1 = 0; sump2 = 0;
	%if k~=1 && K~=0;                 % Exclude the term in k=1, K=0 as included in HsooI
        for pp = find( Uk{iK}(i,:).*V1k{ik}(:,j)' )
          sump1 = sump1 + ( Uk{iK}(i,pp)*V1k{ik}(pp,j) * sixj([L(i) K L(pp); k L(j) 1]) );
        end
	%end
        for pp = find( Uk{ik}(i,:).*V1k{iK}(:,j)' )
          sump2 = sump2 + ( Uk{ik}(i,pp)*V1k{iK}(pp,j) * sixj([L(i) k L(pp); K L(j) 1]) );
        end
	sumK = sumK + ( (-1)^(k+K+L(i)+L(j)) * C(iK) * U{iK}(ik) * (sump1 + 2*sump2) );
      end
      Hcell{4}{ik/2}(i,j) = -2*sqrt(3) * sumK;
    end
  end
end

for ik = 1:3
  Hcell{5}{ik} = Hcell{3}{ik}+ Hcell{4}{ik};
end

mkmat = Hcell;

%index = 0;
%for i = 1:lnst
%  Jmin = abs(L(i)-S(i)); Jmax = L(i)+S(i);
%  for iJ = Jmin:Jmax
%    index = index + 1;
%    statesJ{index} = {S(i) Ls{i} v(i) U{i} iJ i};
%  end
%end

