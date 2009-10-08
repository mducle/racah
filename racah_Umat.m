function Umat = racah_Umat(n,l,k,states)
% Calculates the U^k matrix elements (UvSL||U^k||U'v'SL') for the f^n configuration

nn = n; if nn>7; n=14-nn; end;

if ~exist('states')
states       = racah_states(n,l); end
statesparent = racah_states(n-1,l);

num_states = length(states); Umat = sparse(num_states,num_states);

% Calculates the factors to multiply the lower triangle by: (-1)^((L-S)-(L'-S'))
for i = 1:num_states;
  Ls{i} = states{i}{2}; L(i) = racah_lconv(Ls{i}); S(i) = states{i}{1};
  v(i) = states{i}{3}; U{i} = states{i}{4}; LSdiff(i) = L(i)-S(i);
end
for i = 1:num_states; lTriFact(:,i) = LSdiff(i)-LSdiff; end; lTriFact = tril((-1).^lTriFact,-1);

for i = 1:num_states
  [coefs,parents,ind_parents]= racah_parents(n,l,v(i),U{i},S(i),Ls{i});
  for j = i:num_states    
    if abs(L(i)-L(j))>k; continue; end                                 % Condition on 6j / W symbol
    if S(i)==S(j)
      [coefp,parentp,ind_parentp]= racah_parents(n,l,v(j),U{j},S(j),Ls{j});
      ind_p = []; 
      for z = 1:length(ind_parents)                                    % Finds index of states with
        ind_p = [ind_p ind_parentp(find(ind_parents(z)==ind_parentp))];% cfp(psi) and cfp(psi') 
      end                                                              % both non-zero.
      sumcfp = 0;
      nonzprod = (-1)^(-l-L(j)) * sqrt((2*L(i)+1)*(2*L(j)+1));
      for z = ind_p
        Lparent = racah_lconv(statesparent{z}{2});
        cfp   = coefs(find(z==ind_parents));                           % Looks up the coef. frac.
        cfp_p = coefp(find(z==ind_parentp));                           % parentage for psi and psi'
        W = racahW([l L(i) l L(j) Lparent k]); 
	sumcfp = sumcfp + ( cfp * cfp_p * (-1)^(Lparent+k) * nonzprod * W );
      end
      Umat(i,j) = n * sumcfp;
    end
  end
end

Umat = Umat + Umat'.*lTriFact;

if nn>7; Umat = -(-1)^k .* Umat; end   % Cf. Racah 3.
