function Umat = racah_Umat(n,l,k)
% Calculates the U^k matrix elements (UvSL||U^k||U'v'SL') for the f^n configuration

states       = racah_states(n,l);
statesparent = racah_states(n-1,l);

num_states = length(states); Umat = zeros(num_states);

for i = 1:num_states
  S = states{i}{1}; Ls = states{i}{2}; L = racah_lconv(Ls); v  = states{i}{3}; U  = states{i}{4};
  for j = 1:num_states
    Sp = states{j}{1}; Lsp = states{j}{2}; Lp = racah_lconv(Lsp); vp = states{j}{3}; Up = states{j}{4};
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
      Umat(i,j) = n * sumcfp;
    end
  end
end
