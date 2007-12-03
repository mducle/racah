function X = rajnak_Xeq9(n,l,k,kp,lp)

states       = racah_states(n,l); 
statesparent = racah_states(n-1,l);

num_states = length(states); X = zeros(num_states);

for i = 1:num_states
  L(i) = racah_lconv(states{i}{2});
end 

for i = 1:num_states    
  S = states{i}{1}; Ls = states{i}{2}; v  = states{i}{3}; U  = states{i}{4};
  [coefs,parents,ind_parents]= racah_parents(n,l,v,U,S,Ls);
  cfp = zeros(1,17); cfp(ind_parents) = coefs;
  for j = 1:num_states
    Sp = states{j}{1}; Lsp = states{j}{2}; vp = states{j}{3}; Up = states{j}{4};
    [coefp,parentp,ind_parentp]= racah_parents(n,l,vp,Up,Sp,Lsp);
    cfpp = zeros(1,17); cfpp(ind_parentp) = coefp;
    for kpp = 2:2:6
      sumbt = 0;
      for b = 1:num_states             % Sums over psi_bar
        for t = 1:num_states           % Sums over psi_tilde
          sumbt = sumbt + ( (-1)^(L(i)+l+L(t)) * sqrt((2*L(b)+1)*(2*L(t)+1)) * cfp(b)*cfpp(t) ...
	                    * sixj([L(b) L(t) kp; l l L(i)]) * ninej([l l k; l l kpp; L(b) L(t) kp]) );
        end
      end
      X(i,j) = X(i,j) + ( (2*kpp+1) * sixj([k kp kpp; l l lp]) * sumbt );
    end
  end
end
