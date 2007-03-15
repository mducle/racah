function V11 = racah_V11(n,l)
% Calculates the V11 matrix for l^n electrons

%statesLS     = racah_states(n,l);
states       = racah_states(n,l);
statesparent = racah_states(n-1,l);

%index = 0;
%for i = 1:length(statesLS)
%  L = racah_lconv(statesLS{i}{2});
%  S = statesLS{i}{1};
%  Jmin = abs(L-S); Jmax = L+S;
%  for iJ = Jmin:Jmax
%    index = index + 1;
%    states{index} = {S L statesLS{i}{3:4} iJ i};
%  end
%end

for i = 1:length(states)
  S = states{i}{1}; Ls = states{i}{2}; L = racah_lconv(Ls); v = states{i}{3}; U = states{i}{4}; 
  for j = 1:length(states)
    Sp = states{j}{1}; Lsp = states{j}{2}; Lp = racah_lconv(Lsp); vp = states{j}{3}; Up = states{j}{4}; 
    [coefp,parentp,ind_parentp]= racah_parents(n,l,vp,Up,Sp,Lsp);
    [coefs,parents,ind_parents]= racah_parents(n,l,v,U,S,Ls);
    ind_p = []; 
    for k = 1:length(ind_parents)                                      % Finds index of states with
      ind_p = [ind_p ind_parentp(find(ind_parents(k)==ind_parentp))];  % cfp(psi) and cfp(psi') 
    end                                                                % both non-zero.
    sumcfp = 0;
    % Using new method with values of cfp(psi) and cfp(psi') both nonzero, calculation time was
    % reduced from 515s to 178s for n=4, l=3 and giving the same results. NB. takes ~26min for n=5.
    for k = ind_p
      Sparent = statesparent{k}{1}; Lparent = racah_lconv(statesparent{k}{2});
      vparent = statesparent{k}{3}; Uparent = statesparent{k}{4};
      cfp   = coefs(find(k==ind_parents));                             % Looks up the coef. frac.
      cfp_p = coefp(find(k==ind_parentp));                             % parentage for psi and psi'
      WS = racahW([Sparent S 1/2 1 1/2 Sp]); 
      WL = racahW([Lparent L l 1 l Lp]);
      sumcfp = sumcfp + ( cfp * cfp_p * WS * WL );
    end
    V11(i,j) = n * sqrt( 3*(2*L+1)*(2*Lp+1)*(2*S+1)*(2*Sp+1)/2 ) * sumcfp;
  end
end
