function [cfp,parents,ind_parent] = racah_parents(n,l,v,U,S,L)
% Finds the parents of an l^n configuration.

if nargin<4
  if nargin==3 && ischar(l) && ischar(v)
    conf = [l n+48]
    state = v;
  else
    conf = n;
    state = l;
    n = double(conf(2))-48;
    l = conf(1);
  end
  stc = racah_states(conf,state);  % Child state
  U  = stc{4}; v  = stc{3}; S  = stc{1}; L  = stc{2};
end

if ~ischar(L)
  errmsg = sprintf(['You must specify the orbital angular momentum quantum number as S,P,D,F,G, etc. \n' ...
         'This is to specify certain combinations where the |vUSL> basis is not enough to distinguish \n' ...
	 'between states, eg. 2F6 and 2F7, specified by S=1/2 v=5 U=(31) and L = "F" or "Fp"'],'');
  error(errmsg);
end

Ls = L; L = racah_lconv(L);

states_parent = racah_states(n-1,l);

index = 0; parents = {}; cfp = []; ind_parent = [];
for i = 1:length(states_parent)
  Sp = states_parent{i}{1};
  if abs(Sp-S)<=1/2             % Difference between S_child and S_parent <= S_single_electron=1/2
    Lsp = states_parent{i}{2}; Lp = racah_lconv(Lsp);
    if abs(Lp-L)<=l             % Difference between L_child and L_parent <= l_single_electron=l
      cfp_tmp = racah_cfp(n,U,v,S,Ls,states_parent{i}{4},states_parent{i}{3},Sp,Lsp);
      if cfp_tmp~=0
        index = index + 1;
        parents{index} = states_parent{i};
        cfp(index) = cfp_tmp;
        ind_parent(index) = i;
      end
    end
  end
end

%if nargout==1
%  cfp = {cfp parents};
%end
