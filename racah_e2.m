function e2 = racah_e2(n,l,states)
% Calculates the Racah electrostatic energy operator e2, after Racah IV.

if ~exist('states')
  states = racah_states(n,l);
end

for i = 1:length(states)
  for j = 1:length(states)

    Si = states{i}{1}; Li = states{i}{2}; vi = states{i}{3}; Ui = states{i}{4};
    Sj = states{j}{1}; Lj = states{j}{2}; vj = states{j}{3}; Uj = states{j}{4};
    Wi = racah_vtow(vi,Si,l); Wj = racah_vtow(vj,Sj,l);
    
    if Li==Lj & Si==Sj & vi==vj
      e2(i,j) = racah_e2sign(Si,vi) * sum( racah_xwu(Wi,Ui,Uj) .* racah_chi(Ui,Uj,Li) );
    else
      e2(i,j) = 0;
    end
  end
end

% Adds in the lower triangle to make e2 hermitian
e2diag = eye(size(e2)).*e2;
e2 = e2 + e2' - e2diag;
