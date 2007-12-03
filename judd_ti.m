function ti = judd_ti(n,l)
% Calculates the three-particle configuration interaction operators t_i, after Judd 1966

K  = {[2 2 2] [2 2 4] [2 4 4] [2 4 6] [4 4 4] [4 4 6] [2 6 6] [4 6 6] [6 6 6]};
WU = {[0 0 0] [0 0];
      [2 2 0] [2 2];
      [2 2 2] [0 0];
      [2 2 2] [4 0];
      [4 0 0] [4 0];
      [4 2 0] [2 2];
      [4 2 0] [4 0];
      [4 2 0] [4 2];
      [6 0 0] [6 0]};

for i = 1:length(K)
  Vkkk{i} = judd_Vkkk(n,l,K{i});
end

for i = 1:length(WU)
  ti{i} = zeros(size(Vkkk{1}));
  tcoef = judd_tcoef(WU{i,1},WU{i,2});
  for j = 1:length(K)
    ti{i} = ti{i} + tcoef(j).*Vkkk{j};
  end
end
