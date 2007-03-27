function out = cfp_orthog_test(conf,LS)

if ischar(conf) & length(conf)==2
  n = double(conf(2))-48;
  l = racah_lconv(lower(conf(1)));
elseif isnumeric(conf) & isvector(conf) & length(conf)==2
  n = conf(1);
  l = conf(2);
  conf = [lower(racah_lconv(l)) char(n+48)];
elseif isnumeric(conf) & isscalar(conf)
  n = conf;
  l = 3
  conf = ['f' char(n+48)];
else
  error('conf must be "fn", [n 3] or [n]');
end

if isnumeric(LS) & length(LS)==2
  LS = [char((2*LS(2)+1)+48) racah_lconv(LS(1))];
elseif ~( ischar(LS) & length(LS)==2 )
  error('LS must be a string: "(2S+1)L" e.g. "6P" ');
end

[states,st_id] = racah_states(n,l);

ist = 0;
for i = 1:length(st_id)
  if strncmp(LS,st_id{i},2)
    ist = ist + 1;
    stcmp{ist} = st_id{i};
    [cfp{ist},par{ist},ipar{ist}] = racah_parents(conf,st_id{i});
  end
end

out = zeros(length(stcmp));

for i = 1:length(stcmp)
  for j = 1:length(stcmp)
    tcfpi = zeros(1,length(st_id)); tcfpj = tcfpi;
    tcfpi(ipar{i}) = cfp{i};
    tcfpj(ipar{j}) = cfp{j};
    out(i,j) = sum(tcfpi .* tcfpj);
  end
end

out(find(abs(out)<1e-5)) = 0;
