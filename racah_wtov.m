function v = racah_wtov(W,S,l)
% Converts Racah's W=(w1,w2,w3) quantum numbers to Racah's seniority quantum numbers

if ~isvector(W) | length(W)~=3
  error('W=(w1,w2,w3) must be a vector of length 3');
elseif sum(W<3)~=3
  error('Components w1,w2,w3 must each be less than 2');
end

a = max(find(W==2));
if ~isempty(a)
  v = 2*(a+S);
else
  b1 = max(find(W==1));  % = a+b
  b0 = min(find(W==0));  % = a+b+1
  b = [b1 b0-1];
  v = 2*l+1 - b(1);
end
