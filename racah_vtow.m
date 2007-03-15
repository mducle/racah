function W = racah_vtow(v,S,l)
% Converts Racah's W=(w1,w2,w3) quantum numbers to Racah's seniority quantum numbers

if ~isscalar(v) 
  error('The seniority number v must be a scalar');
end

W = zeros(1,3); 

a = v/2 - S;
b = min([2*S 2*l+1-v]);

if a<0 
  W = [0 0 0];
else
  W(1:a) = 2;
  if (b>0)
    W((a+1):(a+b)) = 1;
  end
  if (a+b)<3
    W((a+b+1):3) = 0;
  end
end
