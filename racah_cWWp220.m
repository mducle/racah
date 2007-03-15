function c = racah_cuup22(W,Wp)
% Looks up the constant of proportionality c(U,U'(22)), after Racah IV

if ~isvector(U) || length(U)~=2 || ~isvector(Up) || length(Up)~=2
  error('U and Up must be vectors of length 2.');
end

t = [0  0  0  0  0  0  0  1  0  0;
     0  0  0  0  0  1  0  0  1  0;
     0  0  1  0  0  0  1  1  1  0;
     0  0  0  1  0  0  1  1  0  1;
     0  0  0  0  1  1  1  1  1  0;
     0  1  0  0  1  2  1  0  2  1;
     0  0  1  1  1  1  3  1  2  1;
     1  0  1  1  1  0  1  2  1  1;
     0  1  1  0  1  2  2  1  3  1;
     0  0  0  1  0  1  1  1  1  1];

      if W==[0 0 0]; i=1;
  elseif W==[1 0 0]; i=2;
  elseif W==[1 1 0]; i=3;
  elseif W==[2 0 0]; i=4;
  elseif W==[1 1 1]; i=5;
  elseif W==[2 1 0]; i=6;
  elseif W==[2 1 1]; i=7;
  elseif W==[2 2 0]; i=8;
  elseif W==[2 2 1]; i=9;
  elseif W==[2 2 2]; i=10;
  else; i = 0;
end
      if Wp==[0 0 0]; j=1;
  elseif Wp==[1 0 0]; j=2;
  elseif Wp==[1 1 0]; j=3;
  elseif Wp==[2 0 0]; j=4;
  elseif Wp==[1 1 1]; j=5;
  elseif Wp==[2 1 0]; j=6;
  elseif Wp==[2 1 1]; j=7;
  elseif Wp==[2 2 0]; j=8;
  elseif Wp==[2 2 1]; j=9;
  elseif Wp==[2 2 2]; j=10;
  else; j = 0;
end

if i~=0 && j~=0
  c = t(i,j);
else
  c = 0;
end
