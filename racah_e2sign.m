function e2sign = racah_e2sign(S,v)

    if S==0   && v==0; e2sign = 1;
elseif S==1/2 && v==1; e2sign = 1;
elseif S==1   && v==2; e2sign = 1;
elseif S==0   && v==2; e2sign = 1;
elseif S==3/2 && v==3; e2sign = 1;
elseif S==1/2 && v==3; e2sign = 1;
elseif S==1   && v==4; e2sign = 1;
elseif S==0   && v==4; e2sign = 1;
elseif S==1/2 && v==5; e2sign = 1;
elseif S==0   && v==6; e2sign = 1;
elseif S==7/2 && v==7; e2sign = -1;
elseif S==3   && v==6; e2sign = -1;
elseif S==5/2 && v==5; e2sign = -1;
elseif S==5/2 && v==7; e2sign = -1;
elseif S==2   && v==4; e2sign = -1;
elseif S==2   && v==6; e2sign = -1;
elseif S==3/2 && v==5; e2sign = -1;
elseif S==3/2 && v==7; e2sign = -1;
elseif S==1   && v==6; e2sign = -1;
elseif S==1/2 && v==7; e2sign = -1;
else; e2sign = 0; end %error('Invalid combination of S and v - check your code!'); end
