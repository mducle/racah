function icstate = ic_showstates(conf,state)

if ~iscell(state) | length(state)<6
  error('Invalid state cell')
end

%for i_st = 1:length(state)
 
if ischar(state{2})
  L = state{2};
else
  L = racah_lconv(state{2});
end

if length(state)==6
  if mod(state{6},1)~=0
    if mod(state{5},1)~=0
      icstate = sprintf('%1g%s%s_%1g/2,m=%1g/2',2*state{1}+1,L,alpha(conf,state),state{5}*2,state{6}*2);
    else
      icstate = sprintf('%1g%s%s_%1g,m=%1g/2',2*state{1}+1,L,alpha(conf,state),state{5},state{6}*2);
    end
  else
    if mod(state{5},1)~=0
      icstate = sprintf('%1g%s%s_%1g/2,m=%1g',2*state{1}+1,L,alpha(conf,state),state{5}*2,state{6});
    else
      icstate = sprintf('%1g%s%s_%1g,m=%1g',2*state{1}+1,L,alpha(conf,state),state{5},state{6});
    end
  end
else
  if mod(state{5},1)~=0
    icstate = sprintf('%1g%s%s_%1g/2',2*state{1}+1,L,alpha(conf,state),state{5}*2);
  else
    icstate = sprintf('%1g%s%s_%1g',2*state{1}+1,L,alpha(conf,state),state{5});
  end
end

function alpha = alpha(conf,state)

S = state{1}; v = state{3}; U = state{4};
if ischar(state{2}); L = state{2}; else; L = racah_lconv(state{2}); end

switch conf
  case 'f5'
    if S==2.5 | (S==1.5 & (L=='S' | L=='L' | L=='M') ) | (S==0.5 & (L=='N' | L=='O') )
      alpha = '';
    elseif S==1.5
      switch L
        case 'P'
	  if U==[1 1]; alpha='1'; else; alpha='2'; end
	case 'D'
	      if v==3; alpha='1';
	  elseif U==[2 0]; alpha='2'; else; alpha='3'; end
	case 'F'
	      if v==3; alpha='1';
	  elseif U==[1 0]; alpha='2'; elseif U==[2 1]; alpha='3'; else; alpha='4'; end
	case 'G'
	      if v==3; alpha='1';
	  elseif U==[2 0]; alpha='2'; elseif U==[2 1]; alpha='3'; else; alpha='4'; end
	case 'H'
	      if U==[1 1]; alpha='1'; elseif U==[2 1]; alpha='2'; else; alpha='3'; end
	case 'I'
	      if v==3; alpha='1';
	  elseif U==[2 0]; alpha='2'; else; alpha='3'; end
        case 'K'
	  if U==[2 1]; alpha='1'; else; alpha='2'; end
	otherwise
	  alpha='';
      end
    else
      switch L
	case 'P'
	      if v==3; alpha='1';
	  elseif U==[1 1]; alpha='2'; elseif U==[3 0]; alpha='3'; else; alpha='4'; end
	case 'D'
	      if v==3; if U==[2 0]; alpha='1'; else; alpha='2'; end
	  elseif v==5; if U==[2 0]; alpha='3'; elseif U==[2 1]; alpha='4'; else; alpha='5'; end; end
	case 'F'
	      if v==1; alpha='1';
	  elseif v==3; alpha='2';
	  elseif v==5;     if U==[1 0]; alpha='3'; elseif U==[2 1]; alpha='4'; 
	               elseif U==[3 0]; alpha='5'; else; alpha='6'; end
	  end
	case 'Fp'
	  alpha='7';
	case 'G'
	      if v==3; if U==[2 0]; alpha='1'; else; alpha='2'; end
	  elseif v==5;     if U==[2 0]; alpha='3'; elseif U==[2 1]; alpha='4'; 
	               elseif U==[3 0]; alpha='5'; else; alpha='6'; end
	  end
	case 'H'
	      if v==3; if U==[1 1]; alpha='1'; else; alpha='2'; end
	  elseif v==5;     if U==[1 1]; alpha='3'; elseif U==[2 1]; alpha='4'; 
	               elseif U==[3 0]; alpha='5'; else; alpha='6'; end
	  end
	case 'Hp'
	  alpha='7';
	case 'I'
	      if v==3; alpha='1';
	  elseif v==5; if U==[2 0]; alpha='2'; elseif U==[3 0]; alpha='3'; else; alpha='4'; end; end
	case 'Ip'
	  alpha='5';
	case 'K'
	      if v==3; alpha='1';
	  elseif v==5; if U==[2 1]; alpha='2'; elseif U==[2 0]; alpha='3'; else; alpha='4'; end; end
	case 'Kp'
	  alpha='5';
	case 'L'
	      if v==3; alpha='1';
	  elseif U==[2 1]; alpha='2'; else; alpha='3'; end
        case 'M'
	  if U==[3 0]; alpha='1'; else; alpha='2'; end
      end	
    end;
  otherwise
    alpha='';
end   % case 'f5'

%f1	2F	1 100 10
%
%f2	3P	2 110 11
%	3F	2 110 10
%	3H	2 110 11
%	1S	0 000 00
%	1D	2 200 20
%	1G	2 200 20
%	1I	2 200 20
%
%f3	4S	3 111 00
%	4D	3 111 20
%	4F	3 111 10
%	4G	3 111 20
%	4I	3 111 20
%	2P	3 210 11
%	2D1	3 210 20
%	2D2	3 210 21
%	2F1	1 100 10
%	2F2	3 210 21
%	2G1	3 210 20
%	2G2	3 210 21
%	2H1	3 210 11
%	2H2	3 210 21
%	2I	3 210 20
%	2K	3 210 21
%	2L	3 210 21
%
%f4	5S	4 111 00
%	5D	4 111 20
%	5F	4 111 10
%	5G	4 111 20
%	5I	4 111 20
%	3P1	2 110 11
%	3P2	4 211 11
%	3P3	4 211 30
%	3D1	4 211 20
%	3D2	4 211 21
%	3F1	2 110 10
%	3F2	4 211 10
%	3F3	4 211 21
%	3F4	4 211 30
%	3G1	4 211 20
%	3G2	4 211 21
%	3G3	4 211 30
%	3H1	2 110 11
%	3H2	4 211 11
%	3H3	4 211 21
%	3H4	4 211 30
%	3I1	4 211 20
%	3I2	4 211 30
%	3K1	4 211 21
%	3K2	4 211 30
%	3L	4 211 21
%	3M	4 211 30
%	1S1	0 000 00
%	1S2	4 220 22
%	1D1	2 200 20
%	1D2	4 220 20
%	1D3	4 220 21
%	1D4	4 220 22
%	1F	4 220 21
%	1G1	2 200 20
%	1G2	4 220 20
%	1G3	4 220 21
%	1G4	4 220 22
%	1H1	4 220 21
%	1H2	4 220 22
%	1I1	2 200 20
%	1I2	4 220 20
%	1I3	4 220 22
%	1K	4 220 21
%	1L1	4 220 21
%	1L2	4 220 22
%	1N	4 220 22
%
%
%f6	7F	6 100 10
%	5S	4 111 00
%	5P	6 210 11
%	5D1	4 111 20
%	5D2	6 210 20
%	5D3	6 210 21
%	5F1	4 111 10
%	5F2	6 210 21
%	5G1	4 111 20
%	5G2	6 210 20
%	5G3	6 210 21
%	5H1	6 210 11
%	5H2	6 210 21
%	5I1	4 111 20
%	5I2	6 210 20
%	5K	6 210 21
%	5L	6 210 21
%	3P1	2 110 11
%	3P2	4 211 11
%	3P3	4 211 30
%	3P4	6 221 11
%	3P5	6 221 30
%	3P6	6 221 31
%	3D1	4 211 20
%	3D2	4 211 21
%	3D3	6 221 20
%	3D4	6 221 21
%	3D5	6 221 31
%	3F1	2 110 10
%	3F1	4 211 10
%
%	3F3	4 211 21
%	3F4	4 211 30
%	3F5	6 221 10
%	3F6	6 221 21
%	3F7	6 221 30
%	3F8	6 221 31A
%	3F9	6 221 31B
%	3G1	4 211 20
%	3G2	4 211 21
%	3G3	4 211 30
%	3G4	6 221 20
%	3G5	6 221 21
%	3G6	6 221 30
%	3G7	6 221 31
%	3H1	2 110 11
%	3H2	4 211 11
%	3H3	4 211 21
%	3H4	4 211 30
%	3H5	6 221 11
%	3H6	6 221 21
%	3H7	6 221 30
%	3H8	6 221 31A
%	3H9	6 221 31B
%	3I1	4 211 20
%	3I2	4 211 30
%	3I3	6 221 20
%	3I4	6 221 30
%	3I5	6 221 31A
%	3I6	6 221 31B
%	3K1	4 211 21
%
%	3K2	4 211 30
%	3K3	6 221 21
%	3K4	6 221 30
%	3K5	6 221 31A
%	3K6	6 221 31B
%	3L1	4 211 21
%	3L2	6 221 21
%	3L3	6 221 31
%	3M1	4 211 30
%	3M2	6 221 30
%	3M3	6 221 31
%	3N	6 221 31
%	3O	6 221 31
%	1S1	0 000 00
%	1S2	4 220 00
%	1S3	6 222 40
%	1S4	6 222 30
%	1P	6 222 30
%	1D1	2 200 20
%	1D2	4 220 20
%	1D3	4 220 21
%	1D4	4 220 22
%	1D5	6 222 20
%	1D6	6 222 40
%	1F1	4 220 21
%	1F2	6 222 10
%	1F3	6 222 30
%	1F4	6 222 40
%	1G1	2 200 20
%	1G2	4 220 20
%
%	1G3	4 220 21
%	1G4	4 220 22
%	1G5	6 222 20
%	1G6	6 222 30
%	1G7	6 222 40A
%	1G8	6 222 40B
%	1H1	4 220 21
%	1H2	4 220 22
%	1H3	6 222 30
%	1H4	6 222 40
%	1I1	2 200 20
%	1I2	4 220 20
%	1I3	4 220 22
%	1I4	6 222 20
%	1I5	6 222 30
%	1I6	6 222 40A
%	1I7	6 222 40B
%	1K1	4 220 21
%	1K2	6 222 30
%	1K3	6 222 40
%	1L1	4 220 21
%	1L2	4 220 22
%	1L3	6 222 40A
%	1L4	6 222 40B
%	1M1	6 222 30
%	1M2	6 222 40
%	1N1	4 220 22
%	1N2	6 222 40
%	1Q	6 222 40
%
%7f	8S	7 000 00
%	6P	5 110 11
%	6D	7 200 20
%	6F	5 110 10
%	6G	7 200 20
%	6H	5 110 10
%	6I	7 200 20
%	4S1	3 111 00
%	4S2	7 220 22
%	4P1	5 211 11
%	4P2	5 211 30
%	4D1	3 111 20
%	4D2	5 211 20
%	4D3	5 211 21
%	4D4	7 220 20
%	4D5	7 220 21
%	4D6	7 220 22
%	4F1	3 111 10
%	4F2	5 211 10
%	4F3	5 211 21
%	4F4	5 211 30
%	4F5	7 220 21
%	4G1	3 111 20
%	4G2	5 211 20
%	4G3	5 211 21
%	4G4	5 211 30
%	4G5	7 220 20
%	4G6	7 220 21
%	4G7	7 220 22
%	4H1	5 211 11
%
%	4H2	5 211 21
%	4H3	5 211 30
%	4H4	7 220 21
%	4H5	7 220 22
%	4I1	3 111 20
%	4I2	5 211 20
%	4I3	5 211 30
%	4I4	7 220 20
%	4I5	7 220 22
%	4K1	5 211 21
%	4K2	5 211 30
%	4K3	7 220 21
%	4L1	5 211 21
%	4L2	7 220 21
%	4L3	7 220 22
%	4M	5 211 30
%	4N	7 220 22
%	2S1	7 222 00
%	2S2	7 222 40
%	2P1	3 210 11
%	2P2	5 221 11
%	2P3	5 221 30
%	2P4	5 221 31
%	2P5	7 222 30
%	2D1	3 210 20
%	2D2	3 210 21
%	2D3	5 221 20
%	2D4	5 221 21
%	2D5	5 221 31
%	2D6	7 222 20
%
%	2D7	7 222 40
%	2F1	1 100 10
%	2F2	3 210 21
%	2F3	5 221 10
%	2F4	5 221 21
%	2F5	5 221 30
%	2F6	5 221 31A
%	2F7	5 221 31B
%	2F8	7 222 10
%	2F9	7 222 30
%	2F10	7 222 40
%	2G1	3 210 20
%	2G2	3 210 21
%	2G3	5 221 20
%	2G4	5 221 21
%	2G5	5 221 30
%	2G6	5 221 31
%	2G7	7 222 20
%	2G8	7 222 30
%	2G9	7 222 40A
%	2G10	7 222 40B
%	2H1	3 210 11
%	2H2	3 210 21
%	2H3	5 221 11
%	2H4	5 221 21
%	2H5	5 221 30
%	2H6	5 221 31A
%	2H7	7 222 31B
%	2H8	7 222 30
%	2H9	7 222 40
%
%	2I1	3 210 20
%	2I2	5 221 20
%	2I3	5 221 30
%	2I4	5 221 31A
%	2I5	5 221 31B
%	2I6	7 222 20
%	2I7	7 222 30
%	2I8	7 222 40A
%	2I9	7 222 40B
%	2K1	3 210 21
%	2K2	5 221 21
%	2K3	5 221 30
%	2K4	5 221 31A
%	2K5	5 221 31B
%	2K6	7 222 30
%	2K7	7 222 40
%	2L1	3 210 21
%	2L2	5 221 21
%	2L3	5 221 31
%	2L4	7 222 40A
%	2L5	7 222 40B
%	2M1	5 221 30
%	2M2	5 221 31
%	2M3	7 222 30
%	2M4	7 222 40
%	2N1	5 221 31
%	2N2	7 222 40
%	2O	5 221 31
%	2Q	7 222 40
