function Ln = racah_lconv(L)

if isnumeric(L)
  switch L
    case 0;  Ln = 'S';
    case 1;  Ln = 'P';
    case 2;  Ln = 'D';
    case 3;  Ln = 'F';
    case 4;  Ln = 'G';
    case 5;  Ln = 'H';
    case 6;  Ln = 'I';
    case 7;  Ln = 'K';
    case 8;  Ln = 'L';
    case 9;  Ln = 'M';
    case 10; Ln = 'N';
    case 11; Ln = 'O';
    case 12; Ln = 'Q';
  otherwise
    error('Angular momentum quatumn number must be L<=12 (term symbol Q)');
  end
elseif ischar(L)
  switch upper(L(1))
    case 'S'; Ln = 0;
    case 'P'; Ln = 1;
    case 'D'; Ln = 2;
    case 'F'; Ln = 3;
    case 'G'; Ln = 4;
    case 'H'; Ln = 5;
    case 'I'; Ln = 6;
    case 'K'; Ln = 7;
    case 'L'; Ln = 8;
    case 'M'; Ln = 9;
    case 'N'; Ln = 10;
    case 'O'; Ln = 11;
    case 'Q'; Ln = 12;
  otherwise
    error('If using string for total angular momentum number must be letters SPDF G-Q');
  end
end
