function F = Filter2D()

  Globals2D;

  flag = ones(Np,1);
  sk = 1;
  for i=0:N
    for j=0:N
      if(i+j>N)
	flag(sk) = 0;
      end
      sk = sk+1;
    end
  end
  F = V*(diag(flag)/V);
