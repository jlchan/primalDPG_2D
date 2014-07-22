function V = Vandermonde2D(N, r, s)

  sk = 1;
  V = zeros(length(r), (N+1)*(N+1));
  for i=0:N
    for j=0:N
      V(:,sk) = JacobiP(r, 0, 0, i).*JacobiP(s, 0, 0, j);
      sk = sk+1;
    end
  end
  
