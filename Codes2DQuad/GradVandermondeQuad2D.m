function [Vr, Vs] = GradVandermonde2D(N, r, s)

  sk = 1;
  Vr = zeros(length(r), (N+1)*(N+1));
  Vs = zeros(length(r), (N+1)*(N+1));
  for i=0:N
    for j=0:N
      Vr(:,sk) = GradJacobiP(r, 0, 0, i).*JacobiP(s, 0, 0, j);
      Vs(:,sk) = JacobiP(r, 0, 0, i).*GradJacobiP(s, 0, 0, j);
      sk = sk+1;
    end
  end
  
