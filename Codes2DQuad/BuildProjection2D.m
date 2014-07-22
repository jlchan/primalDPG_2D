function [TToQ, QToT] = BuildProjection2D()

  Globals2D;
  Globals2D;
  
  tol = 1e-13;
  
  Corder = 4*N;
  [cubR,cubS,cubW, Ncub] = Cubature2D(Corder);
  cubW = diag(cubW);

  entsTToQ = 1:Np*Np;
  entsToT = 1:Np*Np;

  TToQ = zeros(K*Np*Nfaces*Np, 3);
  QToT = zeros(K*Np*Nfaces*Np, 3);

  kQ = 0;
  for k1=1:K
    for f1=1:Nfaces
      k2 = EToE(k1,f1);
      f2 = EToF(k1,f1);
      if(k1<k2)
	% found that these two elements map to Qcnt
	kQ = kQ+1;

	% coordinates of kQ vertices
	XQ = VX(EToV(kQ,:)); 
	YQ = VY(EToV(kQ,:));

	% coordinates of k1 vertices
        XT1 = VX(EToV(k1,:));
        YT1 = VY(EToV(k1,:));

	% coordinates of k2 vertices
        XT2 = VX(EToV(k2,:));
        YT2 = VY(EToV(k2,:));

	% build cubature coordinates on sub triangles
	cubXk1 = 0.5*(-XQ(1)*(cubR+cubS) + XQ(2)*(1+cubR) + XQ(4)*(1+cubS));
	cubYk1 = 0.5*(-YQ(1)*(cubR+cubS) + YQ(2)*(1+cubR) + YQ(4)*(1+cubS));
			 		    		    
	cubXk2 = 0.5*(-XQ(3)*(cubR+cubS) + XQ(4)*(1+cubR) + XQ(2)*(1+cubS));
	cubYk2 = 0.5*(-YQ(3)*(cubR+cubS) + YQ(4)*(1+cubR) + YQ(2)*(1+cubS));

	% need to compute (phi_Q, phi_Q)^{-1} (phi_Q, phi_T)
	% integrate by evaluating phi_Q at cubature nodes in triangle k1 and k2
	[r1,s1] = InvertQuadCoords2D(XQ, YQ, cubXk1, cubYk1, tol);
	[r2,s2] = InvertQuadCoords2D(XQ, YQ, cubXk2, cubYk2, tol);

	[rT1, sT1] = InvertTriCoords2D(XT1, YT1, cubXk1, cubYk1);
	[rT2, sT2] = InvertTriCoords2D(XT2, YT2, cubXk2, cubYk2);

	% now build phiQ and phiT on cubature nodes 
	phiQ1 = InterpMatrixuad2D(r1, s1);
	phiQ2 = InterpMatrixuad2D(r2, s2);

	phiT1 = InterpMatrix2D(rT1, sT1);
	phiT2 = InterpMatrix2D(rT2, sT2);

	if(0==1)
	maxrT1 = max(abs(rT1))
	maxsT1 = max(abs(sT1))
	maxrT2 = max(abs(rT2))
	maxsT2 = max(abs(sT2))
	maxr1 = max(abs(r1))
	maxs1 = max(abs(s1))
	maxr2 = max(abs(r2))
	maxs2 = max(abs(s2))
      end

	% compute areas of triangles
	A1 = TriArea2D(XT1, YT1);
	A2 = TriArea2D(XT2, YT2);

	MQ =    transpose(phiQ1)*cubW*(A1/2)*phiQ1;
	MQ = MQ+transpose(phiQ2)*cubW*(A2/2)*phiQ2;

	MT1 = MassMatrix*J(1,k1);
	MT2 = MassMatrix*J(1,k2);

	% L2 projection matrices 
	idsT1 = ((k1-1)*Np+1):k1*Np;
	idsT2 = ((k2-1)*Np+1):k2*Np;
	ids  = ((kQ-1)*Np+1):kQ*Np;

	rows = idsT1'*ones(1,Np);
	cols = ones(Np,1)*ids;
	vals = MT1\(transpose(phiT1)*cubW*phiQ1)*(A1/3);
	QToT(entsToT,:) = [rows(:), cols(:), vals(:)];
	entsToT = entsToT + Np*Np;

	rows = idsT2'*ones(1,Np);
	cols = ones(Np,1)*ids;
	vals = MT2\(transpose(phiT2)*cubW*phiQ2)*(A2/3);
	QToT(entsToT,:) = [rows(:), cols(:), vals(:)];
	entsToT = entsToT + Np*Np;	

	rows = ids'*ones(1,Np);
	cols = ones(Np,1)*idsT1;
	vals = MQ\(transpose(phiQ1)*cubW*phiT1)*(A1/3);
	TToQ(entsTToQ,:) = [rows(:), cols(:), vals(:)];
	entsTToQ = entsTToQ + Np*Np;	

	rows = ids'*ones(1,Np);
	cols = ones(Np,1)*idsT2;
	vals = MQ\(transpose(phiQ2)*cubW*phiT2)*(A2/3);
	TToQ(entsTToQ,:) = [rows(:), cols(:), vals(:)];
	entsTToQ = entsTToQ + Np*Np;	

      end
    end
  end 

  TToQ = myspconvert(TToQ, Np*K,   Np*K, 1e-14);
  QToT = myspconvert(QToT,   Np*K, Np*K, 1e-14);

