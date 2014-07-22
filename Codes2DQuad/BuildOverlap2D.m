function [F] = BuildOverlap2D()

  Globals2D;
  Globals2D;
  
  N = N;
  Np = (N+1)*(N+1);
  r = JacobiGL(0, 0, N)*ones(1,N+1);
  s = r';
  r = r(:); s = s(:);
  V = Vandermonde2D(N, r, s);
  
  tol = 1e-13;

  vnums = [ 1 2; 2 3 ; 3 1];
  opnum = [ 3; 1; 2];
  anums = [ 1 2 3; 2 3 1; 3 1 2];

  Corder = 6*N;
  [cubR,cubS,cubW, Ncub] = Cubature2D(Corder);
  cubW = diag(cubW);

  cx = sum(VX(EToV),2)/Nfaces;
  cy = sum(VY(EToV),2)/Nfaces;

  ents = 1:Np*Np;

  F = zeros(K*Np*Nfaces*Np, 3);

  for k1=1:K
    for f1=1:Nfaces
      k2 = EToE(k1,f1);
      f2 = EToF(k1,f1);
      if(k1<k2)

	% now build projection operator
        XT1 = VX(EToV(k1,:));
        YT1 = VY(EToV(k1,:));

	% coordinates of k2 vertices
        XT2 = VX(EToV(k2,:));
        YT2 = VY(EToV(k2,:));
	
	% build quad verts in place
	XQ = [XT1(opnum(f1)),XT1(vnums(f1,1)),XT2(opnum(f2)),XT1(vnums(f1,2))];
	YQ = [YT1(opnum(f1)),YT1(vnums(f1,1)),YT2(opnum(f2)),YT1(vnums(f1,2))];


	% build cubature coordinates on sub triangles
	cubXk1 = 0.5*(-XQ(1)*(cubR+cubS) + XQ(2)*(1+cubR) + XQ(4)*(1+cubS));
	cubYk1 = 0.5*(-YQ(1)*(cubR+cubS) + YQ(2)*(1+cubR) + YQ(4)*(1+cubS));
			 		    		    
	cubXk2 = 0.5*(-XQ(3)*(cubR+cubS) + XQ(4)*(1+cubR) + XQ(2)*(1+cubS));
	cubYk2 = 0.5*(-YQ(3)*(cubR+cubS) + YQ(4)*(1+cubR) + YQ(2)*(1+cubS));
	
	cubXk1sub = 0.5*(-cx(k1)*(cubR+cubS) + XQ(2)*(1+cubR) + XQ(4)*(1+cubS));
	cubYk1sub = 0.5*(-cy(k1)*(cubR+cubS) + YQ(2)*(1+cubR) + YQ(4)*(1+cubS));

	XT1sub = [ cx(k1), XT1(vnums(f1,:))];
	YT1sub = [ cy(k1), YT1(vnums(f1,:))];

	A1sub = TriArea2D(XT1sub,YT1sub);

	% need to compute (phi_Q, phi_Q)^{-1} (phi_Q, phi_T)
	% integrate by evaluating phi_Q at cubature nodes in triangle k1 and k2
	[r1,s1] = InvertQuadCoords2D(XQ, YQ, cubXk1, cubYk1, tol);
	[r2,s2] = InvertQuadCoords2D(XQ, YQ, cubXk2, cubYk2, tol);
	[r1sub,s1sub] = InvertQuadCoords2D(XQ,  YQ,  cubXk1sub, cubYk1sub, tol);

	[rT1, sT1] = InvertTriCoords2D(XT1, YT1, cubXk1, cubYk1);
	[rT2, sT2] = InvertTriCoords2D(XT2, YT2, cubXk2, cubYk2);
	[rT1sub,sT1sub] = InvertTriCoords2D (XT1, YT1, cubXk1sub, cubYk1sub);

	if(k1<10)
	subplot(1,2,1);
	scatter(cubXk1, cubYk1, 'ro');
	hold on
	scatter(cubXk2, cubYk2, 'k*');
	scatter(cubXk1sub, cubYk1sub, 'bs');
	hold off;
	
	subplot(1,2,2);
	scatter(r1sub,s1sub);
	hold on;
	scatter(r1,s1, 'r*')
	hold off
	drawnow;
	pause(.02);
	end

	% now build phiQ and phiT on cubature nodes 
	phiQ1 = InterpMatrixuad2D(r1, s1);
	phiQ2 = InterpMatrixuad2D(r2, s2);
	phiQ1sub = InterpMatrixuad2D(r1sub, s1sub);

	phiT1    = InterpMatrix2D(rT1, sT1);
	phiT2    = InterpMatrix2D(rT2, sT2);
	phiT1sub = InterpMatrix2D(rT1sub, sT1sub);

%	[min(r1sub),max(r1sub)]
%	[min(s1sub),max(s1sub)]

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

        subP = MassMatrix\(transpose(phiT1sub)*cubW*phiQ1sub*(A1sub/A1));

	% L2 projection matrices 
	idsT1 = ((k1-1)*Np+1):k1*Np;
	idsT2 = ((k2-1)*Np+1):k2*Np;

	rows = idsT1'*ones(1,Np);
	cols = ones(Np,1)*idsT1;
	vals = subP*(MQ\(transpose(phiQ1)*cubW*phiT1)*(A1/2));
	F(ents,:) = [rows(:), cols(:), vals(:)];
	ents = ents + Np*Np;	

	rows = idsT1'*ones(1,Np);
	cols = ones(Np,1)*idsT2;
	vals = subP*(MQ\(transpose(phiQ2)*cubW*phiT2)*(A2/2));
	F(ents,:) = [rows(:), cols(:), vals(:)];
	ents = ents + Np*Np;	

      end
    end
  end 

  F = myspconvert(F, Np*K,   Np*K, 1e-14);

