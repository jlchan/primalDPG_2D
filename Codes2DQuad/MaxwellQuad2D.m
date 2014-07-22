function [HQ,EQ,time] = MaxwellQuad2D(refc, HQ, EQ, FinalTime)

% function [HQ,EQ] = MaxwellQuad2D(refc, HQ, EQ, FinalTime, F)
% Purpose  : Integrate TM-mode Maxwell's until FinalTime starting with
%            initial conditions Hx,Hy,Ez

Globals2D;
time = 0;

% Runge-Kutta residual storage  
resHQ = zeros(Np,K,2); resEQ = zeros(Np,K); 

% compute time step size
hmin = 2*min(min(1./Fscale));
dt = hmin/((N+1)^2);

tstep=1;

% outer time step loop 
while (time<FinalTime)
  
  if(time+dt>FinalTime), dt = FinalTime-time; end

   for INTRK = 1:5    
      % compute right hand side of TM-mode Maxwell's equations
      [rhsHQ, rhsEQ]   =  MaxwellQuadStrongWeakRHS2D(refc,HQ,EQ);
      
      % initiate and increment Runge-Kutta residuals
      resHQ = rk4a(INTRK)*resHQ + dt*(rhsHQ);  
      resEQ = rk4a(INTRK)*resEQ + dt*(rhsEQ); 
        
      % update fields
      HQ = HQ+rk4b(INTRK)*resHQ;  
      EQ = EQ+rk4b(INTRK)*resEQ;        
   end;
   % Increment time
   time = time+dt;

   tstep = tstep+1;

   if(mod(tstep, 400)==0)
     tstep
     
     PlotFieldQuad2D(N, x, y, EQ); axis([-1 1 -1 1 -1 1])
     drawnow; pause(.01);
   end
 end

return
