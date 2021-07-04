function [matrixu,errormatrix] = W2Q1CNneumann(xbeginning,xend,tbeginning,tend,numberofpointsinx,numberofpointsint,eta,zeta1,zeta2,f,actualsol)


%Dirichlet conditions



% for solving equations of the form u_t - u_xx = f(x,t)




matrixu = zeros(numberofpointsint+1, numberofpointsinx+1); % X's across t down (columns)

dx = (xend-xbeginning)/numberofpointsinx;
dt = (tend-tbeginning)/numberofpointsint;
lambda = dt/(dx^2);
gamma = 1/dx;


% Mesh points (uniform)

pointx = [];
pointt = [];
pointx(1) = xbeginning;
pointt(1) = tbeginning;

        for i=1:(numberofpointsinx)
        pointx(i+1) = xbeginning + (i)*dx;
        end
        
        for j=1:(numberofpointsint)
        pointt(j+1) = tbeginning + (j)*dt;
        end

        %Building the matrix of linear coefficients to be inverted on each
        %loop - LHS.
        
        LHSmatrix = zeros(numberofpointsinx+1,numberofpointsinx+1);
        
        LHSmatrix(1,1) = -gamma;
        LHSmatrix(1,2) = gamma;
        
        for i = 2:numberofpointsinx
                     LHSmatrix(i,i-1) = -lambda;
                     LHSmatrix(i,i) = 2*(1+ lambda);
                     LHSmatrix(i,i+1) = -lambda;
        end
        
        LHSmatrix(numberofpointsinx+1,numberofpointsinx) = -gamma;
        LHSmatrix(numberofpointsinx+1,numberofpointsinx+1) = gamma; 
        
        LHSinverted = inv(LHSmatrix);
        
        
        %Building the matrix of linear coefficients to be inverted on each
        %loop - RHS.
        
        RHSmatrix = zeros(numberofpointsinx+1,numberofpointsinx+1);
        for i = 2:numberofpointsinx
                     RHSmatrix(i,i-1) = lambda;
                     RHSmatrix(i,i) = 2*(1- lambda);
                     RHSmatrix(i,i+1) = lambda;
        end
        
       %initial conditions with first boundary point

         for i = 1:numberofpointsinx+1
        matrixu(1,i) = eta(pointx(i));
         end

       %Create loop for each time iteration and also create our values of f for the right hand side. 
        

        for j = 1:numberofpointsint
      
        
             %Setting up values for f RHS.
             
              fvaluesrhs = zeros(numberofpointsinx+1,1);
              
                    for i=2:numberofpointsinx
                        fvaluesrhs(i) = 0.5*(f(pointx(i),pointt(j)) + f(pointx(i),pointt(j+1))); 
                    end
                
            % Set up Neumann boundary points for the first and last point
            
            fvaluesrhs(1) = zeta1(pointx(1),pointt(j) + dt/2);
            fvaluesrhs(numberofpointsinx+1) = zeta2(pointx(numberofpointsinx+1),(pointt(j) + dt/2));
         
         
         
         %setting up rhs known values for each j iteration from previous
         %loop.
         
             rhsknownvalues = zeros(numberofpointsinx+1,1);
         
             for i =1:numberofpointsinx+1
                rhsknownvalues(i) = matrixu(j,i);
             end
            
             
             unknownvalues = LHSinverted*((RHSmatrix*rhsknownvalues) + fvaluesrhs);
         
            for i = 1:numberofpointsinx+1
                matrixu(j+1,i) = unknownvalues(i);
            end
         
        end
        
        
        matrixofcorrectsolutions = zeros(numberofpointsint+1, numberofpointsinx+1);
        for j = 1:numberofpointsint +1
           for i = 1:numberofpointsinx +1
            matrixofcorrectsolutions(j,i) = actualsol(pointx(i),pointt(j));   
           end
        end
   
        
        errormatrix = abs(matrixofcorrectsolutions-matrixu);
        
  
end