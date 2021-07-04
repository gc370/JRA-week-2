function [matrixu,errormatrix] = W2Q1CNperiodic(xbeginning,xend,tbeginning,tend,numberofpointsinx,numberofpointsint,eta,f,actualsol)

% for solving equations of the form u_t - u_xx = f(x,t)

matrixu = zeros(numberofpointsint+1, numberofpointsinx+1); % X's across t down (columns)

dx = (xend-xbeginning)/numberofpointsinx;
dt = (tend-tbeginning)/numberofpointsint;
lambda = dt/(dx^2);

i=1;
j=1;

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
        %loop.
        
        LHSmatrix = zeros(numberofpointsinx,numberofpointsinx);
        LHSmatrix(1,1) = 2*(1+lambda);
        LHSmatrix(1,2) = -lambda;
        LHSmatrix(1,numberofpointsinx) = -lambda;
        
        LHSmatrix(numberofpointsinx,1) = -lambda;
        LHSmatrix(numberofpointsinx,numberofpointsinx-1) = -lambda;
        LHSmatrix(numberofpointsinx,numberofpointsinx) = 2*(1+lambda);
        
                for i = 2:numberofpointsinx-1
                     LHSmatrix(i,i-1) = -lambda;
                     LHSmatrix(i,i) = 2*(1+lambda);
                     LHSmatrix(i,i+1) = -lambda;
                end
                
        LHSinverted = inv(LHSmatrix); %invert to move to rhs.
        
        RHSmatrix = zeros(numberofpointsinx,numberofpointsinx);
        RHSmatrix(1,1) = 2*(1-lambda);
        RHSmatrix(1,2) = lambda;
        RHSmatrix(1,numberofpointsinx) = lambda;
        
        RHSmatrix(numberofpointsinx,1) = lambda;
        RHSmatrix(numberofpointsinx,numberofpointsinx-1) = lambda;
        RHSmatrix(numberofpointsinx,numberofpointsinx) = 2*(1-lambda);
        
                for i = 2:numberofpointsinx-1
                     RHSmatrix(i,i-1) = lambda;
                     RHSmatrix(i,i) = 2*(1-lambda);
                     RHSmatrix(i,i+1) = lambda;
                end
        
        % initial conditions and boundary combined for first condition
         for i = 1:numberofpointsinx
        matrixu(1,i) = eta(pointx(i));
        end
       

        for j = 1:numberofpointsint
      
        
             %Setting up values for f RHS.
             
              fvaluesrhs = zeros(numberofpointsinx,1);
                    for i=1:numberofpointsinx-1
                        fvaluesrhs(i) = 0.5*(f(pointx(i),pointt(j)) + f(pointx(i),pointt(j+1))); 
                    end
                
        % Pull known values from matrix u from last iteration into new set
        % of equations
         
             rhsknownvalues = zeros(numberofpointsinx,1);
         
             for i =1:numberofpointsinx
                rhsknownvalues(i) = matrixu(j,i);
            end
         
             unknownvalues = LHSinverted*((RHSmatrix*rhsknownvalues) + fvaluesrhs);
         
            for i = 1:numberofpointsinx
                matrixu(j+1,i) = unknownvalues(i);
            end
         
            %Enforcing boundary conditions
            matrixu(j+1,numberofpointsinx+1) =  matrixu(j+1,1);
         
        end
         
       
        
        
        matrixofcorrectsolutions = zeros(numberofpointsint+1, numberofpointsinx+1);
        for j = 1:numberofpointsint +1
           for i = 1:numberofpointsinx +1
            matrixofcorrectsolutions(j,i) = actualsol(pointx(i),pointt(j));   
           end
        end
   
        
        errormatrix = abs(matrixofcorrectsolutions-matrixu);
        
  
end