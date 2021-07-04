function [matrixu] = W2Q1CN(xbeginning,xend,tbeginning,tend,numberofpointsinx,numberofpointsint,eta,gx1,gxend,f,actualsol)


%Dirichlet conditions



% for solving equations of the form u_t - u_xx = f(x,t)


%inputs

%xbeginning =0 ;
%xend = 1;
%tbeginning =0 ;
%tend = 1;
%numberofpointsinx = 10;
%numberofpointsint= 5;
%eta = @(x) sin(pi*x) ; %.... ;
%alpha1 =0;
%beta1 =1;
%alpha2 =0;
%beta2 =1;
%zeta1 = @(x,t) 0;
%zeta2 = @(x,t) 0;
%f = @(x,t) 0;
%actualsol = @(x,t) exp(-(pi^2)*t)*sin(pi*x);
%gx1 = 0;
%gxend =0;
 

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
        
        A = zeros(numberofpointsinx-1,numberofpointsinx-1);
        I = eye(numberofpointsinx-1); % identity matrix
     
       
        A(1,1) = -2*lambda;
        A(1,2) = lambda;
        A(numberofpointsinx-1,numberofpointsinx-2) = lambda;
        A(numberofpointsinx-1,numberofpointsinx-1) = -2*lambda;
        
                for i = 2:numberofpointsinx-2
                     A(i,i-1) = lambda;
                     A(i,i) = -2*lambda;
                     A(i,i+1) = lambda;
                end
                
        LHSmatrix = (2*I) - A;
        LHSinverted = inv(LHSmatrix); %invert to move to rhs.
        RHSmatrix = (2*I) + A;
        
        
        % Create matrixu to contain our results and input initial
        % conditions and boundary conditions using ghost point
        
        
         for i = 2:numberofpointsinx
        matrixu(1,i) = eta(pointx(i));
        end
        % Starting Boundary conditions 
        matrixu(1,1) = gx1;
        matrixu(1,numberofpointsinx+1) = gxend;
        

        for j = 1:numberofpointsint
      
        
             %Setting up values for f RHS.
             
              fvaluesrhs = zeros(numberofpointsinx-1,1);
                    for i=1:numberofpointsinx-1
                        fvaluesrhs(i) = 0.5*(f(pointx(i),pointt(j)) + f(pointx(i),pointt(j+1))); 
                    end
                
         %include boundary values into a new vector and then add it to fvaluesrhs
         
             boundaryvaluesrhs = zeros(numberofpointsinx-1,1);
             boundaryvaluesrhs(1) = 2*lambda*gx1;                       %two since we have the boundary at the j+1 and j'th step known
             boundaryvaluesrhs(numberofpointsinx-1)=2*lambda*gxend;     %two since we have the boundary at the j+1 and j'th step known
         
             fvaluesandboundaries = boundaryvaluesrhs + fvaluesrhs;     % plus since they are added from the LHS
         
         
         %setting up rhs known values for each j iteration from previous
         %loop.
         
             rhsknownvalues = zeros(numberofpointsinx-1,1);
         
             for i =1:numberofpointsinx-1
                rhsknownvalues(i) = matrixu(j,i+1);
            end
         
             unknownvalues = LHSinverted*((RHSmatrix*rhsknownvalues) + fvaluesandboundaries);
         
            for i = 2:numberofpointsinx
                matrixu(j+1,i) = unknownvalues(i-1);
            end
         
            matrixu(j+1,1) = gx1;
            matrixu(j+1,numberofpointsinx+1) = gxend;
         
        end
         
       
        
        
        matrixofcorrectsolutions = zeros(numberofpointsint+1, numberofpointsinx+1);
        for j = 1:numberofpointsint +1
           for i = 1:numberofpointsinx +1
            matrixofcorrectsolutions(j,i) = actualsol(pointx(i),pointt(j));   
           end
        end
   
        
        errormatrix = abs(matrixofcorrectsolutions-matrixu);
        
  
end