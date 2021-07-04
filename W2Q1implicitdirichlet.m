function [matrixu,errormatrix] = W2Q1implicitdirichlet(xbeginning,xend,tbeginning,tend,numberofpointsinx,numberofpointsint,eta,gx1,gxend,f,actualsol)


% for solving equations of the form u_t - u_xx = f(x,t)


matrixu = zeros(numberofpointsint+1, numberofpointsinx+1); % X's across t down (columns)

dx = (xend-xbeginning)/numberofpointsinx;
dt = (tend-tbeginning)/numberofpointsint;
lambda =  dt/(dx^2);

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
        
        matrix2 = zeros(numberofpointsinx-1,numberofpointsinx-1);
        matrix2(1,1) = (1 + 2*(lambda));
        matrix2(1,2) = -lambda;
        matrix2(numberofpointsinx-1,numberofpointsinx-2) = -lambda;
        matrix2(numberofpointsinx-1,numberofpointsinx-1) = (1 + 2*(lambda));
        
                for i = 2:numberofpointsinx-2
                     matrix2(i,i-1) = -lambda;
                     matrix2(i,i) = (1 + 2*(lambda));
                     matrix2(i,i+1) = -lambda;
                end
   
        matrix2 = inv(matrix2);
        
         % - Initial conditions for eta (t = 0) for the main matrix of
         % results
        
        for i = 1:numberofpointsinx+1
        matrixu(1,i) = eta(pointx(i));
        end

        %Build the RHS matrix with known values.
        rhsmatrix = zeros(numberofpointsinx-1,1);
  
        for j=1:numberofpointsint
            
            
        for i = 2: numberofpointsinx-2
        rhsmatrix(i) = dt*f(pointx(i),pointt(j+1)) + matrixu(j,i+1); % dt*f(pointx(i),pointt(j+1)) + matrixu(j,i) for j=1
        end
        
        rhsmatrix(1) = matrixu(j,2) + dt*f(pointx(1),pointt(j+1)) + lambda*gx1;
        rhsmatrix(numberofpointsinx-1) = matrixu(j,numberofpointsinx-1) + dt*f(pointx(numberofpointsinx),pointt(j+1)) + lambda*gxend;
        newmatrix = transpose(matrix2*rhsmatrix);
        
        for i = 2: numberofpointsinx
        matrixu(j+1,i) = newmatrix(i-1);
        end
        
        matrixu(j+1,1) = gx1;
        matrixu(j+1,numberofpointsinx+1) =gxend;
        
        end
        
        
        
        matrixofcorrectsolutions = zeros(numberofpointsint+1, numberofpointsinx+1);
        for j = 1:numberofpointsint +1
           for i = 1:numberofpointsinx +1
            matrixofcorrectsolutions(j,i) = actualsol(pointx(i),pointt(j));   
           end
        end
   
        
        errormatrix = abs(matrixofcorrectsolutions-matrixu);
        
  
end