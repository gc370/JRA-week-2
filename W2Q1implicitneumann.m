function [matrixu,errormatrix] = W2Q1implicitneumann(xbeginning,xend,tbeginning,tend,numberofpointsinx,numberofpointsint,eta,zeta1,zeta2,f,actualsol)


% for solving equations of the form u_t - u_xx = f(x,t)



 

matrixu = zeros(numberofpointsint+1, numberofpointsinx+1); % X's across t down (columns)

dx = (xend-xbeginning)/numberofpointsinx;
dt = (tend-tbeginning)/numberofpointsint;
lambda =  dt/(dx^2);
gamma = (1/dx);

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
        
        matrix2 = zeros(numberofpointsinx+1,numberofpointsinx+1);
        matrix2(1,1) = -gamma;
        matrix2(1,2) = gamma;
        matrix2(numberofpointsinx+1,numberofpointsinx) = -gamma;
        matrix2(numberofpointsinx+1,numberofpointsinx+1) = (gamma);
        
                for i = 2:numberofpointsinx
                     matrix2(i,i-1) = -lambda;
                     matrix2(i,i) = (1 + 2*(lambda));
                     matrix2(i,i+1) = -lambda;
                end
                
                
                
        % Inverting the matrix ready to multiply it onto the known
        % functions.
                
        matrix2 = inv(matrix2);
        
        
         % - Initial conditions for eta (t = 0) for the main matrix of
         % results including boundary conditions
        
        for i = 1:numberofpointsinx+1 
        matrixu(1,i) = eta(pointx(i));
        end

       
        %Build the RHS matrix with known values.
        rhsmatrix = zeros(numberofpointsinx+1,1);
        
        
        
        
        for j=1:numberofpointsint
            
            for i = 2: numberofpointsinx
            rhsmatrix(i) = dt*f(pointx(i),pointt(j+1)) + matrixu(j,i); % dt*f(pointx(i),pointt(j+1)) + matrixu(j,i) for j=1
            end
        
            rhsmatrix(1) = zeta1(pointx(1),pointt(j)); 
            rhsmatrix(numberofpointsinx+1) = zeta2(pointx(numberofpointsinx),pointt(j)); 
            newmatrix = transpose(matrix2*rhsmatrix);
        
           for i = 1: numberofpointsinx+1
            matrixu(j+1,i) = newmatrix(i);
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