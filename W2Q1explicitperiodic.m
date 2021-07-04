function [matrixu,errorsmatrix,matrixofcorrectsolutions] = W2Q1explicitperiodic(xbeginning,xend,tbeginning,tend,numberofpointsinx,numberofpointsint,eta,f,actualsol)


% for solving equations of the form u_t - u_xx = f(x,t)

matrixu = zeros(numberofpointsint+1, numberofpointsinx+1); % X's across t down (columns)

dx = (xend-xbeginning)/numberofpointsinx;
dt = (tend-tbeginning)/numberofpointsint;
lambda = dt/(dx^2);


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
        
        % - Initial conditions for eta (t = 0)
        
        for i = 1:numberofpointsinx
        matrixu(1,i) = eta(pointx(i));
        end
        matrixu(1,numberofpointsinx+1) = matrixu(1,1); % Periodic condition boundary

        % Build the matrix of linear coefficients
        
        RHSmatrix = zeros(numberofpointsinx,numberofpointsinx); %since we are reducing dimensionality for bc.
        
        for i=2: numberofpointsinx-1
            RHSmatrix(i,i-1) = lambda;
            RHSmatrix(i,i) = 1 - 2*lambda;
            RHSmatrix(i,i+1) = lambda;
        end
        
        %building the first row
        
        RHSmatrix(1,1) = 1 - 2*lambda;
        RHSmatrix(1,2) = lambda;
        RHSmatrix(1,numberofpointsinx) = lambda;
        
        %building the last row
        
        RHSmatrix(numberofpointsinx,1) = lambda;
        RHSmatrix(numberofpointsinx,numberofpointsinx-1) = lambda;
        RHSmatrix(numberofpointsinx,numberofpointsinx) = 1 - 2*lambda;
        
        RHSknownvalues = zeros(numberofpointsinx,1);
        
        for j = 1:numberofpointsint
            
        for i = 1:numberofpointsinx    
        RHSknownvalues(i) = matrixu(j,i); 
        valuesoff(i) = dt*f(pointx(i),pointt(j));
        end
        
        %calculate the next iteration
        
        unknownvalues = (RHSmatrix*RHSknownvalues) + valuesoff;
        
        for i = 1:numberofpointsinx
        matrixu(j+1,i) = unknownvalues(i);
        end
        %Enforicing BC
        matrixu(j+1,numberofpointsinx+1) = matrixu(j+1,1);  
        
        
        end
        
        
        

        
        matrixofcorrectsolutions = zeros(numberofpointsint+1, numberofpointsinx+1);
        for j = 1:numberofpointsint +1
           for i = 1:numberofpointsinx +1
            matrixofcorrectsolutions(j,i) = actualsol(pointx(i),pointt(j));   
           end
        end
        
        errorsmatrix = abs(matrixu-matrixofcorrectsolutions);
        
        

 
end









