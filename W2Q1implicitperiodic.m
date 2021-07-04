function [matrixu,errormatrix,matrixofcorrectsolutions] = W2Q1implicitperiodic(xbeginning,xend,tbeginning,tend,numberofpointsinx,numberofpointsint,eta,f,actualsol)


% for solving equations of the form u_t - u_xx = f(x,t)


matrixu = zeros(numberofpointsint+1, numberofpointsinx+1); % X's across t down (columns)

dx = (xend-xbeginning)/numberofpointsinx;
dt = (tend-tbeginning)/numberofpointsint;
lambda =  dt/(dx^2);

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
        
        LHSmatrix = zeros(numberofpointsinx,numberofpointsinx); %since we are reducing dimensionality for bc.
        
        for i=2: numberofpointsinx-1
            LHSmatrix(i,i-1) = -lambda;
            LHSmatrix(i,i) = 1 + 2*lambda;
            LHSmatrix(i,i+1) = -lambda;
        end
        
        %building the first row
        
        LHSmatrix(1,1) = 1 + 2*lambda;
        LHSmatrix(1,2) = -lambda;
        LHSmatrix(1,numberofpointsinx) = -lambda;
        
        %building the last row
        
        LHSmatrix(numberofpointsinx,1) = -lambda;
        LHSmatrix(numberofpointsinx,numberofpointsinx-1) = -lambda;
        LHSmatrix(numberofpointsinx,numberofpointsinx) = 1 + 2*lambda;
        
        LHSinverted = inv(LHSmatrix);
        RHSknownvalues = zeros(numberofpointsinx,1);
        
        for j=1:numberofpointsint
            
            
        for i = 1: numberofpointsinx
        RHSknownvalues(i) = matrixu(j,i); 
        valuesoff(i) = dt*f(pointx(i),pointt(j));
        end
        
        RHSside = RHSknownvalues + valuesoff;
        Newvalues = LHSinverted*RHSside;
        
        for i = 1: numberofpointsinx
        matrixu(j+1,i) = Newvalues(i);
        end
        matrixu(j+1,numberofpointsinx+1) = matrixu(j+1,1);
        
        end
        
        
        
        matrixofcorrectsolutions = zeros(numberofpointsint+1, numberofpointsinx+1);
        for j = 1:numberofpointsint +1
           for i = 1:numberofpointsinx +1
            matrixofcorrectsolutions(j,i) = actualsol(pointx(i),pointt(j));   
           end
        end
   
        
        errormatrix = abs(matrixofcorrectsolutions-matrixu);
        
  
end