function [matrixu,errormatrix] = W2Q1explicitneumann(xbeginning,xend,tbeginning,tend,numberofpointsinx,numberofpointsint,eta,zeta1,zeta2,f,actualsol)


% for solving equations of the form u_t - u_xx = f(x,t)


%inputs

%xbeginning =0 ;
%xend = 2*pi;
%tbeginning =0 ;
%tend = 1;
%numberofpointsinx = 10;
%numberofpointsint= 10;
%eta = @(x) sin(x) ; %.... ;
%zeta1 = @(x,t) exp(-t);
%zeta2 = @(x,t) exp(-t);
%f = @(x,t) 0;
%actualsol = @(x,t) exp(-t)*sin(x);

matrixu = zeros(numberofpointsint+1, numberofpointsinx+1); % X's across t down (columns)

dx = (xend-xbeginning)/numberofpointsinx;
dt = (tend-tbeginning)/numberofpointsint;

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
        
        % - Initial conditions for eta (t = 0) including boundary points
        % for t=0
        
        for i = 1:numberofpointsinx+1
        matrixu(1,i) = eta(pointx(i));
        end

        
        %Build the Euler explicit algorithm.
        
        for j = 1:numberofpointsint
                for i = 2:numberofpointsinx
                    matrixu(j+1,i) = matrixu(j,i) + (dt)*f(pointx(i),pointt(j)) + (dt/(dx^2))*(matrixu(j,i-1) - 2*matrixu(j,i) + matrixu(j,i+1));
                    matrixu(j+1,1) = matrixu(j+1,2) - dx*zeta1(pointx(1),pointt(j+1));
                    matrixu(j+1,numberofpointsinx+1) = matrixu(j+1,numberofpointsinx) + dx*zeta2(pointx(numberofpointsinx+1),pointt(j+1));
                end
        end
        
        
        
        
        matrixofcorrectsolutions = zeros(numberofpointsint+1, numberofpointsinx+1);
        for j = 1:numberofpointsint +1
           for i = 1:numberofpointsinx +1
            matrixofcorrectsolutions(j,i) = actualsol(pointx(i),pointt(j));   
           end
        end
        
   
        errormatrix = abs(matrixu - matrixofcorrectsolutions);

 
end









