%   2D Helmholtz Equation
%   Daniel Cariel
%       This code 


%SURFACE REGION FOR HELMHOLTZ EQUATION 
clc; clear all;
n=input('Input value of n:'); %Mesh size for Gauss Seidel approximation. Number of nodes
a_x= -pi; %Lower boundary in the x-axis
b_x= pi;  %Upper boundary in the y-axis 
a_y= -pi; %Lower boundary in the x-axis
b_y= pi;  %Upper boundary in the y-axis
x= linspace(a_x,b_y,n); %x vector describing nodes in the x-axis
y= linspace(a_y,b_y,n); %y vector describing nodes in the y-axis

%Discretization Coefficients and initial guess 
gamma=-pi;      %Wave constant 
h= (2*b_x)/n; %Discretization step
u= zeros(n);  %Initial guess for gauss seidel approximation 

%%BOUNDARY CONDITIONS
%Dirishlet:

u(1,:)= b_y.*((b_y-a_y).^2)+ ((x(:)-a_x)./(b_x-a_x)).*(((b_y-a_y).^2).*cos(pi.*(b_y/a_y))-b_y.*(b_y-a_y).^2); %Boundary condition at the top edge of the prescribed region

u(:,1)= y(:).*(y(:)-a_y).^2; %Boundary condition at the left edge of the prescribed region 

u(:,n)= ((y(:)-a_y).^2).*cos(pi.*(y(:)/a_y)); %Boundary condition at the right edge of the prescribed region 
 

%GAUSS SEIDEL ITERATIONS
F=zeros(n); %Initia matrix for the forcing function of in the Helmholtz Equation 
error=1; %Initial value of the error for the Gauss Seidel numerical solver
iteration=0; %Iteration counter intended to keep track of the convergence rate of the code 

% CHECKPOINTING 
% Sometimes files take a long time to run to completion. As a result, sometimes they crash due to a variety of reasons: power failure, walltime limit, scheduled shutdown, etc. 
% Checkpoint/Restarting has long been a common technique to tackle this issue. Checkpointing/Restarting essentially means saving data to disk periodically so that, if need be, 
% you can restart the job from the point at which your data was last saved. 

frequency=10;
tic; %measure time performance of the code

% Before the start of the iteration loop, "check-in" each variable
% that should be checkpointed in the event of restarting the job

matfile = 'PoissonEquationSolution.mat';     % mandatory; name of checkpoint mat-file
s = struct();                                % mandatory; create struct for checkpointing
s = chkin(s,{'iteration'});                       % mandatory; iter is iteration loop index
s = chkin(s,{'frequency'});                  % mandatory; frequency is checkpointing period 
                                             % i.e., how often to perform a save

% continue until all variables are checked in. Note that you are only
% checking in the variables, they don't need to have been already defined

chkNames = fieldnames(s);    % the full list of variables to checkpoint
nNames = length(chkNames);   % number of variables in list


while max(error(:))>=1e-6   %Gauss Seidel iterations will continue until the error is within the desired tolerance 1e-6
  iteration=iteration+1;    %Iteration step
  
  
  % If you want to test the restart script, use the function pause(1) to slow down the while loop.
    % This will slow down the while loop to 1 sec per iteration so that ctrl + C can used be to
    % "kill" the code to simulate a computer crash. From there, use the restart script to restart the loop.  
    % pause(.05)
    if mod(iteration, frequency) == 0 % If statement, checkpoints periodically (determined by the frequency)
        chkpt                    % chkpt script performs checkpointing (save) every *frequency* iterations
        fprintf(1, ['Checkpointing frequency is every %2d iterations.' ...
          'Data updated at iteration %3d\n'], ...
          frequency, iteration);      % Confirm after each checkpointing event 
    end
  
   
  u_0=u;                    %Solution of u from previous iteration. It is used to calculate the error
  
    for i=2:n-1
        for j=2:n-1
          F(i,j)= cos((pi/2)*(2*((x(i)-a_x)/(b_x-a_x))+1))*sin((pi*y(j)-a_y)/(b_y-a_y));    %Forcing function describing the problem 
          u(i,j)= (1/((gamma*h^2)-4))*((h^2)*F(i,j)-(u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1))); % Discretization solution of the 2D Helmholtz equation
          u(n,j)= (1/((gamma*h^2)-4))*((h^2)*F(i,j)-(u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j+1)));  %Neumann boundary condition applied to the bottom edge of the prescribed region 
        end 
    end
     u_f=u;                             %Solution of u from most recent iteration. It is compared with the previous iteration to calculate the error
     error= abs((u_f-u_0)./(u_f));      %Error calculation. The error drives the Gauss Seidel solver until tolerance is reached
end 
timedoc=toc

%PLOT OF THE RESULTS
% 2D contour of the surface
figure      
contourf(u) 
colorbar('location','eastoutside','fontSize',11);
xlabel('X Number of Nodes in X-direction','fontSize',11); %3D figure lable on the x axis
ylabel('Y Number of Nodes in Y-direction','fontSize',11); %3D figure lable on the y axis
title('Gauss Seidel for Helmhotlz') 
% 3D surface graph of the solution 
figure
surf(x,y,u,'EdgeColor','none')
xlabel('X Number of Nodes in X-direction','fontSize',11); %3D figure lable on the x axis
ylabel('Y Number of Nodes in Y-direction','fontSize',11); %3D figure lable on the y axis
zlabel('Position U','fontSize',12);                       %3D figure lable on the z axis
title('2D Helmholtz Solution -  Gauss Seidel Numerical Solver');