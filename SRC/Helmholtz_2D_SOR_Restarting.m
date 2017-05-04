%% Restarting
matfile = 'PoissonEquationSolution';   % should match that in test_checkpoint.m
load(matfile);        % retrieve data from matfile
iter1 = iteration+1;       % iter is the last time test_checkpoint issued
                      % a save; we start computing on the next step

%% Iterative Looping- Gauss-Seidel Method                      
%% 
%GAUSS SEIDEL ITERATIONS
F=zeros(n); %Initia matrix for the forcing function of in the Helmholtz Equation 
error=1; %Initial value of the error for the Gauss Seidel numerical solver


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
          u(i,j)= (B/((gamma*h^2)-4))*((h^2)*F(i,j)-(u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1)))+(1-B).*u(i,j); % Discretization solution of the 2D Helmholtz equation
          u(n,j)= (1/((gamma*h^2)-4))*((h^2)*F(i,j)-(u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j+1)));  %Neumann boundary condition applied to the bottom edge of the prescribed region 
        end 
    end
     u_f=u;                             %Solution of u from most recent iteration. It is compared with the previous iteration to calculate the error
     error= abs((u_f-u_0)./(u_f));      %Error calculation. The error drives the Gauss Seidel solver until tolerance is reached
end 
%%
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
title('2D Helmholtz Solution -  Gauss Seidel with Successive Over Relaxation');