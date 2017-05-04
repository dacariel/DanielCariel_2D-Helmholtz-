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
gamma=-1;      %Wave constant 
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

while max(error(:))>=1e-6   %Gauss Seidel iterations will continue until the error is within the desired tolerance 1e-6
  iteration=iteration+1;    %Iteration step
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