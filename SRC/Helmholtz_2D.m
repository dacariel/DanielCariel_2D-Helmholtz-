%   2D Helmholtz Equation
%   Daniel Cariel
%       This codes..... Add description 


%Region for the Helmholtz Equation
clc; clear all;
n=input( 'Input value of n:');
a_x= -pi;
b_x= pi;
a_y= -pi;
b_y= pi;
x= linspace(a_x,b_y,n);
y= linspace(a_y,b_y,n);

%Discretization Coefficients
   gamma= pi;
   h= (b_x)/n;
   
 
   
 
 
 
% Boundary Conditions
 
 
 
%  GAUSS SEIDEL NUMERICAL SOLVER 
 
 u= zeros(n); % Initial guess for gauss seidel approximation 
 
 %Boundary Conditions
    %Dirishlet 
    u(:,1)= b_y.*((b_y-a_y).^2)+ ((x(:)-a_x)/(b_x-a_x)).* (((b_y-a_y).^2).*cos(pi.*(b_y/a_y))-b_y.*(b_y-a_y).^2);

    u(1,:)= ((y(:)-a_y).^2).*cos(pi.*(y(:)/a_y));
 
    u(n,:)= y(:).*(y(:)-a_y).^2;


 
 %Gauss Seidel Iterations
 error=1;
while max(max(error)) > 1e-6 
 u_0=u;
 for i=2:n-1
     for j=2:n-1
          F(i,j)= cos((pi/2)*(2*((x(i)-a_x)/(b_x-a_x))+1))*sin((pi*y(j)-a_y)/(b_y-a_y));
    %     F(i,j)=0;
         u(i,j)= (1/((gamma*h^2)-4))*((h^2)*F(i,j)-(u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1)));
     end 
     u(i,n)= (1/((gamma*h^2)-4))*((h^2)*F(i,j)-(u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j+1))); 
 end 
     error= abs((u_0-u)./(u));
end 

surf(x,y,u)
% figure
% contourf(u)
%      
 
 
 
 