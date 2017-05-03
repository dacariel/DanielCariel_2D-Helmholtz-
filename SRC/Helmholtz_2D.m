%   2D Helmholtz Equation
%   Daniel Cariel
%       This codes..... Add description 


%Region for the Helmholtz Equation
n=input( 'Input value of n:');
a_x= -pi;
b_x= pi;
a_y= -pi;
b_y= pi;
x= linspace(a_x,b_y,n);
y= linspace(a_y,b_y,n);
 
 
 
% Boundary Conditions
 
 
 
%  GAUSS SEIDEL NUMERICAL SOLVER 
 
 u= zeros(n); % Initial guest for gauss seidel approximation 
 
 %Boundary Conditions
 for i=1:n
    u(i,1)= b_y*((b_y-a_y)^2)+ ((x(i)-a_x)/(b_x-a_x))* (((b_y-a_y)^2)*cos(pi*(b_y/a_y))-b_y*(b_y-a_y)^2);
 
    u(i,n)= a_y;
 end 
 for j=1:n
    u(1,j)= ((y(j)-a_y)^2)*cos(pi*(y(j)/a_y));
 
     u(n,j)= y(j)*(y(j)-a_y)^2;
 end 

 
 
 
 
 