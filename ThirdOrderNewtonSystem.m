% Solve for x1 and x2 to minimize
% f(x1, x2) =sum_{j=1}^{N} [y_j - exp(x1 t_j) - exp(x2 u_j)]^2
clear all;
 N = 200;
 bvect = zeros(N,1);
 seed = 101;
 rng(seed);
 Amat = 5*rand(N,2);
 for j = 1:N
 t = Amat(j,1);
 u = Amat(j,2);
 bvect(j) = exp(0.51*t)+exp(0.51*u);
 end
% Solve the nonlinear system of equations
 initial = [0.75; 0.25];
 NumIter = 10; % number of iterations
xvect = initial
for iter = 1:NumIter

 x1 = xvect(1);
 x2 = xvect(2);

 sum1 = 0;
 for j=1:N
 y = bvect(j);
 t = Amat(j,1);
 u = Amat(j,2);
 sum1 = sum1 + ((y - exp(x1*t) - exp(x2*u))*(-1*exp(x1*t)*t));
 end

 sum2 = 0;
 for j=1:N
 y = bvect(j);
 t = Amat(j,1);
 u = Amat(j,2);
 sum2 = sum2 + ((y - exp(x1*t) - exp(x2*u))*(-1*exp(x2*u)*u));
 end
 Fvect_x = [sum1; sum2];


 sum3 = 0;
 for j=1:N
 y = bvect(j);
 t = Amat(j,1);
 u = Amat(j,2);
 sum3 = sum3 + (-1*exp(x1*t)*t)*(-1*exp(x1*t)*t)+(-1*exp(x1*t)*t^2)*(y - exp(x1*t) - exp(x2*u));
 end

 sum4 = 0;
 for j=1:N
 y = bvect(j);
 t = Amat(j,1);
 u = Amat(j,2);
 sum4 = sum4 + (exp(x2*u)*u)*(exp(x1*t)*t);
 end

 sum5 = 0;
 for j=1:N
 y = bvect(j);
 t = Amat(j,1);
 u = Amat(j,2);
 sum5 = sum5 + (-1*exp(x2*u)*u)*(-1*exp(x2*u)*u)+(-1*exp(x2*u)*u^2)*(y - exp(x1*t) - exp(x2*u));
 end

 F_prime_x = [ sum3 sum4; ...
 sum4 sum5];

 yvect = xvect - inv(F_prime_x)*Fvect_x ;

 y1 = yvect(1);
 y2 = yvect(2);
 sum6 = 0;
 for j=1:N
 y = bvect(j);
 t = Amat(j,1);
 u = Amat(j,2);
 sum6 = sum6 + (y - exp(y1*t) - exp(y2*u))*(-1*exp(y1*t)*t);
 end

 sum7 = 0;
 for j=1:N
 y = bvect(j);
 t = Amat(j,1);
 u = Amat(j,2);
 sum7 = sum7 + (y - exp(y1*t) - exp(y2*u))*(-1*exp(y2*u)*t);
 end

 Fvect_y = [sum6; sum7];

 zvect = xvect - inv(F_prime_x)*( Fvect_x + Fvect_y );
 xvect = zvect;

 it = num2str(iter);
 message = ['at iteration ', it, ', we have [x, y] ='];
 disp(message)
 disp(round(xvect,10))
end
