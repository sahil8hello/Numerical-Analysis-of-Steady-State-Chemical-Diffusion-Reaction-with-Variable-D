function [a,b,c,d] = findabcd(x,y)

[~,n] = size(x);

hi = zeros(1,n-1);
for i=1:n-1
    hi(i) = x(i+1) - x(i);
end

A = zeros(n);
p = zeros(n,1);

% Natural boundary conditions
A(1,1) = 1;
A(n,n) = 1;

for i = 2:n-1
    A(i,i-1) = hi(i-1);
    A(i,i) = 2*(hi(i) + hi(i-1));
    A(i,i+1) = hi(i);
    
    p(i) = 3*((y(i+1)-y(i))/hi(i) - (y(i)-y(i-1))/hi(i-1));
end

% Solve for c coefficients
c = Gauss(A, p);


% a coefficients (same as y)
a = y(:);

% compute b and d
b = zeros(n,1);
d = zeros(n,1);

for i = 1:n-1
    b(i) = (y(i+1)-y(i))/hi(i) - hi(i)*( 2*c(i) + c(i+1) )/3;
    d(i) = (c(i+1) - c(i)) / (3*hi(i));
end