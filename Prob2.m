function Prob2()

%% Numerical Solution

L=5;
h=0.05;
D = @(x) 0.05*x + 0.84;
x=linspace(0,L,L/h+1);
A=zeros(L/h-1,L/h+1);
C=zeros(L/h+1,1);
d=zeros(L/h-1,1);

% eqn reduced to 
% ai*C(i-1) + bi*C(i) + ci*C(i+1) = 0;
% where

ai= @(x) D(x)/(h^2)-0.05/h;
bi= @(x) 0.05/h-2*D(x)/(h^2)-0.1;
ci= @(x) D(x)/(h^2);

% boundary conditions;
% at x=0 C=100;
% at x=L C'=0 => C'=( C(L)-C(L-h) )/h =0 => C(L) = C(L-h)
C(1)=100;

for i = 1:(L/h-1)
    A(i, i) = ai(i*h);
    A(i, i+1) = bi(i*h);
    A(i, i+2) = ci(i*h);
end

d(1,1)=A(1,1)*-100;
A(L/h-1,L/h)=A(L/h-1,L/h)+A(L/h-1,L/h+1);
C(2:L/h)=Gauss(A(:,2:L/h), d);
C(L/h+1)=C(L/h);


plot(x,C);
hold on;

%% Cubic Spline
[a,b,c,d]=findabcd(x,C);
p(x,a,b,c,d);

%% Analytical Solution

L = 5;
k = 0.1;
CA0 = 100;
D0 = 0.84;
DL = 1.09;
a = (DL - D0) / L;

alpha = 2 * sqrt(k) / a;

z0 = alpha * sqrt(D0);
zL = alpha * sqrt(DL);

I0_z0 = besseli(0, z0);
K0_z0 = besselk(0, z0);
I1_zL = besseli(1, zL);
K1_zL = besselk(1, zL);

C2_over_C1 = I1_zL / K1_zL;

C1_denominator = I0_z0 + C2_over_C1 * K0_z0;
C1 = CA0 / C1_denominator;

C2 = C1 * C2_over_C1;

x_fine = linspace(0, L, L/h+1);

z_fine = alpha * sqrt(D0 + a * x_fine);

CA_analy = C1 * besseli(0, z_fine) + C2 * besselk(0, z_fine);

figure;
title('Numerical solution Cubic Spline');
hAnaly=plot(x_fine, CA_analy,'r');
legend(hAnaly,'Analytical Solution (Red)');
xlabel('L (cm)');
ylabel('Concentration');
