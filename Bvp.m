function  boundaryvalueproblem
%% Author: Muhammad Asad Ullah
 clear all; clc; close all;
 global  Nb Nt We L1 L2 M L3  lambda Sc Kp EC beta Q gamma
 % lambda = 0.1;
We = 0.1; L1 = 0; L2 = 0;  L3 = 0; gamma =0.1;Nt = 0.1;  M = 0; Sc = 0.2; Kp = 0.5; 
Pr = 0.3;  Nb = 0.1; EC = 0.2; beta = 0; Q = -1.5;

%% Initial values, solution and plot

for i=1:4   
if (i==1)
lambda=0.1;
elseif (i==2)
lambda=0.3;
elseif (i==3)
lambda=0.6;
elseif (i==4)
lambda=1;
end
           
sol = bvpinit(linspace(0,20,1000), [1e-8 1e-8 1e-8 1e-8 1e-8 1e-8 1e-8]);
sol1 = bvp4c(@bvpexam2, @bcexam2, sol);
x1 = sol1.x;
y1 = (sol1.y);

figure (1)
plot(x1,y1(2,:),'LineWidth',1.6)
xlabel('\eta')
ylabel('f^{\prime}(\eta)')
axis([0 20 -inf inf])
legend('\lambda=0.1','\lambda=0.3','\lambda=0.6','\lambda=1')
hold on

figure (2)
plot(x1,y1(4,:),'LineWidth',1.6)
xlabel('\eta')
ylabel('\theta(\eta)')
axis([0 10 -inf inf])
legend('\lambda=0.1','\lambda=0.3','\lambda=0.6','\lambda=1')
hold on

figure (3)
plot(x1,y1(6,:),'LineWidth',1.6)
xlabel('\eta')
ylabel('\phi(\eta)')
axis([0 20 -inf inf])
legend('\lambda=0.1','\lambda=0.3','\lambda=0.6','\lambda=1')
hold on


end

%% Define boundary conditions
    function res = bcexam2(y0, yinf)
    res = [y0(1); y0(2)-L1*(y0(3))-1; y0(4)-L2*(y0(5))-1; y0(6)-L3*(y0(7))-1; yinf(2)-gamma; yinf(4); yinf(6);];
    end
 
%% First order ODEs are define here
 
    function ysol = bvpexam2(x,y)
     yy1 = (lambda*(y(2)+x/2*y(3))+y(2).*y(2)-y(1).*y(3)+(M*M+Kp)*y(2)-(M*M+lambda)*gamma-gamma*gamma)/(1+We*y(3));
        yy2 = -Pr*(y(1).*y(5)-2*y(2).*y(4)-(lambda/2)*(x*y(5)+3*y(4))+Nb*y(5).*y(7)+Nt*y(5).*y(5)+Q*y(4)+EC*(y(3).*y(3)+We*y(3).*y(3).*y(3)+M*M*y(2).*y(2)));
        yy3 =(-Nt/Nb)*yy2-Sc*(y(1).*y(7)-2*y(2).*y(6)-lambda/2*(x*y(7)+3*y(6))-beta*y(6));
        ysol = [y(2);y(3);yy1;y(5);yy2;y(7);yy3;];   
     end
end 