clc; clear all; close all;


N=20;  % total number of stepsizes
h = 1/(N+1); % stepsize
global N h
tic

u0(1)=1; % initial value for xi=0
for i=2:N
    u0(i)=0; % initial values for xi>0
end

t0 = 0; tend = 2; % start and end time
options = odeset('RelTol',1e-6, 'AbsTol',1e-6);
[t,y] = ode23s(@matsolv,[t0,tend],u0,options); % solving with ode23 or ode23s
b = length(t);

toc

function [U] = matsolv(t,u)

global N h
U = zeros(N,1);
U(1) = 1/h^2*(-2*u(1)+u(2));
U(2:N-1) = 1/h^2*(u(1:N-2)-2*u(2:N-1)+u(3:N));
U(N) = 1/h^2*(2*u(N-1)-2*u(N));

end