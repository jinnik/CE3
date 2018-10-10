clc; clear all; close all;

dx=0.2; % stepsize in space
n=1/dx; % total number of stepsizes
x=linspace(0,1,n+1); % space discretization

dt=0.01; % timestep
t=200; % total number of timesteps
c=(dt)/((dx)^2); % courant number

[U] = mol(t,dt,dx); 

surf(U)
xlabel('time')
ylabel('x')
zlabel('u')

tau = [0.5 1 1.5 2];
nt = tau/dt;
for j = 1:4
    y = zeros(n+1,1);
    y(:,j) = U(:,nt(j));
    figure()
    plot(x,y(:,j),'.-')
end

function [U] = mol(T,dt,h)

N=1/h; % total number of stepsizes
c=(dt)/(h^2); % courant number
u0(1)=1; % initial value for xi=0
u0(2:N+1)=0; % initial values for xi>0
u=u0'; % initial value matrix
b0(1:T)=1;  % boundary values for tau<1
b=b0; % upper boundary condition matrix

uin = zeros(1,N+1); % solving at each stepsize
U=zeros(N+1,T); % solution matrix U
U(:,1)=u; % add initial condition to matrix U
U(1,:)=b; % add boundary condition to matrix U

for i=2:T
    uin = u; % saving the backwards time line
    if (dt*i)>1 % boundary values for tau>1
        u(1)=0;
    end
    for j=2:N % the rest of the solution matrix
        u(j)=c*uin(j+1)-(2*c-1)*uin(j)+c*uin(j-1);
        U(:,i)=u;
    end
    u(N+1)=u(N-1); %setting the lower boundary condition
end

end






