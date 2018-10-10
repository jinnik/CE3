clear all;
close all;

N_x = 10;

h_x = 1/N_x;

tau = 0;

U =    ones(N_x, 1);
D = -2*ones(N_x, 1);
L =    ones(N_x, 1);
L(end) = 2;

% How did spdiags work again?!
% A = spdiags([U D L], [1 0 -1])

A = zeros(N_x,N_x);
A(1,1) = -2;
A(1,2) = 1;
A(2,1) = 1;
for i = 2:N_x-1
    A(i,i) = -2;
    A(i,i+1) = 1;
    A(i+1,i) = 1;
end
A(end,end-1) = 2;
A(end,end) = -2;

A = 1/h_x^2*A;

alpha = @(t) t <= 1;
b = @(t) [alpha(t); zeros(N_x-1, 1)];

u_0 = zeros(N_x,1);

N_t = 1000;
h_t = 0.5*h_x^2;
U = zeros(N_x,N_t+1);
U(:,1) = u_0;

C = h_t/h_x^2;

for k = 1:N_t
    ux_0(k) = alpha(k*h_t);
    U(:,k+1) = U(:,k) + h_t * A * U(:,k) + h_t / h_x^2 * b(k*h_t);
end

U = [[1 ux_0]; U];

if 1
%% Comparing with matlab solvers

du = @(t, u) [A*u + 1/h_x^2 * b(t)];

[T,Y] = ode23(du, [0 2], u_0);
ux_0 = zeros(length(T),1);
for i = 1:length(T)
    if T(i) > 1
        break
    end
    ux_0(i) = alpha(T(i));
end

Y = [ux_0 Y];

end
