N_x = 20;

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
for i = 2:N_x-1
    A(i,i) = -2;
    A(i,i+1) = 1;
    A(i+1,i) = 1;
end
A(end,end-1) = 2;
A(end,end) = -2;

A = 1/h_x^2*A;

alpha = @(t) t <= 1;
b = zeros(N_x,1); % Remember to replace b(1) in time loop
                  % b(1) = alpha(t)/h^2;

u_0 = zeros(N_x,1);
u_0(1) = 1;

N_t = 20;
h_t = 0.1
U = zeros(N_x,N_t+1);
U(:,1) = u_0;

for k = 1:N_t
    b(1) = alpha(k*h_t)/h_x^2;
    U(:,k+1) = (1 + h_t * A)*U(:,k) + h_t*b;
    disp(b)
end

