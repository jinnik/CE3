clear all;
close all;

N_x = 10;
h_x = 1/N_x;

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

%% Stable plot:
h_t = 0.5*h_x^2;
N_t = round(2/h_t);
U = zeros(N_x,N_t+1);
U(:,1) = u_0;

t = linspace(0, 2, N_t+1);
x = linspace(0, 1, N_x+1);
C = h_t/h_x^2;
ux_0 = arrayfun(alpha, t);
for k = 1:N_t
    U(:,k+1) = U(:,k) + h_t * A * U(:,k) + h_t / h_x^2 * b(t(k));
end

U = [ux_0; U];
surf(U)
disp('Stable:')
disp(['h_x = ' num2str(h_x)])
disp(['h_t = ' num2str(h_t)])
disp(['C = ' num2str(h_t/h_x^2)])

print3_data('part3_stable_plot.csv', x, t, U)

%% 2D plot for some time steps

figure()
leg={};
time_steps = [0.5 1 1.5 2];
Y = zeros(length(U(:,1)),length(time_steps));
for i = 1:length(time_steps)
    Y(:,i) = U(:,find(t == time_steps(i)));
    plot(x, Y(:,i))
    leg = {leg{:} ['\tau = ' num2str(time_steps(i))]};
    hold on
end
legend(leg);

print_data('part3_2D.csv', x, Y);

%% Unstable plot:
h_t = 1*h_x^2;
N_t = round(2/h_t);
U = zeros(N_x,N_t+1);
U(:,1) = u_0;

for k = 1:N_t
    ux_0(k) = alpha(k*h_t);
    U(:,k+1) = U(:,k) + h_t * A * U(:,k) + h_t / h_x^2 * b(k*h_t);
end
figure()
surf(U)
disp('Unstable:')
disp(['h_x = ' num2str(h_x)])
disp(['h_t = ' num2str(h_t)])
disp(['C = ' num2str(h_t/h_x^2)])

t = linspace(0, 2, N_t+1);
x = linspace(0, 1, N_x);
print3_data('part3_unstable_plot.csv', x, t, U)


%% Comparing matlab solvers

du = @(t, u) [A*u + 1/h_x^2 * b(t)];
tic
[T,Y] = ode23(du, [0 2], u_0);
time_ode23 = toc;
ux_0 = arrayfun(alpha, T);

Y = [ux_0 Y];

tic
[Ts,Ys] = ode23s(du, [0 2], u_0);
time_ode23s = toc;
ux_0 = arrayfun(alpha, Ts);
Ys = [ux_0 Ys];
disp('Solver stats: ')
disp(['N_t for ode23: '            num2str(length(T))    ])
disp(['N_t for ode23s: '           num2str(length(Ts))   ])
disp(['time for ode23: '           num2str(time_ode23)   ])
disp(['time for ode23s: '          num2str(time_ode23s)  ])
disp(['max time step for ode23: '  num2str(max(diff(T))) ])
disp(['max time step for ode23s: ' num2str(max(diff(Ts)))])

fprintf('%s %s %s %s\r\n', '$N$', 'timesteps' ,'CPU-time', '$\Delta t_{max}$')

% Work in progress..
if 0
for N_x = [10 20 40]
    du = build_du(N_x);
    % du = @(t, u) [A*u + N_x^2 * b(t)];

    tic
    [T,Y] = ode23(du, [0 2], u_0); % Could use feval to choose between
                                   % ode23 and ode23s
    time_ode23 = toc;
    tic
    [Ts,Ys] = ode23s(du, [0 2], u_0);
    time_ode23s = toc;

    fprintf('%d %f %f %f %f %f %f\r\n', N, length(T), length(Ts), time_ode23, ...
            time_ode23s, max(diff(T)), max(diff(Ts)));
end
end
