clear all;
close all;

N_x = 10;
h_x = 1/N_x;

alpha = @(t) t <= 1;

[A, b, ~] = build_A_du(N_x);
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

disp(' ')
%% Comparing matlab solvers

solvers = {'ode23', 'ode23s'};
delete('part4_solver_comparison.csv')
custom_fprintf('part4_solver_comparison.csv', '%s;%s;%s;%s;%s;%s;%s;%s;%s;%s\r\n', ...
               'N', [solvers{1} ' timesteps'], [solvers{2} ' timesteps'], ['ode23s opt timesteps'], ...
               [solvers{1} ' CPU-time'], [solvers{2} ' CPU-time'], ['ode23s opt CPU-time'], ...
               [solvers{1} ' Delta tmax'], [solvers{2} ' Delta tmax'], ['ode23s opt Delta tmax'])
for N_x = [10 20 40]
    options = odeset('RelTol',1e-6, 'AbsTol',1e-6);
    [~, ~, du] = build_A_du(N_x);
    u_0 = zeros(N_x,1);
    time = zeros(1,3);
    T = {[], [], []}; Y = {[], [], []};
    i = 1;
    for solver = solvers
        tic
        [T{i}, Y{i}] = feval(solver{:}, du, [0 2], u_0, options);
        ux_0 = arrayfun(alpha, T{i});
        Y{i} = [ux_0 Y{i}];
        time(i) = toc;
        i = i+1;
    end
    [A, ~, du] = build_A_du(N_x, 'sparse');
    options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6, 'Jpattern', A);
    tic
    [T{3}, Y{3}] = ode23s(du, [0 2], u_0, options);
    time(3) = toc;

    custom_fprintf('part4_solver_comparison.csv', '%d;%d;%d;%d;%f;%f;%f;%f;%f;%f\r\n', ...
                   N_x, length(T{1}), length(T{2}), length(T{3}), time(1), time(2), time(3), ...
                   max(diff(T{1})), max(diff(T{2})), max(diff(T{3})));
end
