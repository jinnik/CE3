function [A, b, du] = build_A_du(N_x)
    %BUILD_A_DU Create A, b and du
    %   Given N_x, the matrix A, vector b and vector
    %   du can be created.
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
    A(end,end-1) = 2; % BC isolated end of rod
    A(end,end) = -2; % BC isolated end of rod
    A = 1/h_x^2*A;

    alpha = @(t) t <= 1; % Heat pulse ends at t=1
    b = @(t) [alpha(t); zeros(N_x-1, 1)]; % b(1) = alpha(t)
                                          % is the BC at the
                                          % left-most point
    du = @(t, u) [A*u + N_x^2 * b(t)];
end
