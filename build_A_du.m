function [A, b, du] = build_A_du(N_x, varargin)
    %BUILD_A_DU Create A, b and du
    %   Given N_x, the matrix A, vector b and vector
    %   du can be created. If varargin is given then
    %   the matrix is created as sparse.

    h_x = 1/N_x;

    if ~isempty(varargin)
        Diag = -2*ones(1, N_x);
        Upper =   ones(1,N_x);
        Lower =   ones(1,N_x-2);
        Lower = [Lower 2 0];
        A = spdiags([Diag' Upper' Lower'], [0 1 -1], N_x, N_x);
        A = 1/h_x^2*A;
    else
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
    end
    alpha = @(t) t <= 1; % Heat pulse ends at t=1
    b = @(t) [alpha(t); zeros(N_x-1, 1)]; % b(1) = alpha(t)
                                          % is the BC at the
                                          % left-most point
    du = @(t, u) [A*u + N_x^2 * b(t)];
end
