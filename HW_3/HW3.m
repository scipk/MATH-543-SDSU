%% MATH 543 - Homework 3
% By Parham Khodadi
clear; clc; close all;

%% Classical Gram–Schmidt Function
% Based on page 23 of Notes #6
function [Q,R] = qr_cgs(A)
    [m,n] = size(A);
    Q = zeros(m,n,class(A));
    R = zeros(n,n,class(A));

    for k = 1:n
        v = A(:,k);

        for i = 1:k-1
            R(i,k) = Q(:,i)' * A(:,k);
            v = v - R(i,k) * Q(:,i);
        end

        R(k,k) = norm(v,2);

        % Avoid dividing by zero when setting Q(:,k) to v/R(k,k)
        if R(k,k) == 0
            Q(:,k) = zeros(m,1,class(A));
        else
            Q(:,k) = v / R(k,k);
        end
    end
end

%% Validation

sizes = [3,5,251];

for s = sizes
    fprintf('\n------ Test size: %d x %d ------\n', s,s);
    A = rand(s);
    [Q, R] = qr_cgs(A);
    
    % (i) A - Q*R = 0
    err_factor = norm(A - Q*R, 'fro');

    fprintf('||A - Q*R||_F      = %.4e\n', err_factor);

    % (ii) Q'*Q = I
    err_orth = norm(Q'*Q - eye(s), 'fro');
    fprintf('||Q''*Q - I||_F     = %.4e\n', err_orth);

    % (iii) strictly-lower part = 0
    err_upper = norm(tril(R,-1), 'fro');
    fprintf('||tril(R,-1)||_F   = %.4e\n', err_upper);

    % Part 3
    if s <= 5
        fprintf('--- MATLAB qr(A,0) ---\n');

        [Qm,Rm] = qr(A,0);

        err_factor_matlab = norm(A - Qm*Rm, 'fro');
        fprintf('||A - Qm*Rm||_F                        = %.3e\n', err_factor_matlab);
        
        err_Q_diff = norm(abs(Q) - abs(Qm), 'fro');
        fprintf('Sign-insensitive ||abs(Q)-abs(Qm)||_F  = %.3e\n', err_Q_diff);

        err_R_diff = norm(abs(R) - abs(Rm), 'fro');
        fprintf('Sign-insensitive ||abs(R)-abs(Rm)||_F  = %.3e\n', err_R_diff);
    end
end


%% Find a non-zero matrix where CGS breaks
fprintf('\n------ Part 4: CGS break example (Hilbert matrix) ------\n');

n = 12;
A = hilb(n);

[Q,R] = qr_cgs(A);
fprintf('||A - Q*R||_F  = %.3e\n', norm(A - Q*R,'fro'));
fprintf('||Q''*Q - I||_F = %.3e\n', norm(Q'*Q - eye(n),'fro'));
