%% MATH 543 - Homework 4
% By Parham Khodadi
clear; clc; close all;

%% Function for Modified Gram-Schmidt (MGS) QR-factorization
% As found on page 58 (Algorithm 8.1) of Textbook

function [Q,R] = qr_mgs(A)
    [m,n] = size(A);
    Q = zeros(m,n,class(A));
    R = zeros(n,n,class(A));
    V = A; 
    
    for i = 1:n
        R(i,i) = norm(V(:, i)); 
        Q(:,i) = V(:, i) / R(i,i);
        
        for j = (i+1):n
            R(i,j) = Q(:, i)' * V(:, j);
            V(:, j) = V(:, j) - R(i,j) * Q(:, i);
        end
    end
end

%% Classical Gram–Schmidt Function
% Based on page 23 of Notes #6.
% Copied over from HW3
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

%% Experiment #1
x = (-128:128)'/128;
A = [x.^0, x.^1, x.^2, x.^3];
[Q,R] = qr(A,0);

scale = Q(257,:);
Q = Q*diag(1./scale);

% Check if it matches Figure 7.1
plot(Q);
title("Experiment #1")
grid on

%% Experiment #2
[U,X] = qr(randn(80));
[V,X] = qr(randn(80));
S = diag(2.^(-1:-1:-80));

A = U*S*V;

[QC,RC] = clgs(A);
[QM,RM] = mgs(A);