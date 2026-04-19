%% MATH 543 HW 8
% by Parham Khodadi

clear; clc; close all;

%% Problem 1

rng(42);
m    = 10;
nPts = 500;
t    = linspace(0, 20, nPts);

% Pre-generate all 12 matrices so both figures use the same ones
matrices = cell(1,12);
for trial = 1:12
    matrices{trial} = randn(m) - 2*eye(m);
end

% Semilogy plots: ||e^{tA}||_2 vs e^{t*alpha(A)}
f1 = figure('Position', [100 100 1400 900]);
for trial = 1:12
    A       = matrices{trial};
    eigvals = eig(A);
    alphaA  = max(real(eigvals));
    [~,idx] = max(real(eigvals));
    domEig  = eigvals(idx);

    normEtA = arrayfun(@(s) norm(expm(s*A),2), t);
    refLine = exp(t * alphaA);

    if abs(imag(domEig)) > 1e-10, tag = 'COMPLEX'; else, tag = 'real'; end

    subplot(3,4,trial);
    semilogy(t, normEtA, 'b-', 'LineWidth',1.5); hold on;
    semilogy(t, refLine, 'r--','LineWidth',1.2); hold off; grid on;
    title(sprintf('Trial %d: \\alpha=%.2f (%s)', trial, alphaA, tag),'FontSize',9);
    xlabel('t'); ylabel('||e^{tA}||_2');
    legend('||e^{tA}||_2','e^{t\alpha}','Location','best','FontSize',7);
end
sgtitle('24.3: ||e^{tA}||_2 vs e^{t\alpha(A)}','FontSize',14,'FontWeight','bold');
print(f1,'Figures/P1_main','-depsc2');

% Ratio plots: ||e^{tA}||_2 / e^{t*alpha(A)}
f2 = figure('Position', [100 100 1400 900]);
for trial = 1:12
    A       = matrices{trial};
    eigvals = eig(A);
    alphaA  = max(real(eigvals));
    [~,idx] = max(real(eigvals));
    domEig  = eigvals(idx);

    normEtA = arrayfun(@(s) norm(expm(s*A),2), t);
    ratio   = normEtA ./ exp(t*alphaA);

    subplot(3,4,trial);
    plot(t, ratio, 'b-', 'LineWidth',1.2); grid on;
    title(sprintf('dom eig = %.2f%+.2fi', real(domEig), imag(domEig)),'FontSize',8);
    xlabel('t'); ylabel('ratio');
end
sgtitle('Ratio ||e^{tA}||_2 / e^{t\alpha(A)}','FontSize',14,'FontWeight','bold');
print(f2,'Figures/P1_ratio','-depsc2');

%% Problem 2

fprintf('=== Problem 2: Householder Hessenberg ===\n\n');

% --- 5x5 ---
rng(100); A5 = randn(5);
[H5, Q5] = householder_hessenberg(A5);
fprintf('--- 5x5 ---\n');
fprintf('H =\n'); disp(H5);
fprintf('||Q*H*Q''-A|| = %.2e,  ||Q''Q-I|| = %.2e\n\n', ...
    norm(Q5*H5*Q5'-A5), norm(Q5'*Q5-eye(5)));

% --- 7x7 ---
rng(200); A7 = randn(7);
[H7, Q7] = householder_hessenberg(A7);
fprintf('--- 7x7 ---\n');
fprintf('H =\n'); disp(H7);
fprintf('||Q*H*Q''-A|| = %.2e,  ||Q''Q-I|| = %.2e\n\n', ...
    norm(Q7*H7*Q7'-A7), norm(Q7'*Q7-eye(7)));

% --- 9x9 validation vs hess() ---
rng(300); A9 = randn(9);
[H9, Q9]       = householder_hessenberg(A9);
[Q9lib, H9lib] = hess(A9);

fprintf('--- 9x9 Validation vs hess() ---\n');
fprintf('              Ours          hess()\n');
fprintf('||QHQ''-A||   %.2e    %.2e\n', norm(Q9*H9*Q9'-A9), norm(Q9lib*H9lib*Q9lib'-A9));
fprintf('||Q''Q-I||    %.2e    %.2e\n', norm(Q9'*Q9-eye(9)), norm(Q9lib'*Q9lib-eye(9)));
fprintf('max below    %.2e    %.2e\n', max(max(abs(tril(H9,-2)))), max(max(abs(tril(H9lib,-2)))));
fprintf('\neig(A): '); fprintf('%.6f ', sort(real(eig(A9)))); fprintf('\n');
fprintf('eig(H): '); fprintf('%.6f ', sort(real(eig(H9)))); fprintf('\n\n');

%% Problem 3

fprintf('=== Problem 3: Rayleigh Quotient Iteration (11x11) ===\n\n');

rng(543); m = 11;
Araw = randn(m); A = Araw + Araw';
[V_lib, D_lib] = eig(A);
eigvals_lib = diag(D_lib);

fprintf('Library eigenvalues:\n');
for j=1:m, fprintf('  lambda_%2d = %18.14f\n', j, eigvals_lib(j)); end

fprintf('\n--- RQI runs ---\n');
for s = 1:5
    rng(1000+s);
    v0 = randn(m,1); v0 = v0/norm(v0);
    [lam, v, hist] = rayleigh_quotient_iteration(A, v0);
    [err, idx] = min(abs(eigvals_lib - lam));
    fprintf('Start %d: lambda=%.15f  matches #%d  |err|=%.2e  res=%.2e\n', ...
        s, lam, idx, err, norm(A*v - lam*v));
end

% Verification
fprintf('\n--- Detailed eigenpair verification ---\n');
rng(2024);
v0 = randn(m,1); v0 = v0/norm(v0);
[lam, v, hist] = rayleigh_quotient_iteration(A, v0);
[~, idx] = min(abs(eigvals_lib - lam));
v_lib = V_lib(:,idx);
if dot(v,v_lib)<0, v_lib = -v_lib; end

fprintf('  RQI lambda  = %20.15f\n', lam);
fprintf('  Lib lambda  = %20.15f\n', eigvals_lib(idx));
fprintf('  |difference|= %.2e\n', abs(lam - eigvals_lib(idx)));
fprintf('  ||v_rqi - v_lib|| = %.2e\n', norm(v - v_lib));
fprintf('  ||Av - lam*v||    = %.2e\n\n', norm(A*v - lam*v));

fprintf('  Cubic convergence:\n');
true_lam = eigvals_lib(idx);
for k = 2:length(hist)
    ep = abs(hist(k-1)-true_lam); ec = abs(hist(k)-true_lam);
    if ep > 1e-15
        fprintf('    k=%d: err=%.4e  ratio=%.2f\n', k-1, ec, log10(ec)/log10(ep));
    end
end

%% Functions

function [H, Q] = householder_hessenberg(A)
    m = size(A,1);
    H = A;  Q = eye(m);
    for k = 1:m-2
        x  = H(k+1:m, k);
        s  = sign(x(1)); if s==0, s=1; end
        e1 = zeros(length(x),1); e1(1)=1;
        v  = s*norm(x)*e1 + x;
        v  = v/norm(v);
        H(k+1:m, k:m) = H(k+1:m, k:m) - 2*v*(v'*H(k+1:m, k:m));
        H(1:m, k+1:m) = H(1:m, k+1:m) - 2*(H(1:m, k+1:m)*v)*v';
        Q(:, k+1:m)   = Q(:, k+1:m)   - 2*(Q(:, k+1:m)*v)*v';
    end
    H = triu(H,-1);
end

function [lambda, v, hist] = rayleigh_quotient_iteration(A, v0, maxiter, tol)
    if nargin<3, maxiter=100; end
    if nargin<4, tol=1e-14; end
    m = size(A,1);
    v = v0/norm(v0);
    lambda = v'*A*v;
    hist = lambda;
    for k = 1:maxiter
        w = (A - lambda*eye(m)) \ v;
        v = w/norm(w);
        lambda = v'*A*v;
        hist(end+1) = lambda;
        if norm(A*v - lambda*v) < tol, break; end
    end
end