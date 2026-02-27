%% MATH 543 - Homework 2
% By Parham Khodadi
clear; clc; close all;

%% TB-4.1: Singular Value Decompositions
fprintf("TB-4.1\n\n");

matrices = { ...
    [3,0;0,-2], ...
    [2,0;0,3], ...
    [0,2;0,0;0,0], ...
    [1,1;0,0], ...
    [1,1;1,1] ...
};
labels = {'(a)', '(b)', '(c)', '(d)', '(e)'};

for k = 1:numel(matrices)
    A = matrices{k};
    [U,S,V] = svd(A);
    fprintf('SVD of Matrix %s:\n', labels{k});
    fprintf('U:\n'); disp(U);
    fprintf('S:\n'); disp(S);
    fprintf('V:\n'); disp(V);
    fprintf('\n');
end

fprintf("-----------------\n");

%% TB-4.3: Unit Circle / Ellipse Plots for 2x2 Matrices

% Matrix (3.7) from the textbook
A_37 = [1,2;0,2];

matrices_2x2 = { ...
    A_37, ...
    [3,0;0,-2], ...
    [2,0;0,3], ...
    [1,1;0,0], ...
    [1,1;1,1] ...
};
labels_2x2 = {'3.7', 'a', 'b', 'd', 'e'};

for k = 1:numel(matrices_2x2)
    A = matrices_2x2{k};
    lbl = labels_2x2{k};

    [U,S,V] = svd(A);
    sigma1 = S(1,1);
    sigma2 = S(2,2);
    v1 = V(:,1);
    v2 = V(:,2);
    u1 = U(:,1);
    u2 = U(:,2);

    fprintf('\n------ TB-4.3: Matrix (%s) ------\n', lbl);
    fprintf('Sigma1 = %.4f\nSigma2 = %.4f\n', sigma1, sigma2);
    fprintf('v1 = [%.6g; %.6g], v2 = [%.6g; %.6g]\n', v1(1), v1(2), v2(1), v2(2));
    fprintf('u1 = [%.6g; %.6g], u2 = [%.6g; %.6g]\n', u1(1), u1(2), u2(1), u2(2));

    plot_svd(A, sigma1, sigma2, v1, v2, u1, u2, lbl);
end

%% Plot Function
function plot_svd(A, sigma1, sigma2, v1, v2, u1, u2, lbl)
    t = linspace(0, 2*pi, 31415);
    X = [cos(t); sin(t)];
    Y = A * X;

    % Unit circle with right singular vectors
    figure;
    plot(X(1,:), X(2,:), 'k', 'LineWidth', 1.5); hold on; grid on;
    quiver(0, 0, v1(1), v1(2), 0, 'LineWidth', 2, 'MaxHeadSize', 0.25);
    quiver(0, 0, v2(1), v2(2), 0, 'LineWidth', 2, 'MaxHeadSize', 0.25);
    axis equal;
    xlabel('x'); ylabel('y');
    title(sprintf('Unit circle with right singular vectors (%s)', lbl));
    legend('unit circle', '$v_1$', '$v_2$', ...
        'Location', 'best', 'Interpreter', 'latex');
    exportgraphics(gcf, sprintf('%s_unit_circle_V.eps', lbl), 'ContentType', 'vector');

    % Ellipse with left singular vectors
    p1 = sigma1 * u1;
    p2 = sigma2 * u2;

    figure;
    plot(Y(1,:), Y(2,:), 'k', 'LineWidth', 1.5); hold on; grid on;
    quiver(0, 0, p1(1), p1(2), 0, 'LineWidth', 2, 'MaxHeadSize', 0.25);
    quiver(0, 0, p2(1), p2(2), 0, 'LineWidth', 2, 'MaxHeadSize', 0.25);
    axis equal;
    xlabel('x'); ylabel('y');
    title(sprintf('Ellipse with left singular vectors (%s)', lbl));
    legend('ellipse', '$\sigma_1 u_1$', '$\sigma_2 u_2$', ...
        'Location', 'best', 'Interpreter', 'latex');
    exportgraphics(gcf, sprintf('%s_ellipse_U.eps', lbl), 'ContentType', 'vector');
end