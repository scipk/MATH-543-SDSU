%% MATH 543 - Homework 4
% By Parham Khodadi
clear; clc; close all;

%% Part A: Eigenvalues and Spectral Radius
m_vals = [8, 16, 32, 64, 128];
num_matrices = 100;
figure('Position', [100, 100, 1000, 600]);

for i = 1:length(m_vals)
    m = m_vals(i);
    
    % 1. Preallocate a column vector of the exact required size
    all_eigenvalues = zeros(m * num_matrices, 1); 
    avg_rho = 0;
    
    for k = 1:num_matrices
        A = randn(m,m)/sqrt(m);
        eigs_A = eig(A);
        
        % 2. Calculate the start and end indices for the current batch
        start_idx = (k - 1) * m + 1;
        end_idx = k * m;
        
        % 3. Insert the eigenvalues into the preallocated array
        all_eigenvalues(start_idx:end_idx) = eigs_A; 
        
        avg_rho = avg_rho + max(abs(eigs_A)); % Compute spectral radius
    end
    avg_rho = avg_rho / num_matrices;
    
    % Plot eigenvalues in the complex plane
    subplot(2, 3, i);
    scatter(real(all_eigenvalues), imag(all_eigenvalues), 5, 'filled', 'MarkerFaceAlpha', 0.5);
    title(sprintf('m = %d, avg \\rho = %.3f', m, avg_rho));
    axis equal;
    xlim([-2 2]); ylim([-2 2]);
    xlabel('Real'); ylabel('Imaginary');
    
    % Draw unit circle for reference
    hold on;
    theta = linspace(0, 2*pi, 100);
    plot(cos(theta), sin(theta), 'r--');
    hold off;
end
sgtitle('Part A: Superimposed Eigenvalues of 100 Random Matrices');


%% Part B: Norms vs Spectral Radius
m_vals_large = 2.^(3:9); % Test m values from 8 up to 512
trials = 50; 
avg_norms = zeros(size(m_vals_large));
avg_rhos = zeros(size(m_vals_large));

for i = 1:length(m_vals_large)
    m = m_vals_large(i);
    temp_norm = 0;
    temp_rho = 0;
    for k = 1:trials
        A = randn(m,m)/sqrt(m);
        temp_norm = temp_norm + norm(A, 2);
        temp_rho = temp_rho + max(abs(eig(A)));
    end
    avg_norms(i) = temp_norm / trials;
    avg_rhos(i) = temp_rho / trials;
end

figure;
plot(m_vals_large, avg_norms, '-o', 'LineWidth', 1.5, 'DisplayName', '2-Norm ||A||_2');
hold on;
plot(m_vals_large, avg_rhos, '-s', 'LineWidth', 1.5, 'DisplayName', 'Spectral Radius \rho(A)');
yline(2, 'b--', 'Expected ||A||_2 limit \approx 2', 'HandleVisibility', 'off');
yline(1, 'r--', 'Expected \rho(A) limit \approx 1', 'HandleVisibility', 'off');
xlabel('Matrix Size m');
ylabel('Value');
title('Part B: Behavior of 2-Norm and Spectral Radius as m \rightarrow \infty');
legend('Location', 'best');
grid on;


%% Part C: Condition Numbers and Smallest Singular Value (\sigma_min)
m_vals_C = [8, 16, 32, 64];
num_samples = 2000; % Large sample size for distinct tail probabilities
k_vals = 1:10;
thresholds = 2.^(-k_vals); 

figure;
hold on;
for i = 1:length(m_vals_C)
    m = m_vals_C(i);
    sigma_min_vals = zeros(num_samples, 1);
    
    for j = 1:num_samples
        A = randn(m,m)/sqrt(m);
        s = svd(A); % Returns singular values in descending order
        sigma_min_vals(j) = s(end); % Extract smallest singular value
    end
    
    proportions = zeros(length(thresholds), 1);
    for t = 1:length(thresholds)
        proportions(t) = sum(sigma_min_vals < thresholds(t)) / num_samples;
    end
    
    plot(k_vals, proportions, '-o', 'LineWidth', 1.5, 'DisplayName', sprintf('m = %d', m));
end

set(gca, 'YScale', 'log'); % Log scale to see the tail behavior clearly
xlabel('k (where threshold is 2^{-k})');
ylabel('Proportion of matrices with \sigma_{min} < 2^{-k}');
title('Part C: Tail Distribution of Smallest Singular Value \sigma_{min}');
legend('Location', 'southwest');
grid on;