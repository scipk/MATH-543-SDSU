%% MATH 543 HW 7
% By Parham Khodadi
clear; clc; close all;

%% 6.5.1
fprintf('6.5.1\n');
t0 = tic;

nTrials1    = 1024; % min size for speed
mList       = round(logspace(log10(4), log10(2048), nTrials1));
rho1        = zeros(nTrials1,1);    % pre alloc

for k = 1:nTrials1
    m       = mList(k);
    A       = randn(m); % build mxm matrix with N(0,1) entries
    U       = triu(lu(A));  % Extract U from PA=LU
    rho1(k) = max(abs(U(:))) / max(abs(A(:)));  % eq. 22.2
end
fprintf('   done in %.1f s\n', toc(t0));

f1 = figure('Position',[100 100 800 600],'Color','w');
loglog(mList, rho1, '.', 'MarkerSize', 6);
hold on;
mref = logspace(log10(4), log10(2048), 200);
loglog(mref, sqrt(mref), 'k--', 'LineWidth', 1.25); % rho = sqrt(m)

xlabel('m','FontSize',12);
ylabel('growth factor \rho','FontSize',12);
title('GE w/PP growth factor, N(0,1) entries','FontSize',13);
legend({'\rho','m^{1/2}'},'Location','northwest','FontSize',11);
grid on;
axis tight;
print(f1,'Figures/P6_5_1.eps','-depsc2');

%% 6.5.2
fprintf('6.5.2\n');
t0 = tic;

mSet        = [8 16 32 64];
nTrials2    = 1048576;

f2 = figure('Position',[100 100 800 600],'Color','w'); 
hold on;

for j = 1:numel(mSet)
    m    = mSet(j);
    rhos = zeros(nTrials2,1);
    for k = 1:nTrials2
        A       = randn(m); % build mxm matrix with N(0,1) entries
        U       = triu(lu(A));  % Extract U from PA=LU
        rhos(k) = max(abs(U(:))) / max(abs(A(:)));  % eq. 22.2
    end
    histogram(rhos, 'Normalization','pdf', 'DisplayStyle','stairs', 'DisplayName', sprintf('m = %d', m));
    fprintf('   m = %3d done (%.1f s elapsed)\n', m, toc(t0));
end

xlabel('growth factor  \rho','FontSize',12);
ylabel('probability density','FontSize',12);
title('Probability density of \rho, N(0,1) entries','FontSize',13);
legend('show','FontSize',11); grid on;
print(f2,'Figures/P6_5_2.eps','-depsc2');

%% 6.5.3
fprintf('6.5.3\n');
t0 = tic;

rho3 = zeros(nTrials1,1);
for k = 1:nTrials1
    m         = mList(k);
    A         = rand(m);    % diff from 6.5.1 here uses rand(m) instead of randn(m)
    U         = triu(lu(A));
    rho3(k)   = max(abs(U(:))) / max(abs(A(:)));    % eq. 22.2
end
fprintf('   done in %.1f s\n', toc(t0));

f3 = figure('Position',[100 100 800 600],'Color','w');
loglog(mList, rho3, '.', 'MarkerSize', 6);
hold on;

loglog(mref, sqrt(mref), 'k--', 'LineWidth', 1.25);
xlabel('m','FontSize',12);
ylabel('growth factor  \rho','FontSize',12);
title('GE w/PP growth factor, Uniform[0,1] entries','FontSize',13);
legend({'\rho','m^{1/2}'},'Location','northwest','FontSize',11);
grid on;
axis tight;
print(f3,'Figures/P6_5_3.eps','-depsc2');

%% 6.5.4
fprintf('6.5.4\n');
t0 = tic;

f4 = figure('Position',[100 100 800 600],'Color','w');
hold on;

for j = 1:numel(mSet)
    m    = mSet(j);
    rhos = zeros(nTrials2,1);
    for k = 1:nTrials2
        A       = rand(m);  % Same as 6.5.2 but uses rand(m) instead of randn(m)
        U       = triu(lu(A));
        rhos(k) = max(abs(U(:))) / max(abs(A(:)));  % eq. 22.2
    end
    histogram(rhos, 'Normalization','pdf', 'DisplayStyle','stairs', 'DisplayName', sprintf('m = %d', m));
    fprintf('   m = %3d done (%.1f s elapsed)\n', m, toc(t0));
end

xlabel('growth factor  \rho','FontSize',12);
ylabel('probability density','FontSize',12);
title('Probability density of \rho, Uniform[0,1] entries','FontSize',13);
legend('show','FontSize',11);
grid on;
print(f4,'Figures/P6_5_4.eps','-depsc2');

fprintf('\nAll figures in ./Figures/\n');