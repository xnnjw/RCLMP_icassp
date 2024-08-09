clc
clear all
close all
%
addpath functions
addpath utils
%% General Setting
%
MC_iters = 1;   % number of MC iterations
P       = 1;    % total scale (power) of X
N       = 1024; % number of MTD-s
fading.type='uniform';
fading.lower_limit = -15;
fading.upper_limit = 0;
pilot ='bernoulli';

%% 1/2 - M
%
L     = 32;  % number of pilots
Mlist = 8:8:32;
Klist = [8,16,24,32];

Pmd = zeros(length(Mlist),3,length(Klist));
Per = zeros(length(Mlist),3,length(Klist));

t_allstart = tic;
for k = 1:length(Klist)
    rng('default');
    K_a = Klist(k);
    fprintf('\n---- K=%3d',K_a)
    for m = 1:length(Mlist)
        M   = Mlist(m); % choose the number of antennas
        fprintf('\n- M=%3d ',M);
        [Pmd(m,:,k),Per(m,:,k),~,~] = activityDetectionPE(L,N,K_a,M,P,MC_iters,'fading',fading,'pilot',pilot);
        fprintf(' in %.2f mins ~~~ \n', toc(t_allstart) / 60);
    end
end
elapsed_time = toc(t_allstart)/3600;
fprintf(' === All done in %.2f mins, or %.3f hours. === \n', elapsed_time*60, elapsed_time);
fprintf(' ============ M is done, now L ============ \n');

Tensor2Plot_pmd_M = Pmd;
%% 2/2 - L
M       = 32;  % number of antennas
Llist   = 8:8:32;

Pmd = zeros(length(Llist),3,length(Klist));
Per = zeros(length(Llist),3,length(Klist));

t_allstart = tic;
for k = 1:length(Klist)
    rng('default');
    K_a = Klist(k);
    fprintf('\n---- K=%3d',K_a)
    for l = 1:length(Llist)
        L = Llist(l); % choose the number of antennas
        fprintf('\n- L=%3d ',L);
        [Pmd(l,:,k),Per(l,:,k),~,~] = activityDetectionPE(L,N,K_a,M,P,MC_iters,'fading',fading,'pilot',pilot);
        fprintf(' in %.2f mins ~~~ \n', toc(t_allstart) / 60);
    end
end
elapsed_time = toc(t_allstart)/3600;
fprintf(' === All done in %.2f mins, or %.3f hours. === \n', elapsed_time*60, elapsed_time);

Tensor2Plot_pmd_L = Pmd;

%% Plot
figure;
for k = 1:length(Klist)
    subplot(2, 4, k);
    plot(Mlist, Tensor2Plot_pmd_M(:, :, k), 'LineWidth', 1.5);
    set(gca, 'YScale', 'log');
    xlabel('M');
    if k == 1
        ylabel('PMD');
    end
    title(['K=', num2str(Klist(k))]);
    grid on;
    
    subplot(2, 4, k + length(Klist));
    plot(Llist, Tensor2Plot_pmd_L(:, :, k), 'LineWidth', 1.5);
    set(gca, 'YScale', 'log');
    xlabel('L');
    if k == 1
        ylabel('PMD');
    end
    grid on;
end

legend('CWO', 'CL-OMP', 'RCL-MP', 'Location', 'southwest');