% noise_IC_plot.m
% Loads noise_IC.mat and produces figures for objective 4
% Run compute_noiseIC.m first

clear all, close all
load('noise_IC.mat')
DZ = 0.8;

[BT_max, idx_max] = max(BT_mean);
sigma_opt = sigma_list(idx_max);

if idx_max > 1
    fprintf('Stochastic resonance at sigma* = %.4f\n', sigma_opt)
    fprintf('Mean BT at sigma=0  = %.2f s\n', BT_mean(1))
    fprintf('Mean BT at sigma*   = %.2f s\n', BT_mean(idx_max))
    fprintf('Improvement         = %.2f s\n', ...
            BT_mean(idx_max)-BT_mean(1))
    sigma_opt_plot = sigma_opt;
else
    fprintf('No stochastic resonance: BT decreases monotonically\n')
    idx_trans = find(BT_mean < 200, 1, 'first');
    if ~isempty(idx_trans)
        sigma_opt_plot = sigma_list(idx_trans);
    else
        sigma_opt_plot = sigma_list(round(length(sigma_list)/2));
    end
end

idx_half = find(BT_mean < BT_mean(1)*0.5, 1, 'first');
if ~isempty(idx_half)
    sigma_critical = sigma_list(idx_half);
    fprintf('Critical sigma: %.3f\n', sigma_critical)
end

%two trajectories%
figure

subplot(2,1,1)
plot(tv1, phiv1*180/pi, 'b.-')
hold on
plot(tv2, phiv2*180/pi, 'r.-')
line([0,max(tv1)],  0.8*[1 1], 'Color','k','LineStyle','--')
line([0,max(tv1)], -0.8*[1 1], 'Color','k','LineStyle','--')
xlabel('Time (s)')
ylabel('\theta (deg)')
%title('(a) Stick angle: no noise vs with noise')
legend({'\sigma = 0  (no noise)', ...
        sprintf('\\sigma = %.3f', sigma2)}, 'Location','best')
set(gca,'FontSize',25)
xlim([0 1])
grid on

subplot(2,1,2)
plot(tv1, xv1, 'b.-')
hold on
plot(tv2, xv2, 'r.-')
xlabel('Time (s)')
ylabel('x(t) (m)')
legend({'\sigma = 0', ...
        sprintf('\\sigma = %.3f', sigma2)}, 'Location','best')
set(gca,'FontSize',25)
xlim([0 1])
grid on

% mean BT vs sigma%
figure

errorbar(sigma_list, BT_mean, BT_std, 'ko-', ...
         'LineWidth',1.5, 'MarkerSize',7, ...
         'CapSize',6, 'MarkerFaceColor','k')
hold on
if idx_max > 1
    xline(sigma_opt,'r--','LineWidth',1.5)
    text(sigma_opt+0.005, BT_max*0.95, ...
         sprintf('\\sigma^* = %.3f', sigma_opt), ...
         'FontSize',25, 'Color','r')
end
if ~isempty(idx_half)
    xline(sigma_critical,'b--','LineWidth',1.5)
    text(sigma_critical+0.005, BT_mean(1)*0.55, ...
         sprintf('\\sigma_c = %.3f', sigma_critical), ...
         'FontSize',25, 'Color','b')
end
xlabel('\sigma  (m s^{-3/2})')
ylabel('Mean balance time (s)')
%title(sprintf(['Mean balance time vs noise intensity','  (N = %d runs per \\sigma)'], N_MC))
set(gca,'FontSize',25)
grid on

% three trajectories %
rng(42)

sigma_large = sigma_list(end);
sigma_plot  = [0, sigma_opt_plot, sigma_large];
col_list    = {'b','g','r'};

if idx_max > 1
    mid_label = sprintf('\\sigma = %.3f  (optimal \\sigma^*)', ...
                        sigma_opt_plot);
else
    mid_label = sprintf('\\sigma = %.3f  (transition region)', ...
                        sigma_opt_plot);
end

%title_list = {'\sigma = 0  (no noise)',   mid_label,   sprintf('\\sigma = %.3f  (large noise)', sigma_large)};

figure
for k = 1:3
    sigma = sigma_plot(k);
    DZ    = 0.8;
    noise_IC
    subplot(3,1,k)
    plot(tv, phiv*180/pi, col_list{k}, 'LineWidth',1.2)
    xlabel('Time (s)')
    ylabel('\theta (deg)')
  %  title(title_list{k})
    set(gca,'FontSize',25)
    xlim([0 5])
    grid on
end

% PSD at three noise levels%
rng(42)

sigma_before = 0;
sigma_during = 0.20;
sigma_after  = 0.30;
min_BT       = 30;

figure

sigma = sigma_before;
DZ    = 0.8;
noise_IC
[pxx1,fxx1] = periodogram(phiv,[],[],1/dt);

subplot(3,1,1)
plot(fxx1, pxx1, 'b', 'LineWidth',1.2)
ylabel('PSD')
%title(sprintf('(a) \\sigma = 0  (before transition),  BT = %.0f s', ... tv(end)))
set(gca,'FontSize',25)
xlim([0 3])
grid on

sigma   = sigma_during;
DZ      = 0.8;
found   = false;
attempt = 0;
while ~found && attempt < 100
    attempt = attempt + 1;
    noise_IC
    if tv(end) >= min_BT
        found  = true;
        phiv_b = phiv;
        tv_b   = tv;
    end
end

if found
    [pxx2,fxx2] = periodogram(phiv_b,[],[],1/dt);
    [~,idx_d] = min(abs(sigma_list - sigma_during));
    subplot(3,1,2)
    plot(fxx2, pxx2, 'g', 'LineWidth',1.2)
    ylabel('PSD')
   % title(sprintf(['(b) \\sigma = %.3f  (transition),'  '  BT = %.0f s  (mean = %.1f s)'],      sigma_during, tv_b(end), BT_mean(idx_d)))
    set(gca,'FontSize',25)
    xlim([0 3])
    grid on
end

sigma   = sigma_after;
DZ      = 0.8;
found   = false;
attempt = 0;
while ~found && attempt < 200
    attempt = attempt + 1;
    noise_IC
    if tv(end) >= min_BT
        found  = true;
        phiv_c = phiv;
        tv_c   = tv;
    end
end

if found
    [pxx3,fxx3] = periodogram(phiv_c,[],[],1/dt);
    [~,idx_a] = min(abs(sigma_list - sigma_after));
    subplot(3,1,3)
    plot(fxx3, pxx3, 'r', 'LineWidth',1.2)
    xlabel('Frequency (Hz)')
    ylabel('PSD')
  %  title(sprintf(['(c) \\sigma = %.3f  (after transition),' '  BT = %.0f s  (mean = %.1f s)'], sigma_after, tv_c(end), BT_mean(idx_a)))
    set(gca,'FontSize',25)
    xlim([0 3])
    grid on
end


fprintf('DZ            : %.2f deg\n', DZ)
fprintf('sigma tested  : %d values\n', length(sigma_list))
fprintf('N MC runs     : %d per sigma\n', N_MC)
fprintf('BT at sigma=0 : %.2f +/- %.2f s\n', BT_mean(1), BT_std(1))
if idx_max > 1
    fprintf('sigma*        : %.4f\n', sigma_opt)
    fprintf('Result        : STOCHASTIC RESONANCE CONFIRMED\n')
else
    fprintf('Result        : no stochastic resonance\n')
    if ~isempty(idx_half)
        fprintf('sigma_c       : %.3f\n', sigma_critical)
    end
end
