% noise_DZ2_plot.m
% Loads noise_IC.mat (DZ=0.8) and noise_DZ2.mat (DZ=2.0)
% Produces figures connecting O3 and O4
% Run noise_IC_run.m AND noise_DZ2_run.m first

clear all, close all
load('noise_IC.mat')
load('noise_DZ2.mat')

% stochastic resonance check DZ=0.8
[BT_max_08, idx_max_08] = max(BT_mean);
sigma_opt_08 = sigma_list(idx_max_08);
idx_half_08  = find(BT_mean < BT_mean(1)*0.5, 1,'first');

% stochastic resonance check DZ=2.0
[BT_max_20, idx_max_20] = max(BT_mean_dz2);
sigma_opt_20 = sigma_list_dz2(idx_max_20);
idx_half_20  = find(BT_mean_dz2 < BT_mean_dz2(1)*0.5, 1,'first');

fprintf('DZ=0.8:\n')
if idx_max_08 > 1
    fprintf('SR at sigma*=%.4f\n', sigma_opt_08)
else
    fprintf('No stochastic resonance\n')
end
if ~isempty(idx_half_08)
    fprintf('sigma_c = %.3f\n', sigma_list(idx_half_08))
end

fprintf('\n DZ=2.0:\n')
fprintf(' std at sigma=0 = %.3f deg\n', ...
        std(phiv1_dz2*180/pi))
if idx_max_20 > 1
    fprintf('SR at sigma*=%.4f\n', sigma_opt_20)
    fprintf('BT improvement = %.2f s\n', ...
            BT_mean_dz2(idx_max_20)-BT_mean_dz2(1))
else
    fprintf('No stochastic resonance\n')
end
if ~isempty(idx_half_20)
    fprintf('sigma_c = %.3f\n', sigma_list_dz2(idx_half_20))
end

% BT vs sigma for both DZ values%
figure

errorbar(sigma_list, BT_mean, BT_std, 'bo-', ...
         'LineWidth',1.5, 'MarkerSize',7, ...
         'CapSize',5, 'MarkerFaceColor','b')
hold on
errorbar(sigma_list_dz2, BT_mean_dz2, BT_std_dz2, 'ro-', ...
         'LineWidth',1.5, 'MarkerSize',7, ...
         'CapSize',5, 'MarkerFaceColor','r')

if idx_max_08 > 1
    xline(sigma_opt_08,'b--','LineWidth',1.5)
    text(sigma_opt_08+0.005, BT_max_08*0.95, ...
         sprintf('\\sigma^*=%.3f', sigma_opt_08), ...
         'FontSize',25, 'Color','b')
end
if idx_max_20 > 1
    xline(sigma_opt_20,'r--','LineWidth',1.5)
    text(sigma_opt_20+0.005, BT_max_20*0.90, ...
         sprintf('\\sigma^*=%.3f', sigma_opt_20), ...
         'FontSize',25, 'Color','r')
end
if ~isempty(idx_half_08)
    xline(sigma_list(idx_half_08),'b-.','LineWidth',1.2)
    text(sigma_list(idx_half_08)+0.005, BT_mean(1)*0.4, ...
         sprintf('\\sigma_c=%.3f', sigma_list(idx_half_08)), ...
         'FontSize',25, 'Color','b')
end
if ~isempty(idx_half_20)
    xline(sigma_list_dz2(idx_half_20),'r-.','LineWidth',1.2)
    text(sigma_list_dz2(idx_half_20)+0.005, BT_mean_dz2(1)*0.3, ...
         sprintf('\\sigma_c=%.3f', sigma_list_dz2(idx_half_20)), ...
         'FontSize',25, 'Color','r')
end

xlabel('\sigma  (m s^{-3/2})')
ylabel('Mean balance time (s)')
%title(sprintf(['Effect of noise: \\Pi = 0.8° vs \\Pi = 2.0°', '  (N = %d)'], N_MC))
legend({'\Pi = 0.8°', '\Pi = 2.0°'}, 'Location','best')
set(gca,'FontSize',25)
grid on


%figure 2: trajectory and PSD comparison%
rng(42)
sigma_compare = 0.15;

DZ    = 0.8;
sigma = sigma_compare;
noise_IC
phiv_08 = phiv;
tv_08   = tv;

DZ    = 2.0;
sigma = sigma_compare;
noise_IC
phiv_20 = phiv;
tv_20   = tv;

N = min(length(tv_08), length(tv_20));

figure

subplot(2,1,1)
plot(tv_08(1:N), phiv_08(1:N)*180/pi, 'b', 'LineWidth',1.2)
hold on
plot(tv_20(1:N), phiv_20(1:N)*180/pi, 'r', 'LineWidth',1.2)
line([0,tv_08(min(N,end))],  0.8*[1 1],'Color','b','LineStyle','--')
line([0,tv_08(min(N,end))], -0.8*[1 1],'Color','b','LineStyle','--')
line([0,tv_20(min(N,end))],  2.0*[1 1],'Color','r','LineStyle','--')
line([0,tv_20(min(N,end))], -2.0*[1 1],'Color','r','LineStyle','--')
xlabel('Time (s)')
ylabel('\theta (deg)')
%title(sprintf('Trajectories at \\sigma = %.3f: \\Pi=0.8° vs \\Pi=2.0°',  sigma_compare))
legend({'\Pi = 0.8°', '\Pi = 2.0°', ...
        'Dead zone boundaries'}, 'Location','best')
set(gca,'FontSize',25)
xlim([0 min(30, min(tv_08(end),tv_20(end)))])
grid on

subplot(2,1,2)
[pxx_08,fxx_08] = periodogram(phiv_08,[],[],1/dt);
[pxx_20,fxx_20] = periodogram(phiv_20,[],[],1/dt);
plot(fxx_08, pxx_08, 'b', 'LineWidth',1.2)
hold on
plot(fxx_20, pxx_20, 'r', 'LineWidth',1.2)
xlabel('Frequency (Hz)')
ylabel('PSD')
%title(sprintf('PSD at \\sigma = %.3f', sigma_compare))
legend({'\Pi = 0.8°', '\Pi = 2.0°'}, 'Location','best')
set(gca,'FontSize',25)
xlim([0 3])
grid on


fprintf('%-8s %-12s %-12s %-12s\n', ...
        'Pi(deg)','std(sig=0)','sigma_c','SR')

% DZ=0.8
if idx_max_08 > 1
    sr_08 = 'Yes';
else
    sr_08 = 'No';
end
if ~isempty(idx_half_08)
    sc_08 = sigma_list(idx_half_08);
else
    sc_08 = NaN;
end
fprintf('%-8.1f %-12.3f %-12.3f %-12s\n', ...
        0.8, std(phiv1*180/pi), sc_08, sr_08)

% DZ=2.0
if idx_max_20 > 1
    sr_20 = 'Yes';
else
    sr_20 = 'No';
end
if ~isempty(idx_half_20)
    sc_20 = sigma_list_dz2(idx_half_20);
else
    sc_20 = NaN;
end
fprintf('%-8.1f %-12.3f %-12.3f %-12s\n', ...
        2.0, std(phiv1_dz2*180/pi), sc_20, sr_20)


fprintf(' sigma_c reduced from %.3f to %.3f\n', ...
        sc_08, sc_20)
fprintf('Noise tolerance reduced by %.1f%% for larger dead zone\n', ...
        100*(sc_08-sc_20)/sc_08)