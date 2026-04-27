% noise_DZ2_run.m
% Same Monte Carlo noise sweep as noise_IC_run.m
% but with DZ = 2.0 deg
% Tests whether larger dead zone changes noise tolerance

rng(42)
DZ = 2.0;


fprintf('noise_DZ2_run: DZ = %.2f deg\n', DZ)
sigma = 0;
DZ    = 2.0;
noise_IC
fprintf(' std at sigma=0 = %.3f deg \n\n',std(phiv*180/pi))

% path 1: no noise baseline at DZ=2.0
sigma = 0;
DZ    = 2.0;
noise_IC

sigma1_dz2 = sigma;
phiv1_dz2  = phiv;
xv1_dz2    = xv;
tv1_dz2    = tv;
BT1_dz2    = tv(end);

save('noise_DZ2.mat','sigma1_dz2','phiv1_dz2', ...
     'xv1_dz2','tv1_dz2','BT1_dz2','dt')

% path 2: noise near transition at DZ=2.0
sigma = 0.10;
DZ    = 2.0;
noise_IC

sigma2_dz2 = sigma;
phiv2_dz2  = phiv;
xv2_dz2    = xv;
tv2_dz2    = tv;
BT2_dz2    = tv(end);

save('noise_DZ2.mat','sigma2_dz2','phiv2_dz2', 'xv2_dz2','tv2_dz2','BT2_dz2','-append')

% Monte Carlo sweep :same sigma values for direct comparison
sigma_list_dz2 = [0, 0.05, 0.08, 0.10, 0.12, 0.15, 0.155, 0.160, 0.165, 0.170, 0.20, 0.25, 0.30, 0.40, 0.50];
N_MC         = 1000;
BT_mean_dz2  = zeros(1,length(sigma_list_dz2));
BT_std_dz2   = zeros(1,length(sigma_list_dz2));
PSD_peak_dz2 = zeros(1,length(sigma_list_dz2));

fprintf('Monte Carlo sweep DZ=2.0:\n')
for k = 1:length(sigma_list_dz2)
    BT_runs = zeros(1,N_MC);
    for n = 1:N_MC
        sigma = sigma_list_dz2(k);
        DZ    = 2.0;
        noise_IC
        BT_runs(n) = tv(end);
    end
    BT_mean_dz2(k) = mean(BT_runs);
    BT_std_dz2(k)  = std(BT_runs);
    sigma = sigma_list_dz2(k);
    DZ    = 2.0;
    noise_IC
    [pxx,fxx]       = periodogram(phiv,[],[],1/dt);
    [~,idx]         = max(pxx);
    PSD_peak_dz2(k) = fxx(idx);
    fprintf('sigma=%.3f | BT=%.1f+/-%.1fs | PSD=%.3fHz\n',sigma_list_dz2(k), BT_mean_dz2(k),BT_std_dz2(k), PSD_peak_dz2(k))
end

save('noise_DZ2.mat','sigma_list_dz2','BT_mean_dz2', 'BT_std_dz2','PSD_peak_dz2','N_MC','-append')

