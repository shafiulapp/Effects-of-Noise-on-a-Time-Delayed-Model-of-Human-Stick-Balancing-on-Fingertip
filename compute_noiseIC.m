% noise_IC_run.m
%  noise sweep at default dead zone DZ=0.8
% sigma AND DZ set before every call to noise_IC

rng(42)
DZ = 0.8;

% verification
fprintf('compute_noise_IC: DZ = %.2f deg\n', DZ)
sigma = 0;
noise_IC
fprintf('std at sigma=0 = %.3f deg\n\n', std(phiv*180/pi))

% path 1: no noise
sigma = 0;
DZ    = 0.8;
noise_IC

sigma1 = sigma;
phiv1  = phiv;
xv1    = xv;
tv1    = tv;
BT1    = tv(end);

save('noise_IC.mat','sigma1','phiv1','xv1','tv1','BT1','dt')

% path 2: noise near transition
sigma = 0.15;
DZ    = 0.8;
noise_IC

sigma2 = sigma;
phiv2  = phiv;
xv2    = xv;
tv2    = tv;
BT2    = tv(end);

save('noise_IC.mat','sigma2','phiv2','xv2','tv2','BT2','-append')

% Monte Carlo sweep
sigma_list = [0, 0.05, 0.08, 0.10, 0.12, 0.15, 0.155, 0.160, 0.165, 0.170, 0.20, 0.25, 0.30, 0.40, 0.50];
N_MC     = 1000;
BT_mean  = zeros(1,length(sigma_list));
BT_std   = zeros(1,length(sigma_list));
PSD_peak = zeros(1,length(sigma_list));

fprintf('Monte Carlo sweep DZ=0.8: \n')
for k = 1:length(sigma_list)
    BT_runs = zeros(1,N_MC);
    for n = 1:N_MC
        sigma = sigma_list(k);
        DZ    = 0.8;
        noise_IC
        BT_runs(n) = tv(end);
    end
    BT_mean(k) = mean(BT_runs);
    BT_std(k)  = std(BT_runs);
    sigma = sigma_list(k);
    DZ    = 0.8;
    noise_IC
    [pxx,fxx]   = periodogram(phiv,[],[],1/dt);
    [~,idx]     = max(pxx);
    PSD_peak(k) = fxx(idx);
    fprintf('sigma=%.3f | BT=%.1f+/-%.1fs | PSD=%.3fHz\n', sigma_list(k),BT_mean(k),BT_std(k),PSD_peak(k))
end

save('noise_IC.mat','sigma_list','BT_mean','BT_std','PSD_peak','N_MC','-append')

