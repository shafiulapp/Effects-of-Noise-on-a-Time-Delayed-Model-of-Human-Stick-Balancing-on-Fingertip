% deadzone_IC_run.m
% Tests sensitivity to dead zone width
% Calls deadzone_IC.m with different values of DZ


% path 1: default dead zone
DZ = 0.8;
deadzone_IC

DZ1 = DZ;
phiv1 = phiv;
xv1 = xv;
tv1 = tv;
BT1 = tv(end);

save('deadzone_IC.mat','DZ1','phiv1','xv1','tv1','BT1','dt')

% path 2: larger dead zone
DZ = 2.0;
deadzone_IC

DZ2 = DZ;
phiv2 = phiv;
xv2 = xv;
tv2 = tv;
BT2 = tv(end);

save('deadzone_IC.mat','DZ2','phiv2','xv2','tv2','BT2','-append')

% sweep multiple DZ values
DZ_list = [0.5, 0.8, 1.0, 1.2, 1.5, 2.0];
std_list = zeros(1,length(DZ_list));
BT_list  = zeros(1,length(DZ_list));
PSD_peak = zeros(1,length(DZ_list));

for k = 1:length(DZ_list)
    DZ = DZ_list(k);
    deadzone_IC
    std_list(k) = std(phiv*180/pi);
    BT_list(k)  = tv(end);
    [pxx,fxx]   = periodogram(phiv,[],[],1/dt);
    [~,idx]     = max(pxx);
    PSD_peak(k) = fxx(idx);
    fprintf('DZ=%.2f | BT=%.1fs | std=%.3fdeg | PSD=%.3fHz\n', ...
            DZ_list(k),BT_list(k),std_list(k),PSD_peak(k))
end

save('deadzone_IC.mat','DZ_list','std_list','BT_list','PSD_peak','-append')

%% visualize results

clear all, close all
load('deadzone_IC.mat')

% theoretical prediction: std scales with 1.78*Pi
% from free-fall solution theta(tau) = Pi*cosh(omega_n*tau)
% cosh(omega_n*tau) = 1.78 for L=0.56m (Milton et al.)
scale   = std_list(2) / (1.78*DZ_list(2));
std_th  = scale * 1.78 * DZ_list;

% figure 1: overlaid trajectories zoomed in (professor's style)
figure

subplot(2,1,1)
plot(tv1, phiv1*180/pi, 'b.-')
hold on
plot(tv2, phiv2*180/pi, 'r.-')
line([0,max(tv1)], DZ1*[1 1], 'Color','b','LineStyle','--')
line([0,max(tv1)], -DZ1*[1 1], 'Color','b','LineStyle','--')
line([0,max(tv2)], DZ2*[1 1], 'Color','r','LineStyle','--')
line([0,max(tv2)], -DZ2*[1 1], 'Color','r','LineStyle','--')
xlabel('Time (s)')
ylabel('\theta (deg)')
%title('(a) Early stick angle for two dead zone widths')
legend({sprintf('\\Pi = %.1f°', DZ1), ...
        sprintf('\\Pi = %.1f°', DZ2), ...
        'Dead zone boundaries'}, 'Location','best')
set(gca,'FontSize',25)
xlim([0 1])
ylim([-3 4])   % fix scale so blue is visible
grid on

subplot(2,1,2)
plot(tv1, xv1, 'b.-')
hold on
plot(tv2, xv2, 'r.-')
xlabel('Time (s)')
ylabel('x(t) (m)')
legend({sprintf('\\Pi = %.1f°', DZ1), ...
        sprintf('\\Pi = %.1f°', DZ2)}, 'Location','best')
set(gca,'FontSize',25)
xlim([0 1])
grid on

% create two separate figures instead of additional subplots

% figure A: std vs DZ with theoretical line
figure
plot(DZ_list, std_list, 'ko-', 'LineWidth',1.5, 'MarkerSize',7)
hold on
plot(DZ_list, std_th, 'r--', 'LineWidth',1.5)
xlabel('\Pi (deg)')
ylabel('std(\theta) (deg)')
%title('(a) Oscillation amplitude vs dead zone')
legend({'Simulation', 'Theory: 1.78\Pi / \surd2'}, 'Location','best')
set(gca,'FontSize',25)
grid on

% figure B: PSD comparison for the two example runs
figure
[pxx1,fxx1] = periodogram(phiv1,[],[],1/dt);
[pxx2,fxx2] = periodogram(phiv2,[],[],1/dt);
plot(fxx1, pxx1, 'b', 'LineWidth',1.2)
hold on
plot(fxx2, pxx2, 'r', 'LineWidth',1.2)
xlabel('Frequency (Hz)')
ylabel('PSD')
%title('(b) PSD comparison')
legend({sprintf('\\Pi = %.1f°', DZ1), sprintf('\\Pi = %.1f°', DZ2)}, 'Location','best')
set(gca,'FontSize',25)
xlim([0 3])
grid on
