% Tests sensitivity to initial conditions
% Calls sensitivity_data.m with two different values of phi0


%  baseline
phi0 = 0.6*pi/180;
sensitivity_data

IC1 = phi0;
phiv1 = phiv;
xv1 = xv;
tv1 = tv;

save('sensitivity_data.mat','IC1','phiv1','xv1','tv1')

%  perturbed by delta = 1e-3 degrees
delta = (1e-3)*pi/180;
phi0 = 0.6*pi/180 + delta;
sensitivity_data

IC2 = phi0;
phiv2 = phiv;
xv2 = xv;
tv2 = tv;

save('sensitivity_data.mat','IC2','phiv2','xv2','tv2','-append')

% compute and visualize

clear all, close all
load('sensitivity_data.mat')
dt = 0.01;


N = min(length(tv1),length(tv2));
tv = tv1(1:N);
sep = abs(phiv1(1:N) - phiv2(1:N));
sep = max(sep,1e-20);

% velocities for phase portrait
dphiv1 = diff(phiv1)/dt;
dphiv2 = diff(phiv2)/dt;


env_vals = [];
env_time = [];
for k = 2:length(sep)-1
    if sep(k) > sep(k-1) && sep(k) > sep(k+1)
        env_vals(end+1) = log10(sep(k));
        env_time(end+1) = tv(k);
    end
end

% fit lambda on growth phase
fit_start = 17;
fit_end = 35;
idx_fit = env_time >= fit_start & ...
          env_time <= fit_end & ...
          10.^env_vals > 1e-8;
p = polyfit(env_time(idx_fit), ...
            env_vals(idx_fit)*log(10), 1);
lambda = p(1);

t_fit_line = linspace(min(env_time(idx_fit)), ...
                      max(env_time(idx_fit)), 200);
fit_line = exp(polyval(p, t_fit_line));

fprintf('lambda = %.4f s^{-1}\n', lambda)
if lambda > 0
    fprintf('lambda > 0: microchaos confirmed\n')
end


% stick angle trajectories
figure
plot(tv1, phiv1*180/pi, 'b','LineWidth', 1.5)
hold on
plot(tv2, phiv2*180/pi, 'r--','LineWidth', 1.5)
xlabel('Time (s)')
ylabel('\theta (deg)')
%title('(a) Stick angle: \phi_0 vs \phi_0+\delta')
legend({'\phi_0 = 0.6°','\phi_0+\delta,  \delta=10^{-3}^°'},'Location','best')
set(gca,'FontSize',25)
xlim([0 40])
grid on

dphiv1 = gradient(phiv1, tv1);
dphiv2 = gradient(phiv2, tv2);

windows = [0 8; 15 23; 30 38];

figure
for i = 1:3
    t1 = windows(i,1);
    t2 = windows(i,2);

    idx1 = tv1 >= t1 & tv1 <= t2;
    idx2 = tv2 >= t1 & tv2 <= t2;

    subplot(1,3,i)
    plot(phiv1(idx1)*180/pi, dphiv1(idx1)*180/pi, 'b', 'LineWidth', 1.5)
    hold on
    plot(phiv2(idx2)*180/pi, dphiv2(idx2)*180/pi, 'r--', 'LineWidth', 1.5)
    xlabel('\theta (deg)')
    ylabel('d\theta/dt (deg/s)')
    title(sprintf('%d-%d s', t1, t2))
    grid on
    set(gca,'FontSize',16)
end
legend('\phi_0','\phi_0+\delta')



