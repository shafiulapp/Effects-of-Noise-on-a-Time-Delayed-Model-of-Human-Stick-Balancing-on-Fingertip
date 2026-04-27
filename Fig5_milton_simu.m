% pendulum-cart model with feedback delay and predcitor feedback (PF) control
% time-domain simulation obtained using the zeroth-order semidiscretization
% P1, D1: control gains for the angular position of the stick
% P2, D2: control gains for the location of the pivot point (cart)
% DV(1): deadzone for angular position of the pendulum
% DV(2): deadzone for position of the pivot point
% DV(3): deadzone for angular velocity of the pendulum
% DV(4): deadzone for velocity of the pivot point


clear
tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameters to give
L = 0.56;       % [m] stick length
DZ = 0.8;   % dead zone [deg]
m = 0.0247*L;     % mass of stick [kg]
m0 = 1.2;        % mass of the cart (hand) [kg]
c = L/2;    % location of the cente of gravity [m]
dt = 0.01;   % discretization time step [s]
r = 23;       % discrete delay, tau = r*dt
g = 9.81;    % grav. acc. [m/s2]
tmax = 240;   % [s]
acc_max = 50;  % max acceleration of the fingertip
jerk_max = 600;  % max jerk of the fingertip m/s^2
x_check_lim = 0.335; % max fingertip displacement [m]
phi_check_lim = 20;  % max stick angle [deg]
wx = 1200;   % weight for x in the cost function

% limits of control gain: s-starting, e-end, st-step
P1 = 55;
D1 = 14;
P2 = 10;
D2 = 20;

% deadzones DV=[angle, position, ang.velocity, velocity]
% in [deg, m, deg/s, m/s]
% no deadzone for the efferent copies
DVdeg = [DZ; 0; 0; 0];
% DV in rad
DV = [DVdeg(1)*pi/180; DVdeg(2); DVdeg(3)*pi/180; DVdeg(4)];

% initial conditons
phi0 = 0.6*pi/180; % [rad]
dphi0 = 0;   % [rad/s]
x0 = 0;      % [m]
dx0 = 0;     % [m/s]

% end of parameters to give
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parameters, system matrices for semidiscreization
Ic = m*L^2/12;
M = [[m*c, m+m0]; [Ic+m*c^2, m*c]];
K = [[0, 0]; [-m*g*c, 0]];
invM = inv(M);
B = [0; 0; invM*[1; 0]];
A = [[zeros(2), eye(2)]; [-invM*K, zeros(2)]];
P = expm(A*dt);
A_int = @(s)expm(A*(dt-s));
R = quadv(A_int,0,dt);
Phi = zeros(4+r) + diag(ones(r+3,1),-1);
Phi(1:4,1:4) = P;
Phi(1:4,4+r) = R*B;

% max control force
m0corr = 1/([0,1]*inv(M)*[1;0]);


% scale for colorbar
iimax = 10;
for ii = 1:iimax
    im = ii;
    cmp(ii,:) = 1-[im/iimax, im/iimax, im/iimax];
end
colormap(cmp)
% scale for colorbar

D = [P1, P2, D1, D2];
S = D*P^r;
Phi(5,1:4) = S;
for j = 1:r
    V(j) = D*P^(j-1)*R*B;
    Phi(5,4+j) = V(j);
end

% simulation for iicmax different initial conditions
iic = 0;
tiic = 0;
X = zeros(r+4,1);
X(1:4) = [phi0; x0; dphi0; dx0];
i = 0;
x_check = 0;
phi_check = 0;
u_prev = 0;
v_prev = 0;
while (x_check < x_check_lim)*(phi_check < phi_check_lim)*(i*dt<tmax)
    i = i+1;
    X0 = X;
    DX0 = X0(1:4).*(abs(X0(1:4))>DV); % set to zero if within deadzone
    X(1:4) = P*X0(1:4) + R*B*X0(r+4);
    u = S*DX0(1:4) + V*X0(5:(4+r));
    if abs(u)/m0corr < acc_max
        ua = u;
    else
        ua = sign(u)*acc_max*m0corr;
    end
    if abs(ua-u_prev)/m0corr/dt < jerk_max
        uaj = ua;
    else
        uaj = u_prev + sign(ua-u_prev)*jerk_max*m0corr*dt;
    end
    X(5) = uaj;
    X(6:(r+4)) = X0(5:(r+3));
    x_check = abs(X(2));
    phi_check = abs(X(1)*180/pi);
    tv(i) = i*dt;
    phiv(i)=X(1);
    xv(i)=X(2);
    u_prev = uaj;
    v_prev = X(4);
end

subplot(2,1,1)
plot(tv, phiv*180/pi, 'k', 'LineWidth', 1.0)
xlabel('t (s)')
ylabel('\theta (deg)')
%title('Simulation A --- stick angle \theta(t)')
set(gca,'FontSize',13)
xlim([0 60])
grid on

subplot(2,1,2)
[pt,ft] = periodogram(phiv,[],[],1/dt);
plot(ft, pt, 'k', 'LineWidth', 1.0)
xlabel('Frequency (Hz)')
ylabel('PSD')
%title('Simulation A --- power spectral density')
set(gca,'FontSize',13)
xlim([0 3])
grid on

% print key values to compare with Milton et al.
[~,idx] = max(pt);
fprintf('std(theta) = %.3f deg\n', std(phiv*180/pi))
fprintf('PSD peak   = %.3f Hz\n', ft(idx))
fprintf('BT         = %.0f s\n',  tv(end))
fprintf('toc: %.1f s\n', toc)

exportgraphics(gcf,'fig_reproduction.png','Resolution',300)