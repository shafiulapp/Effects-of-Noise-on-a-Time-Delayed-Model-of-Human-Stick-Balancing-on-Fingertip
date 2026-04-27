% noise_IC.m
% Base script for noise effects study
% sigma AND DZ must both be defined before calling this script
% sigma: noise intensity [m s^{-3/2}]
% DZ:    dead zone width [deg]

L = 0.56;
m = 0.0247*L;
m0 = 1.2;
c = L/2;
dt = 0.01;
r = 23;
g = 9.81;
tmax = 240;
acc_max = 50;
jerk_max = 600;
x_check_lim = 0.335;
phi_check_lim = 20;
P1 = 55; D1 = 14; P2 = 10; D2 = 20;

DVdeg = [DZ; 0; 0; 0];
DV = [DVdeg(1)*pi/180; DVdeg(2); DVdeg(3)*pi/180; DVdeg(4)];

phi0  = 0.6*pi/180;
dphi0 = 0;
x0    = 0;
dx0   = 0;

Ic   = m*L^2/12;
M    = [m*c, m+m0; Ic+m*c^2, m*c];
K    = [0, 0; -m*g*c, 0];
invM = inv(M);
B    = [0; 0; invM*[1; 0]];
A    = [zeros(2), eye(2); -invM*K, zeros(2)];
P    = expm(A*dt);
A_int = @(s)expm(A*(dt-s));
R    = quadv(A_int,0,dt);
m0corr = 1/([0,1]*invM*[1;0]);
D = [P1,P2,D1,D2];
S = D*P^r;
for j = 1:r
    V(j) = D*P^(j-1)*R*B;
end

X = zeros(r+4,1);
X(1:4) = [phi0; x0; dphi0; dx0];
i = 0; x_check = 0; phi_check = 0; u_prev = 0;
tv = []; phiv = []; xv = [];

while (x_check < x_check_lim)*(phi_check < phi_check_lim)*(i*dt < tmax)
    i  = i+1;
    X0 = X;
    DX0    = X0(1:4).*(abs(X0(1:4))>DV);
    X(1:4) = P*X0(1:4) + R*B*X0(r+4);
    X(4)   = X(4) + sigma*sqrt(dt)*randn;
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
    X(5)       = uaj;
    X(6:(r+4)) = X0(5:(r+3));
    x_check    = abs(X(2));
    phi_check  = abs(X(1)*180/pi);
    tv(i)   = i*dt;
    phiv(i) = X(1);
    xv(i)   = X(2);
    u_prev  = uaj;
end