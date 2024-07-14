% LMS convergence simulation and model
% JIN CHIY
% version 2024-07

 
clear, close all

%% =============== Simulation Parameters ===============
% MC realization number
RLZ =50;
% step size
mu = 0.01; %0.005
% Data length
N = 2e3;


%% =============== System Parameters ===============
% System coefficient
wopt = [0.8 0.6 0.4 0.3 0.2 0.1 0 -0.2 -0.1 0.05]';
% System order
L = length(wopt);
% Input correlation
rho = 0.5;
% noise variance
varn = 0.001;


%% ===============  MC filtering process ===============
% initialization of error vectors
e2 = zeros(N,1);
e2av = e2;
exe2 = e2;
exe2av = e2;
vn = zeros(L,N);
vnav = vn;

for r  = 1 : RLZ
    % Data generation
    xw = sqrt(1-rho^2)*randn(N,1);    % white sitimulator to generate a normalized input
    x = filter(1,[1, -rho], xw);      % input signal
    y = filter(wopt,1,x);             % output signal
    z = sqrt(varn)*randn(N,1);        % observation noise
    d = y + z;                        % noisy output (reference)
    
    % filtering
    w = zeros(L,1);                   % Weight initialization
    for n = L : N
       
       xn = x(n: -1: n-L+1);          % instant input
       en = d(n) - w'*xn;             % estimation error
       
       exen = en - z(n);              % excess error
       vn(:,n) = w - wopt;            % weight erro
       
       
       w = w + mu * xn * en;          % filter update
       e2(n) = en^2;
       exe2(n) = exen^2;
    end
    e2av = e2av + e2;
    exe2av = exe2av + exe2;
    vnav = vnav + vn;
end

e2av = e2av/RLZ;
exe2av = exe2av/RLZ;
vnav = vnav/RLZ;

%% ===============  Theoretical model ===============
% Correlation matrix
R = toeplitz(rho.^[0: L-1]);
% initialization
vt = zeros(L,N);
vt(:,L) = -wopt;
K = wopt*wopt';
et = zeros(N,1);
for n = L+1 : N
    
    % First-order model
    vt(:,n) = (eye(L)-mu*R)*vt(:,n-1);
    
    % Second-order model
    % Either the following two update relation for K can be used:
%    K = K+ mu^2*(trace(R*K)*R+2*R*K*R) - 2*mu*R*K + mu^2*varn*R;
     K = (eye(L)-mu*R)*K*(eye(L)-mu*R) + mu^2*varn*R;
    
    exe2t(n) =  trace(K*R);
    e2t(n) = varn + exe2t(n); 
end


%% ===============  Illustration ===============
% Weight evolution
figure, hold on, grid on
for i = 1 : L
    plot(vnav(i,L:end) + wopt(i),'linewidth',2);
    plot(vt(i,L:end) + wopt(i),'r','linewidth',2);
end
xlabel('iteration', 'fontsize', 14,'interpreter', 'Latex')
ylabel('Weights','fontsize', 14,'interpreter', 'Latex')

% Mean square error
figure, semilogy(e2av(L:end))
grid on, hold on, semilogy(e2t,'r', 'linewidth',2)
xlabel('iteration', 'fontsize', 14,'interpreter', 'Latex')
ylabel('MSE','fontsize', 14,'interpreter', 'Latex')

% Excess mean square error
figure, semilogy(exe2av(L:end))
grid on, hold on, semilogy(exe2t,'r', 'linewidth',2)
xlabel('iteration', 'fontsize', 14,'interpreter', 'Latex')
ylabel('EMSE','fontsize', 14,'interpreter', 'Latex')