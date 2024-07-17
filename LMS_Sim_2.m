% LMS convergence simulation and model
% JIN CHIY
% version 2024-07

clear, close all

%% =============== Simulation Parameters ===============
RLZ = 50;          % MC realization number
mu = 0.01;         % Step size
N = 2000;          % Data length

%% =============== System Parameters ===============
wopt = [0.8 0.6 0.4 0.3 0.2 0.1 0 -0.2 -0.1 0.05]'; % System coefficient
L = length(wopt);    % System order
rho = 0.5;           % Input correlation
varn = 0.001;        % Noise variance

%% ===============  MC filtering process ===============
[e2av, exe2av, vnav] = mcFilteringProcess(RLZ, N, L, wopt, mu, rho, varn);

%% ===============  Theoretical model ===============
[e2t, exe2t, vt] = theoreticalModel(N, L, mu, wopt, rho, varn);

%% ===============  Illustration ===============
plotResults(vnav, vt, wopt, e2av, e2t, exe2av, exe2t, L);

%% Functions
function [e2av, exe2av, vnav] = mcFilteringProcess(RLZ, N, L, wopt, mu, rho, varn)
    e2 = zeros(N, 1);
    e2av = e2;
    exe2 = e2;
    exe2av = e2;
    vn = zeros(L, N);
    vnav = vn;

    for r = 1:RLZ
        [x, d, z] = generateData(N, rho, wopt, varn);
        [e2, exe2, vn] = lmsFiltering(N, L, mu, x, d, wopt, z, vn);
        e2av = e2av + e2;
        exe2av = exe2av + exe2;
        vnav = vnav + vn;
    end

    e2av = e2av / RLZ;
    exe2av = exe2av / RLZ;
    vnav = vnav / RLZ;
end

function [x, d, z] = generateData(N, rho, wopt, varn)
    xw = sqrt(1 - rho^2) * randn(N, 1); % White noise generator
    x = filter(1, [1, -rho], xw);       % Input signal
    y = filter(wopt, 1, x);             % Output signal
    z = sqrt(varn) * randn(N, 1);       % Observation noise
    d = y + z;                          % Noisy output (reference)
end

function [e2, exe2, vn] = lmsFiltering(N, L, mu, x, d, wopt, z, vn)
    w = zeros(L, 1); % Weight initialization
    e2 = zeros(N, 1);
    exe2 = zeros(N, 1);

    for n = L:N
        xn = x(n:-1:n-L+1);         % Instant input
        en = d(n) - w' * xn;        % Estimation error
        exen = en - z(n);           % Excess error
        vn(:, n) = w - wopt;        % Weight error

        w = w + mu * xn * en;       % Filter update
        e2(n) = en^2;
        exe2(n) = exen^2;
    end
end

function [e2t, exe2t, vt] = theoreticalModel(N, L, mu, wopt, rho, varn)
    R = toeplitz(rho .^ [0:L-1]);  % Correlation matrix
    vt = zeros(L, N);
    vt(:, L) = -wopt;
    K = wopt * wopt';
    exe2t = zeros(N, 1);
    e2t = zeros(N, 1);

    for n = L+1:N
        vt(:, n) = (eye(L) - mu * R) * vt(:, n-1);  % First-order model
        K = (eye(L) - mu * R) * K * (eye(L) - mu * R) + mu^2 * varn * R; % Second-order model
        exe2t(n) = trace(K * R);
        e2t(n) = varn + exe2t(n);
    end
end

function plotResults(vnav, vt, wopt, e2av, e2t, exe2av, exe2t, L)
    % Weight evolution
    figure, hold on, grid on
    for i = 1:L
        plot(vnav(i, L:end) + wopt(i), 'linewidth', 2);
        plot(vt(i, L:end) + wopt(i), 'r', 'linewidth', 2);
    end
    xlabel('iteration', 'fontsize', 14, 'interpreter', 'Latex')
    ylabel('Weights', 'fontsize', 14, 'interpreter', 'Latex')

    % Mean square error
    figure, semilogy(e2av(L:end))
    grid on, hold on, semilogy(e2t, 'r', 'linewidth', 2)
    xlabel('iteration', 'fontsize', 14, 'interpreter', 'Latex')
    ylabel('MSE', 'fontsize', 14, 'interpreter', 'Latex')

    % Excess mean square error
    figure, semilogy(exe2av(L:end))
    grid on, hold on, semilogy(exe2t, 'r', 'linewidth', 2)
    xlabel('iteration', 'fontsize', 14, 'interpreter', 'Latex')
    ylabel('EMSE', 'fontsize', 14, 'interpreter', 'Latex')
end
