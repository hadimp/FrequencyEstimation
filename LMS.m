%% Implementation of the LMS algorithm

function [theta_n, thetas] = LMS(mu, M, N, r, f1, fs, theta_n, noise, SNR)
    thetas = zeros(N,1);
    for s = 1:N
        yn = CalcY(M, N, r, f1, fs, theta_n, noise, SNR);
        Beta = CalcBeta(theta_n, M, N, r, yn);
        Beta = Beta(s);
        y = yn(M, s);
        theta_n = theta_n - (2 * mu * y * Beta);
        thetas(s,1) = theta_n;
    end
end
