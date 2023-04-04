%%
% WRITTEN BY: HADI MOHAMMADPOUR
% Prepared for: DR. C.Moloney
% ENGI 9821 - DSP
% Winter 2021 - Memorial University of Newfoundland
% 
% This file is the main code for the design assignment which aims to
% implement the LMS algorithm for frequency estimation and tracking. This
% is based on the paper title "Novel Adaptive IIR Filter for Frequency
% Estimation and Tracking" by Li Tan and Jean Jiang.
%
%% INPUTS
% noise: Indicates whether consider the SNR factor and add noise to the
%        signal or not.
% SNR: Signal to noise ratio for the input
% M: Number of subfilters (one for fundamental frequency and M-1 for
%    harmonics)
% r: Radius of the poles of transfer function
% fs: Sampling frequency of the signal
% N: Number of samples to consider
% f1: Fundamental frequency of the input
% n_points: Number of points in range [0 pi/m] to consider for theta to
%          calculate MSE, MSE1
% mu: Step size for LMS algorithm
%
%%
clear
close all

%% Parameters' list
noise = 0;
SNR = 18;
M = 3;
r = 0.95;
fs = 8000;
N = 400;
f1 = 1000;
n_points = 1400;
mu = 0.0001;

%% Initializing MSE, MSE1, and average
MSE = zeros(1, n_points);
MSE1 = zeros(1, n_points);
average = 0;

%% MSE, and MSE1 calculation
%This part of the code calls different functions, and tries to 
%calculate and plot MSE and MSE1 based on the y values. 
theta = linspace(0, pi/M, n_points);
figure(1)
for k = 1:n_points
    f = theta * fs / (2*pi);
    y = CalcY(M, N, r, f1, fs, theta(1,k),noise, SNR);
    [MSE(1, k), MSE1(1,k)] = mse(M, N, y);
    
    plot(f, MSE);
    hold on
    plot(f, MSE1);
    hold on
    average = average + (MSE(1,k)/n_points);
    yline(average);
    hold off
    legend('MSE', 'MSE1', 'Average MSE');
    title('MSE, MSE1, and Average MSE for 0 < theta < pi/M');
    xlabel('Frequency (Hz)');
    ylabel('MSE');
end
saveas(1,'MSE.png')
%% Approximating initial Theta value
[Min, Index] = min(MSE);

% Finding the intersection points of MSE & average and MSE1 & average and
% making sure that it is in global minimum valley

Capture_range = theta((MSE1 - average < 0.0001) & (MSE - average < 0.0001) & (MSE - Min < 0.2));

init_theta = Capture_range(1);
%init_theta = 900 * 2*pi / fs;


%% Applying LMS for frequency estimation and tracking

theta_n = init_theta;
[theta_n, thetas] = LMS(mu, M, N, r, f1, fs, theta_n, noise, SNR);


%% Calculating and plotting filter's output

y_final = CalcY(M, N, r, f1, fs, theta_n, noise, SNR);
freqs = thetas * fs / (2*pi);

n = linspace(1, N+1, N+1);
figure(2)
for i = 1:M+1
    subplot(M+1,1, i)
    plot(n, y_final(i,:));
    xlabel('Frequency (Hz)');
    if i==1
        ylabel('x(n)');
    else
        ylabel('y(n)');
    end
    saveas(2,'signals.png')

end
i = linspace(1, N, N);
figure(3)
plot(i, freqs);
title('Frequency Tracking');
xlabel('Iterations (n)');
ylabel('Frequency (Hz)');
ylim([freqs(1)-100 freqs(N)+100]);
saveas(3,'freq_tracking.png')


%% Plotting magnitude of the frequency response

figure(4)
H_final = 1;
for m = 1:M+1
    b = [1 -2*cos(m*theta_n) 1];
    a = [1 -2*r*cos(m*theta_n) r^2];
    w = linspace(-pi, pi, n_points*3);
    [H, ~] = freqz(b, a, w);
    A = 1/abs(H(n_points));
    H = H * A;
    H_final = H_final .* H;
    subplot(M+1, 1, m)
    plot(w/pi, abs(H));
    xlabel('w / pi');
    ylabel('|H(f)|');
    ylim([0 1.2])
end

A = 1/abs(H_final(n_points));
H_final = H_final * A;
subplot(M+1, 1, M+1)
plot(w/pi, abs(H_final));
title('Overall System');
xlabel('w / pi');
ylabel('|H(f)|');
ylim([0 1.2])
saveas(4,'system_response.png')

figure(5)
f = w * fs / (2*pi);
plot(f, abs(H_final));
title('Overall System for 0 < f < fs/2');
xlabel('w / pi');
ylabel('|H(f)|');
ylim([0 1.2])
saveas(5,'Overall_system.png')


H_overall_db = mag2db(H_final);
figure(6)
plot(f, H_overall_db);
title('Overall System for 0 < f < fs/2');
xlabel('w / pi');
ylabel('|H(f)| (dB)');
saveas(6,'Overall_system_db.png')