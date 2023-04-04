%% Calculation of the output of subfilters
function y = CalcY(M, N, r, f1, fs, theta, noise, SNR)    
    
    N = N+1;
    M = M+1;
    x = zeros(1, N);
    % input to the filter bank
    for n = 1:N
        x(n) = x(n) + sin(2*pi*f1/fs*n) + 0.5 * cos(2*pi*2*f1/fs*n) - 0.25 * cos(2*pi*3*f1/fs*n);
    end

    
    if noise
        x = awgn(x, SNR);
    end
    % initializing y as a matrix of output values for M notch filters and N
    % samples of the input. Assumption is that if the index of matrix is
    % below the range of acceptable values(<1) we set y(m,u)=0 (Is this
    % assumption correct?)

    y = zeros(M, N);
    for m = 1:M
        if m == 1
            y(1,:) = x;
        else
            b = [1 -2*cos((m-1)*theta) 1];
            a = [1 -2*r*cos((m-1)*theta) r^2];
            x = y(m-1,:);
            y(m,:) = filter(b, a, x);
        end
    end
end