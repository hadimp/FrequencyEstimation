%% Calculate MSE and MSE1 based on outputs of subfilters
function [MSE, MSE1] = mse(M, N, y)
    N = N+1;
    M = M+1;

    MSE = 0;
    MSE1 = 0;

    for i = 1:N
        MSE = MSE + y(M,i)^2;
    end
    MSE = MSE / N;
    
    for i = 1:N
        MSE1 = MSE1 + y(2,i)^2;
    end
    MSE1 = MSE1 / N;
end
