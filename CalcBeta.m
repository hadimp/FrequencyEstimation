%% Caculation of Gradient terms for LMS algorithm
function Beta = CalcBeta(theta_n, M, N, r, y)

    beta = zeros(M, N);
    for m = 2:M
        for n = 1:N
            if n==1
                beta(m, 1) = beta(m-1, 1);
            elseif n==2
                beta(m, 2) = beta(m-1, 2) - 2*cos((m-1)*theta_n)*beta(m-1, 1) + 2*(m-1)*sin((m-1)*theta_n)*y(m-1, 1) + 2*r*cos((m-1)*theta_n)*beta(m, 1) - 2*r*(m-1)*sin((m-1)*theta_n)*y(m, 1);
            else    
                beta(m, n) = beta(m-1, n) - 2*cos((m-1)*theta_n)*beta(m-1, n-1) + 2*(m-1)*sin((m-1)*theta_n)*y(m-1, n-1) + beta(m-1, n-2) + 2*r*cos((m-1)*theta_n)*beta(m, n-1) - (r^2)*beta(m, n-2) - 2*r*(m-1)*sin((m-1)*theta_n)*y(m, n-1);
            end
        end
    end
    Beta = beta(M,:); 
end