function f_hat = dst1_scalar(f)
    % Compute DST-I coefficients using scalar operations
    % f: column vector of length N-1, representing f(x_j) on interior grid
    % Returns f_hat: column vector of DST-I coefficients (N-1 x 1)

    f = f(:);                % Ensure column
    N = length(f) + 1;       % N-1 interior points

    f_hat = zeros(N-1, 1);   % Preallocate

    for l = 1:N-1
        sum_val = 0;
        for j = 1:N-1
            sum_val = sum_val + f(j) * sin(pi * l * j / N);
        end
        f_hat(l) = (2 / N) * sum_val;
    end
end