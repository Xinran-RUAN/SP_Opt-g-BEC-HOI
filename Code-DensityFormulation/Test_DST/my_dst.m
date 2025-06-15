% compute via odd extension
function[hat_X] = my_dst(X)
X_OddExt = OddExtension(X);
N = length(X);
hat_X_OddExt = fft(X_OddExt) / (2 * N + 2);
hat_X = -2 * imag(hat_X_OddExt(2:N+1));
