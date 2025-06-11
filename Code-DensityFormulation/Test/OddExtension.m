function[X_OddExt] = OddExtension(X)
[~,n] = size(X);
if n == 1
    X_OddExt = [0; X; 0; -X(end:-1:1)];
else
    X_OddExt = [0, X, 0, -X(end:-1:1)];
end
