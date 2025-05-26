function[X_OddExt] = OddExtension(X)
[~,n] = size(X);
if n == 1
    X_OddExt = [-X(end:-1:1); 0; X; 0];
else
    X_OddExt = [-X(end:-1:1), 0, X, 0];
end
