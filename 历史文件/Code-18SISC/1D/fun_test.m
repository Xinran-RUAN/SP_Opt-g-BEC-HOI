function[E,grad_E]=fun_test(X)
E=sum(X.^2); grad_E=2*X;