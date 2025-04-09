function[QL] = myQ_L(X_temp, y, L, data)
[E, G] = myDiscreteEnergy(y, data);
QL = E + sum((X_temp - y) .* G) + 0.5 * L * sum((X_temp - y).^2);
