function[FX, dFX] = F_X_Func(X, vep)
% FX = (X>vep) .* sqrt(X) ...
%    + (X<=vep) .* (X / sqrt(vep)) .* (1 + 1/2 * (1 - X ./ vep));
% dFX = (X>vep) ./ sqrt(X) / 2 ...
%     + (X<=vep) .* (3/2 - X / vep) / sqrt(vep);

FX = sqrt(X + vep);
dFX = 1 ./ sqrt(X + vep) / 2;
