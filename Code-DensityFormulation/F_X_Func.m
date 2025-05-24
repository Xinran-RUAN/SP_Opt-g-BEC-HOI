function[FX, dFX] = F_X_Func(X, vep)
FX = sqrt(X + vep);
dFX = 1 ./ sqrt(X + vep) / 2;

% =========================================================
% % C1拼接，不够光滑
% FX = (X>vep) .* sqrt(X) ...
%    + (X<=vep) .* (X / sqrt(vep)) .* (1 + 1/2 * (1 - X ./ vep));
% dFX = (X>vep) ./ sqrt(X) / 2 ...
%     + (X<=vep) .* (3/2 - X / vep) / sqrt(vep);
% =========================================================
% % 包含导数爆炸部分
% delta = vep / 10;
% 
% % 平滑拼接函数
% S = 0.5 * (1 + tanh((X - vep)/delta));
% dS = (1 - tanh((X - vep)/delta).^2) / (2*delta);
% 
% % 多项式 P(X) 近似 sqrt(X) 在 [0, vep]
% P = (X / sqrt(vep)) .* (1 + 0.5 * (1 - X / vep));
% dP = (1 / sqrt(vep)) .* (3/2 - X / vep);
% 
% % 构造 F 和导数
% FX  = S .* sqrt(X) + (1 - S) .* P;
% dFX = S .* (0.5 ./ sqrt(X)) + (1 - S) .* dP + (sqrt(X) - P) .* dS;

% figure(1)
% plot(FX,'b'); hold on;
% plot(FX1,'r--'); hold off
% 
% figure(2)
% plot(dFX,'b'); hold on;
% plot(dFX1,'r--'); hold off
% pause(0.1)
