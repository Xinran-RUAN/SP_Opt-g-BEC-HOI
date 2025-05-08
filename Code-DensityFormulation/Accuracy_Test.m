% check accuracy
N_test = 4;
load MGPE-FP1d-Bet-10-Del-10-Vep-0.001-dx-0.001.mat;
rho_ex = [Rho; Rho(1)]; x_ex = [data.x; data.xmax]; dx = data.dx;
RHO = zeros(length(rho_ex), N_test);

load MGPE-FP1d-Bet-10-Del-10-Vep-0.001-dx-0.5.mat;
rho1 = [Rho; Rho(1)]; x1 = [data.x; data.xmax];
RHO(:,1) = interp1(x1, rho1, x_ex, 'linear');
load MGPE-FP1d-Bet-10-Del-10-Vep-0.001-dx-0.25.mat;
rho2 = [Rho; Rho(1)]; x2 = [data.x; data.xmax];
RHO(:,2) = interp1(x2, rho2, x_ex, 'linear');
load MGPE-FP1d-Bet-10-Del-10-Vep-0.001-dx-0.125.mat;
rho3 = [Rho; Rho(1)]; x3 = [data.x; data.xmax];
RHO(:,3) = interp1(x3, rho3, x_ex, 'linear');
load MGPE-FP1d-Bet-10-Del-10-Vep-0.001-dx-0.0625.mat;
rho4 = [Rho; Rho(1)]; x4 = [data.x; data.xmax];
RHO(:,4) = interp1(x4, rho4, x_ex, 'linear');

err_list = zeros(1,4);
for kk =1:4
    err_list(kk) = sqrt(dx * sum((rho_ex - RHO(:,kk)).^2));
end
