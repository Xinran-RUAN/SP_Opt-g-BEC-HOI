% interpolation via SP
function[Rho_f] = finer_SP(Rho,NN)
for jj = 1 : NN 
    N = length(Rho);
    Rho_shift = 2/(N+1)*dst_iii(dst(Rho),N+1);
    Rho_mat = [Rho_shift'; [Rho;0]'];
    Rho = Rho_mat(1:end-1)';
end
Rho_f=Rho;