% linear interpolation
function[Rho_f]=finer(Rho,NN)
Rho_f=100*ones(2^NN*(length(Rho)-1)+1,1);
Rho_f(1)=Rho(1);
for jj=1:2^NN
%     Rho_f(jj+1:2^NN:end)=((2^NN-jj)*Rho(1:end-1)+jj*Rho(2:end))/2^NN;
    Rho_f(jj+1:2^NN:end)=Rho(1:end-1)+jj/2^NN*(Rho(2:end)-Rho(1:end-1));
end

