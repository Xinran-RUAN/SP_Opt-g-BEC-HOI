# SP_Opt-g-BEC-HOI
This is the joint work with Bo Lin

We aim to design a spectrally accurate method for computing ground states of BEC with higher order interactions. 

The idea is to apply density-based gradient flow / optimization technique, combining with a postive preserving spectrally accurate projection technique.

## Optimization of Energy Functionals Based on Regularized Density Formulations

For 1D case, we consider the energy functional 
$$
	E(\rho) = \int_a^b \left[ \frac{|\nabla \rho|^2}{\rho+\varepsilon} + V(x) \rho + \frac{\beta}{2}\rho^2 + \frac{\delta}{2}|\nabla\rho|^2\right]dx, x\in[a,b].
$$
Denote $x_j=a+jh$, $j=0,1,...,N$, where $h=\frac{b-a}{N}$ and $\rho_j$ to be the approximation of $\rho(x_j)$, with the periodic boundary condition $\rho_0 = \rho_N$.

The energy functional can be discretized as 
$$
	E(\rho) = h\sum_{j=0}^{N-1} \left[ \frac{|\delta_x^s \rho_j|^2}{\rho_j+\varepsilon} + V(x) \rho_j + \frac{\beta}{2}\rho_j^2 + \frac{\delta}{2}|\delta_x^s \rho_j|^2\right], x\in[a,b].
$$