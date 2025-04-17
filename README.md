# SP_Opt-g-BEC-HOI
This is the joint work with Bo Lin

We aim to design a spectrally accurate method for computing ground states of BEC with higher order interactions. 

The idea is to apply density-based gradient flow / optimization technique, combining with a postive preserving spectrally accurate projection technique.

## Code-DensityFormulation

Optimization of Energy Functionals Based on Regularized Density Formulations

### Discrete energy functional
For 1D case, we consider the energy functional 
$$
	E(\rho) = \int_a^b \left[ \frac{|\nabla \rho|^2}{\rho+\varepsilon} + V(x) \rho + \frac{\beta}{2}\rho^2 + \frac{\delta}{2}|\nabla\rho|^2\right]dx, x\in[a,b].
$$
Denote $$x_j=a+jh, j=0,1,...,N,$$ where $h=\frac{b-a}{N}$. The energy functional can be discretized as 
$$
	E_h(\rho) = h\sum_{j=0}^{N-1} \left[ \frac{|\delta_x^s \rho_j|^2}{\rho_j+\varepsilon} + V(x_j) \rho_j + \frac{\beta}{2}\rho_j^2 + \frac{\delta}{2}|\delta_x^s \rho_j|^2\right],
$$
where $\rho_j=\rho(x_j)$, satisfying the periodic boundary condition, and $\delta_x^s \rho_j$ is the pseudospectral approximation of $\nabla \rho(x_j)$. 
We aim to find $\rho_g\in S_N^+$ such that 
$$
	\rho_g = {argmin}_{\rho \in S_N^+} E_h(\rho), 
$$
where the discrete function space $S_N^+$ is defined as  
$$
	S_N^+ := span\{\Phi_l(x) = e^{i\frac{2\pi l(x-a)}{b-a}}, l = 0,1,2,\cdots, N-1 \} \cap \{ f(x) : f(x_j) \ge 0, j = 0,1,2,\cdots, N-1 \}.
$$
Here we do not force $\rho_g(x)$ to be positive over $[a,b]$. Instead, we only require $\rho_g(x)$ to be positive on the grids.

### Details of the pseudospectral discretization
For simplicity of notations, we use $\mu_l =\frac{2\pi l}{b-a}$. Since 
$$
	\rho(x) \approx \sum_{l=-\frac{N}{2}}^{\frac{N}{2}-1} \hat{\rho}_l e^{i\frac{2\pi l(x-a)}{b-a}} = \sum_{l=-\frac{N}{2}}^{\frac{N}{2}-1} \hat{\rho}_l e^{i\mu_l(x-a)}，
$$
the first derivatives at grid points can be computed as 
$$
	\frac{d \rho}{dx}(x) \approx \delta_x^s \rho := \sum_{l=-\frac{N}{2}}^{\frac{N}{2}-1} i\mu_l\hat{\rho}_l e^{i\mu_l(x_j-a)}
	= \sum_{l=-\frac{N}{2}}^{\frac{N}{2}-1} i\mu_l\hat{\rho}_l e^{i \frac{2\pi l j}{N}}.
$$
At grid points, 

#### Relations via Matrices
* $\rho = F$
 

### Details of 



### Future work

1. Bad behavior for small $\varepsilon$: new regularization / JKO (possible? uniformly effective)