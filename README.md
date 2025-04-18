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
	\rho(x) \approx \sum_{l=-\frac{N}{2}}^{\frac{N}{2}-1} \hat{\rho}_l e^{i\mu_l(x-a)}，
$$
the first derivatives at grid points can be computed as 
$$
	\frac{d \rho}{dx}(x) \approx \delta_x^s \rho := \sum_{l=-\frac{N}{2}}^{\frac{N}{2}-1} i\mu_l\hat{\rho}_l e^{i\mu_l(x-a)}.
$$
At grid points, 
$$
	\rho_j = \sum_{l=-\frac{N}{2}}^{\frac{N}{2}-1} \hat{\rho}_l e^{i\frac{2\pi l j}{N}}, 
	\quad 
	\delta_x^s \rho_j = \sum_{l=-\frac{N}{2}}^{\frac{N}{2}-1} i\mu_l\hat{\rho}_l e^{i \frac{2\pi l j}{N}}.
$$
Conversely,
$$
	\hat{\rho}_l = \frac{1}{N} \sum_{j=0}^{N-1} \rho_j e^{-i \frac{2\pi jl}{N}}.
$$
It is easy to see that 
$$
	\hat{\rho}_{l+N} = \hat{\rho}_{l}.
$$

**Proposition 1**：If $\rho\in\mathbb{R}^N$, then $\hat{\rho}_{N-k} = \overline{\hat{\rho}}_k$. （共轭对称性）
**Proof**: 
$$
	\hat{\rho}_{N-k} = \frac{1}{N} \sum_{j=0}^{N-1} \rho_j e^{-i \frac{2\pi j(N-k)}{N}} = \frac{1}{N} \sum_{j=0}^{N-1} \rho_j e^{i \frac{2\pi jk}{N}} = \overline{\hat{\rho}}_k.
$$

**Corollary 1**: $\hat{\rho}_{-\frac{N}{2}} = \hat{\rho}_{\frac{N}{2}} \in \mathbb{R}$, $i\mu_{-\frac{N}{2}}\hat{\rho}_{-\frac{N}{2}} e^{-i \pi j} = i\mu_{-\frac{N}{2}}\hat{\rho}_{-\frac{N}{2}}(-1)^j$ is a pure imaginary number.

**Proposition 2**: $Re(\delta_x^s \rho_j) = \sum_{l=-\frac{N}{2}+1}^{\frac{N}{2}-1} i\mu_l\hat{\rho}_l e^{i \frac{2\pi l j}{N}}$, $Im(\delta_x^s \rho_j) = i\mu_{-\frac{N}{2}}\hat{\rho}_{-\frac{N}{2}} e^{-i \pi j} = i\mu_{-\frac{N}{2}}\hat{\rho}_{-\frac{N}{2}}(-1)^j$

**proof**: For $k\neq -\frac{N}{2}$, 
$$
	i\mu_k \hat{\rho}_k e^{i \frac{2\pi k j}{N}} + i\mu_{-k} \hat{\rho}_{-k} e^{i \frac{-2\pi k j}{N}} = i\mu_k (\hat{\rho}_k e^{i \frac{2\pi k j}{N}} - \overline{\hat{\rho}_{k} e^{i \frac{2\pi k j}{N}}}) \in \mathbb{R},
$$
where the periodicity and Proposition 1 are applied. Therefore, $ \sum_{l=-\frac{N}{2}+1}^{\frac{N}{2}-1} i\mu_l\hat{\rho}_l e^{i \frac{2\pi l j}{N}} \in \mathbb{R}$. By corollary, for $k= -\frac{N}{2}$, the pure imaginary. The conclusion follows immediately.

**As a remark, in practice we only need to take the real part for computing $\delta_x^s \rho_j$.**

#### Relations via Matrices
The above relations can be described by using matrices. Define the diagonal matrix $D= \text{diag}\{i\mu_0, i\mu_1, \cdots, i\mu_{\frac{N}{2}-1}, i\mu_{-\frac{N}{2}}, i\mu_{-\frac{N}{2}+1}, \cdots, i\mu_{-1}\}$ and the symmetric matrix $F=(F_{j,k})\in\mathbb{R}^{N\times N}$, where 
$$
	F_{j,k} = e^{i\frac{2\pi jk}{N}},
$$
and define the vectors $\rho = (\rho_0, \rho_1, \cdots, \rho_{N-1})^T$ and $\hat{\rho} = (\hat{\rho}_0, \hat{\rho}_1, \cdots, \hat{\rho}_{N-1})^T$, then

| Relation | Matlab Code |
| ------- | ---------  |
|$\rho = F \hat{\rho}$| `rho = N * ifft(Frho)`|
|$\hat{\rho} = F^{-1} \rho = \frac{1}{N} \overline{F} \rho$ | `Frho = 1 / N * fft(rho)`|
|$\delta_x^s \rho = \frac{1}{N} F D\overline{F}\rho$ | `Frho = fft(rho);`<br>`Fdrho = D * Frho;`<br>`drho = real(ifft(Frho));`|


 

### Gradient of the discrete energy functional 



### Future work

1. Bad behavior for small $\varepsilon$: new regularization / JKO (possible? uniformly effective)