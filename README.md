# SP_Opt-g-BEC-HOI
This is the joint work with Bo Lin

We aim to design a spectrally accurate method for computing ground states of BEC with higher order interactions. 

The idea is to apply density-based gradient flow / optimization technique, combining with a postive preserving spectrally accurate projection technique.

## Code-DensityFormulation

Optimization of Energy Functionals Based on Regularized Density Formulations

### Discrete energy functional
For 1D case, we consider the energy functional 
$$
	E(\rho) = \int_a^b \left[ \frac{|\nabla \rho|^2}{8(\rho+\varepsilon)} + V(x) \rho + \frac{\beta}{2}\rho^2 + \frac{\delta}{2}|\nabla\rho|^2\right]dx, x\in[a,b].
$$
Denote $$x_j=a+jh, j=0,1,...,N,$$ where $h=\frac{b-a}{N}$. The energy functional can be discretized as 
$$
	E_h(\rho) = h\sum_{j=0}^{N-1} \left[ \frac{|\delta_x^s \rho_j|^2}{8(\rho_j+\varepsilon)} + V(x_j) \rho_j + \frac{\beta}{2}\rho_j^2 + \frac{\delta}{2}|\delta_x^s \rho_j|^2\right],
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

**Proposition**: $E_h(\rho)$ is convex in $S_N^+$.

**Proof**: It is sufficient to prove the convexity of 
$$
	E_{kin}(\rho) =  \frac{h}{8}\sum_{j=0}^{N-1} \frac{|\delta_x^s \rho_j|^2}{\rho_j+\varepsilon}.
$$
A direct computation shows that 
$$
\begin{aligned}
	\frac{E_{kin}(a) + E_{kin}(b)}{2} & =  \frac{h}{16}\sum_{j=0}^{N-1} \left[\frac{|\delta_x^s a_j|^2}{a_j+\varepsilon} + \frac{|\delta_x^s b_j|^2}{b_j+\varepsilon}\right]\\
	& = \frac{h}{16}\sum_{j=0}^{N-1} \frac{|\delta_x^s a_j|^2\left(1 + \frac{b_j+\varepsilon}{a_j+\varepsilon}\right) + |\delta_x^s b_j|^2\left(1 + \frac{a_j+\varepsilon}{b_j+\varepsilon}\right)}{a_j + b_j +2\varepsilon}\\
	& \ge \frac{h}{16}\sum_{j=0}^{N-1} \frac{|\delta_x^s a_j|^2 + 2(\delta_x^s a_j)(\delta_x^s b_j) + |\delta_x^s b_j|^2}{a_j + b_j +2\varepsilon}\\
	& = \frac{h}{16}\sum_{j=0}^{N-1} \frac{(\delta_x^s a_j + \delta_x^s b_j)^2}{a_j + b_j +2\varepsilon}\\
	& = \frac{h}{8}\sum_{j=0}^{N-1} \frac{\left(\delta_x^s \frac{a_j+b_j}{2} \right)^2}{\frac{a_j + b_j}{2} + \varepsilon}\\
	& = E_{kin}\left(\frac{a+b}{2}\right)
\end{aligned}
$$

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

**Corollary 1**: $\hat{\rho}_{-\frac{N}{2}} = \hat{\rho}_{\frac{N}{2}} \in \mathbb{R}$.

**Proposition 2**: $Re(\delta_x^s \rho_j) = \sum_{l=-\frac{N}{2}+1}^{\frac{N}{2}-1} i\mu_l\hat{\rho}_l e^{i \frac{2\pi l j}{N}}$, $Im(\delta_x^s \rho_j) = i\mu_{-\frac{N}{2}}\hat{\rho}_{-\frac{N}{2}} e^{-i \pi j} = i\mu_{-\frac{N}{2}}\hat{\rho}_{-\frac{N}{2}}(-1)^j$

**proof**: For $k = -\frac{N}{2}$, $i\mu_{-\frac{N}{2}}\hat{\rho}_{-\frac{N}{2}} e^{-i \pi j} = i\mu_{-\frac{N}{2}}\hat{\rho}_{-\frac{N}{2}}(-1)^j$ is a pure imaginary number.
For $k\neq -\frac{N}{2}$, 
$$
	i\mu_k \hat{\rho}_k e^{i \frac{2\pi k j}{N}} + i\mu_{-k} \hat{\rho}_{-k} e^{i \frac{-2\pi k j}{N}} = i\mu_k (\hat{\rho}_k e^{i \frac{2\pi k j}{N}} - \overline{\hat{\rho}_{k} e^{i \frac{2\pi k j}{N}}}) \in \mathbb{R},
$$
where the periodicity and Proposition 1 are applied. Therefore, $ \sum_{l=-\frac{N}{2}+1}^{\frac{N}{2}-1} i\mu_l\hat{\rho}_l e^{i \frac{2\pi l j}{N}} \in \mathbb{R}$. 

**As a remark, in practice we only need to take the real part for computing $\delta_x^s \rho_j$.**

#### Relations via Matrices
The above relations can be described by using matrices. Define the diagonal matrix $\Lambda= \text{diag}\{i\mu_0, i\mu_1, \cdots, i\mu_{\frac{N}{2}-1}, i\mu_{-\frac{N}{2}}, i\mu_{-\frac{N}{2}+1}, \cdots, i\mu_{-1}\}$ and the symmetric matrix $F=(F_{j,k})\in\mathbb{R}^{N\times N}$, where 
$$
	F_{j,k} = e^{i\frac{2\pi jk}{N}},
$$
and define the vectors $\rho = (\rho_0, \rho_1, \cdots, \rho_{N-1})^T$ and $\hat{\rho} = (\hat{\rho}_0, \hat{\rho}_1, \cdots, \hat{\rho}_{N-1})^T$, then

| Relation | Matlab Code |
| ------- | ---------  |
|$\rho = F \hat{\rho}$| `rho = N * ifft(Frho)`|
|$\hat{\rho} = F^{-1} \rho = \frac{1}{N} \overline{F} \rho$ | `Frho = 1 / N * fft(rho)`|
|$\delta_x^s \rho = D\rho$ | `drho = fourier_diff(rho, Lambda);`|
| $D^T\rho$ | `drho_T = fourier_diff_T(rho, Lambda);`|

where $D = \mathcal{R}(\frac{1}{N} F \Lambda\overline{F}), D^T = \mathcal{R}(\frac{1}{N} \overline{F} \Lambda F)$ and 

```Matlab
function [drho] = fourier_diff(rho, Lambda)
    Frho = fft(rho);
    Fdrho = Lambda .* Frho;
    drho = real(ifft(Fdrho));
end
```

```Matlab
function [drho_T] = fourier_diff_T(rho, Lambda)
    iFrho = ifft(rho);
    iFdrho = Lambda .* iFrho;
    drho_T = real(fft(iFdrho));
end
```

### Gradient of the discrete energy functional 
Vectorization of the energy functional: 
$$
	E_h(\rho) = h\sum_{j=0}^{N-1} 
	\left[
	\frac{(D\rho)_j^2}{8(\rho_j+\varepsilon)} + V_j\rho_j +\frac{\beta}{2}\rho_j^2 + \frac{\delta}{2} (D\rho)_j^2,
	\right]
$$
where $D = \mathcal{R}(\frac{1}{N} F \Lambda\overline{F})$. A detailed computation shows the partial derivative of each $\rho_i$ as follows
$$
	\frac{\partial}{\partial \rho_i} \frac{(D\rho)_j^2}{8(\rho_j+\varepsilon)} 
	=
	\frac{(D\rho)_j D_{ji}}{4(\rho_j+\varepsilon)}
	-\frac{(D\rho)_j^2}{8(\rho_j+\varepsilon)^2}\delta_{ij}
$$
$$
	\frac{\partial}{\partial \rho_i}(V_j\rho_j) = V_j \delta_{ij}
$$
$$
	\frac{\partial}{\partial \rho_i}\left(\frac{\beta}{2}\rho_j^2\right) = \beta \rho_j \delta_{ij}
$$
$$
	\frac{\partial}{\partial \rho_i}\left(\frac{\delta}{2}(D\rho)_j^2\right) = \delta (D\rho)_j D_{ji}.
$$
To summarize 

| Energy | Gradient | Matlab Code |
| ------- | ------- | ---------  |
|$h\sum_{j=0}^{N-1} \frac{(D\rho)_j^2}{\rho_j+\varepsilon}$| $h \left[ D^T\frac{D\rho}{4(\rho+\varepsilon)} - \frac{(D\rho)^2}{8(\rho+\varepsilon)^2}\right]$ | `Drho = fourier_diff(rho, Lambda);` <br> `Grho = Drho ./ (rho + vep);` <br> `DtGrho = fourier_diff_T(Grho, Lambda);` <br> `dE_kin = h * (DtGrho / 4 - Grho.^2 / 8);` |
|$h\sum_{j=0}^{N-1} V_j\rho_j$ | $hV$ | `dE_pot = h * V;`|
|$h\sum_{j=0}^{N-1} \frac{\beta}{2}\rho_j^2$ | $h\beta \rho$ | `dE_beta = h * beta * rho` |
|$h\sum_{j=0}^{N-1} \frac{\delta}{2} (D\rho)_j^2$ | $h\delta D^T D\rho$ | `Drho = fourier_diff(rho, Lambda);` <br> `DtDrho = fourier_diff_T(Drho, Lambda);` <br> `dE_delta = h * delta * DtDrho;` |

**Proposition**: $\nabla E(\rho)$ is Lipschitz continuous in $\rho$, i.e. 
$$
	\|\nabla E(x) - \nabla E(y)\| \le L \|x-y\|.
$$
**Proof**: We only need to show the Lipschitz continuity of the kinetic energy term. 
$$
\begin{aligned}
	\| D^T\frac{D x}{4(x+\varepsilon)} - D^T\frac{D y}{4(y+\varepsilon)} \| & \le \frac{\|D\|}{4} \left\|  \frac{D x}{x+\varepsilon} - \frac{D y}{y+\varepsilon} \right\| \\
	& \le \frac{\|D\|}{4} \left(\left\|\frac{Dx - Dy}{x+\varepsilon}\right\| + \left\|Dy \frac{y-x}{(x+\varepsilon)(y+\varepsilon)}\right\|\right)\\
	& \le \frac{\|D\|}{4} \left(\frac{\|D\|}{\varepsilon}\|y-x\| + \frac{\|Dy\|_\infty}{\varepsilon^2} \|y-x\|\right)\\
	& \le \frac{\|D\|}{4} \left(\frac{\|D\|}{\varepsilon} + \frac{\|D\|_\infty \|y\|_\infty}{\varepsilon^2} \right)\|y-x\|\\
	& \le \frac{\|D\|}{4} \left(\frac{\|D\|}{\varepsilon} + \frac{\|D\|_\infty}{h\varepsilon^2} \right)\|y-x\|
\end{aligned}
$$
where $\|D\|_\infty = max_j \sum_{k=1}^N |D_{jk}|$.

## ISTA/FISTA
Key functions used in the algorithm:

*  $Q_L(x, y)$: Quadratic approximation of $f(x)$ defined as (linear approximation with quadratic penalty)
$$
	Q_L(x, y) := f(y) + <x - y, \nabla f(y)> + \frac{L}{2}\|x-y\|^2.
$$

*  $p_L(y)$: Unique minimizer of $Q_L(x,y)$ in the feasible set $S_N^+(1)$, where $S_N^+(1) = S_N^+ \cap \{X\in\mathbb{R}^N \, | \,  h \sum_{i=1}^N X(i)= 1\}$. 
Simple algebra shows that 
$$
	p_L(y) = argmin_{x\in S_N^+(1)} \left\|x - \left(y - \frac{1}{L} \nabla f(y)\right)\right\|.
$$

**Algorithm 1.1**: ISTA with constant stepsize
Step k: 
$$
	X_k = p_L(X_{k-1})
$$
Here $L$ can be taken as $L(f)$, a Lipschitz constant of $\nabla f$.

**Algorithm 1.2**: ISTA with backtracking
$$
	L_0 > 0, \, \eta > 1
$$
Step k: 
Find smallest $i_k$ such that with $\overline{L} = \eta^{i_k} L_{k-1}$, 
$$
	f(p_{\overline{L}}(X_{k-1})) \le Q_{\overline{L}}(p_{\overline{L}}(X_{k-1}), X_{k-1}).
$$
Set $L_k = \overline{L}$ and compute
$$
	X_k = p_{L_k}(X_{k-1}).
$$

**Algorithm 2.1**: FISTA with constant stepsize
$$
	Y_1 = X_0, \, t_1 = 1
$$
Step k: 
$$
\begin{aligned}
	&X_k = p_L(Y_k)\\
	&t_{k+1} = \frac{1+\sqrt{1+4t_k^2}}{2}\\
	&Y_{k+1} = X_k + \left(\frac{t_k-1}{t_{k+1}}\right)(X_k - X_{k-1})
\end{aligned}
$$

**Algorithm 2.2**: FISTA with backtracking
$$
	L_0 > 0, \, \eta > 1, \, y_1 = x_0, \, t_1 = 1
$$
Step k: 
Find smallest $i_k$ such that with $\overline{L} = \eta^{i_k} L_{k-1}$, 
$$
	f(p_{\overline{L}}(Y_{k})) \le Q_{\overline{L}}(p_{\overline{L}}(Y_{k}), Y_{k}).
$$
Set $L_k = \overline{L}$ and compute
$$
\begin{aligned}
	&x_k = p_L(y_k)\\
	&t_{k+1} = \frac{1+\sqrt{1+4t_k^2}}{2}\\
	&y_{k+1} = x_k + \left(\frac{t_k-1}{t_{k+1}}\right)(x_k - x_{k-1})
\end{aligned}
$$

### Convergence analysis
**Proposition 1**:If $f$ is convex and $\|\nabla f(x) - \nabla f(y) \| \le L_0 \|x - y\|$, then 
$$
\begin{aligned}
	& f(x) \ge f(y) + <x - y, \nabla f(y)> \\
	& f(x) \le f(y) + <x - y, \nabla f(y)> + \frac{L_0}{2}\|x-y\|^2 = Q(x, y).
\end{aligned}
$$
**Proof**: We only show the second inequality here.
$$
\begin{aligned}
	f(x) - f(y) - <x - y, \nabla f(y)> &= \int_0^1 <\nabla f(y+t(x-y)) - \nabla f(y), x-y> dt \\
	& \le L_0 \int_0^1 t \|x-y\|^2 dt \\
	& = \frac{L_0}{2} \|x-y\|^2.
\end{aligned}
$$

**Proposition 2**: For any $x \in \mathbb{R}^N,y \in S_N^+(1)$, if 
$$
	f(p_L(y)) \le Q(p_L(y), y) = f(y) + <p_L(y) - y, \nabla f(y)> + \frac{L}{2}\|p_L(y)-y\|^2,
$$
then 
$$
	F(x) - F(p_L(y)) \ge \frac{L}{2} \|p_L(y)-y\|^2 + L<y-x, p_L(y)-y>, 
$$
where $F(x) := f(x) + \mathbb{I}_{S_N^+(1)}(x)$.
**Proof**: If $x\notin S_N^+(1)$, then $F(x)=+\infty$ and the inequality is trivial. 
If $x\in S_N^+(1)$, then
$$
\begin{aligned}
	F(x) - F(p_L(y)) & = f(x) - f(p_L(y)) \\
	& \ge f(x) - f(y) - <p_L(y) - y, \nabla f(y)> - \frac{L}{2}\|p_L(y)-y\|^2 \\
	& \ge <x - p_L(y), \nabla f(y)> - \frac{L}{2}\|p_L(y)-y\|^2 .
\end{aligned}
$$

If $y - \frac{1}{L} \nabla f(y) \in S_N^+(1)$, then 
$$
	p_L(y) = y - \frac{1}{L} \nabla f(y)
$$
and therefore
$$
	\nabla f(y) = L(y - p_L(y)).
$$
By substitution, we get
$$
\begin{aligned}
	f(x) - f(p_L(y)) &\ge L<x - p_L(y), y - p_L(y)> - \frac{L}{2}\|p_L(y)-y\|^2 \\
	& = \frac{L}{2} \|p_L(y)-y\|^2 + L<y-x, p_L(y)-y>.
\end{aligned}
$$

If $y - \frac{1}{L} \nabla f(y) \notin S_N^+(1)$, by definition of $p_L(y)$, we have 
$$
	\left\|x -\left( y - \frac{1}{L}\nabla f(y)\right)\right\| \ge  \left\|p_L(y) -\left( y - \frac{1}{L}\nabla f(y)\right)\right\|,
$$
which implies 
$$
	<x-p_L(y), \nabla f(y)> \, \ge \frac{L}{2}\|p_L(y)-y\|^2 - \frac{L}{2}\|x-y\|^2 
$$
_**Fails!**_

**Lemma**: If $x^*$ is a local minimizer of $f(x)$ in $\Omega$, then we must have
$$
	-\nabla f \in N_\Omega(x^*),
$$
where $N_\Omega(x^*)$ is the normal cone defined as 
$$
	N_\Omega(x) = \{ v \, | \, v^Tw \le 0, \forall w \in T_\Omega(x) \}.
$$
Here $T_\Omega(x)$ is the tangent cone, the set of all tangents at $x$. We say $d$ is a tangent at $x$ if there exist $\{z_k\}\to x$ and $\{t_k\}\to0$ such that 
$\lim_{k\to\infty} \frac{z_k-x}{t_k} = d$. 

With the above lemma, noticing the definition of $p_L(y)$, we must have 
$$
	p_L(y) - \left(y-\frac{1}{L} \nabla f(y)\right) \in N_\Omega(p_L(y)).
$$
Here $\Omega = S_N^+(1)$ is the feasible set. Define 
$$
	\gamma = -\nabla f(y) + L(y - p_L(y)),
$$
then 
$$
	<\gamma, w> \le 0, \, \forall w\in T_\Omega(p_L(y)).
$$
It is easy to check that for all $x\in\Omega$, $x-p_L(y) \in T_\Omega(p_L(y))$, which implies 
$$
	<x-p_L(y), \gamma> \le 0.
$$
Then 
$$
\begin{aligned}
	f(x) - f(p_L(y)) \ge f(x) - f(p_L(y)) + <x-p_L(y), \gamma> \\
\end{aligned}
$$


## Numerical Test
$V(x)=\frac{x^2}{2}$, $\beta = 10$, $\delta = 10$, $h_0=1/2$

|   | $h=h_0$ |$h = h_0/2$|$h=h_0/4$|$h=h_0/8$| $h = 10^{-3}$|
|-------|-------|-------|-------|-------|-------|
| $\varepsilon = 10^{-4}$ | 2.228056825443657 | 2.228060264817093 |2.228060336933893 | 2.228060335770537
| $\varepsilon = 10^{-3}$ | 2.204129497462118 | 2.204062679240165 | 2.204062800910641 | 2.204062789539389 | 2.204062787481187|

We take the solution with $h=10^{-3}$ as the exact solution.

|  $\varepsilon = 10^{-3}$   | $h=h_0$ |$h = h_0/2$|$h=h_0/4$|$h=h_0/8$| 
|-------|-------|-------|-------|-------|
| $Abs(E(\rho)-E(\rho_{ex}))$ | 6.67E-5 | 1.08E-7 | 1.34E-8 | 2.06E-9 |
| L2-norm | 3.87E-3 | 1.02E-3 | 2.52E-4 | 5.91E-5 |

## Future work

1. Bad behavior for small $\varepsilon$: new regularization / JKO (possible? uniformly effective)