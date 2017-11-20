$$
\begin{align}
	\mathbf{m} &= \frac{1}{N_k \Omega} \sum_k \sum_n^{occ} \sum_\mu \mathbf{m}^\mu_{\mathbf{k}n}\\
	                      &= \frac{1}{N_k \Omega} \sum_k \sum_n \sum_\mu f(E_n) \cdot \mathbf{m}^\mu_{\mathbf{k}n}
\end{align}
$$

using $\Psi_m = \sum_n \phi_n c_n$:
$$
\begin{align}
\mathbf{m}^\mu_{\mathbf{k}n}&= -\frac{e}{2m_e} \left<\Psi_{\mathbf{k}n}| \mathbf{L}^\mu|\Psi_{\mathbf{k}n} \right> \\
                                                         &= -\frac{e}{2m_e} \sum_{m, o} c_m^n c_o^{n*} \left<\phi_m^n| \mathbf{L}^\mu| \phi_o^n\right> \\
                                                         &= -\frac{e}{2m_e} \sum_{m} |c_m^n|^2 \left<\phi_m^n| \mathbf{L}^\mu| \phi_m^n\right> 
 \end{align}
$$

