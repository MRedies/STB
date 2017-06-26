In this document I use:
$$ m_e = \hbar = e = c =1$$
Using the general tight-binding Hamiltonian:
$$ H_{i\mu}^{j\nu}(k) = \sum_n e^{i \vec{k}(\mathbf{r^{i\leftrightarrow j}})} H_{0i\mu}^{j\nu}$$
Which in my case boils down to:

$$ H_{i\mu}^{j\nu}(k)= \sum_n \left(  \varepsilon_i \delta_{ij} - J(...) \delta_{ij}\right)e^{i \vec{k} \vec{0}} + t_{ij} e^{i\vec{k} \vec{r_{ij}}} $$

The derivative with respect to the l-th k-space component:
$$ \partial_{k_l} H_{i\mu}^{j\nu}(k) = i \: r_{ij}^l \: t_{ij} e^{i \vec{k} \: \vec{r_{ij}}} $$

We want to calculate
$$
\sigma_{xy} = \sum_n \frac{1}{2\pi} \int_{BZ} d \mathbf{k} \Omega_z^n(\mathbf{k}) f(E_{n\mathbf{k}}) = \frac{V_k}{N_k} \frac{1}{2\pi} \sum_{n,\mathbf{k_i}} \Omega_z^n(\mathbf{k_i})f(E_{n\mathbf{k_i}})
$$
with 
$$
\begin{aligned}
 \mathbf{\Omega}^n_i &= \frac{1}{2} \epsilon_{ijk} \Omega^n_{jk}\\
\Rightarrow \mathbf{\Omega}_z^n &= \frac{1}{2} \left( \Omega^n_{xy} -\Omega^n_{yx}\right)
\end{aligned}
$$
$$ \Omega_{ij}^n = -2 Im \sum_{m \neq n} \frac{\langle u_{n\mathbf{k}}|\partial_{k_i}H(\mathbf{k})|u_{m\mathbf{k}}\rangle \langle u_{m\mathbf{k}}|\partial_{k_j}H(\mathbf{k})|u_{n\mathbf{k}} \rangle}{\left(\varepsilon_{n\mathbf{k}} -\varepsilon_{m\mathbf{k}} \right)^2}$$.



For this we need the matrix derivative $\partial_{k_i} H(\mathbf{k})  $. With the Hamiltonian


$$
\begin{aligned}
H^{j\nu}_{l\mu}(\mathbf{k}) &= \sum_n e^{i \mathbf{k} \cdot \mathbf{r_{lj}} }  H^{nj\nu}_{0_{l\mu}}  \\
&= \sum_n(\varepsilon_i \delta_{lj} - J_{lj}(...))e^{i \mathbf{k} \cdot \mathbf{0}} + t_{lj} e^{i \mathbf{k} \mathbf{r^{l\leftrightarrow j}}}
\end{aligned}$$
We get
$$ \partial_{k_i} H^{j\nu}_{l\mu}(\mathbf{k}) = i \: r^{l \leftrightarrow j}_i  t_{lj} e^{i \mathbf{k} \mathbf{r^{l\leftrightarrow j}}}.$$
With this the elements

\begin{aligned}
\langle u_{n\mathbf{k}}|\partial{k_i}H(\mathbf{k})|u_{m\mathbf{k}}\rangle &= i\langle u_{n\mathbf{k}}| r^{l \leftrightarrow j}  t_{lj} e^{i \mathbf{k} \mathbf{r^{l\leftrightarrow j}}}|u_{m\mathbf{k}}\rangle \\

&= i\langle \psi_{n\mathbf{k}}| e^{i\mathbf{k}\cdot \mathbf{r}}r^{l \leftrightarrow j}  t_{lj} e^{i \mathbf{k} \mathbf{r^{l\leftrightarrow j}}} e^{-i\mathbf{k}\cdot \mathbf{r}}|\psi_{m\mathbf{k}}\rangle \\

&= i\langle \psi_{n\mathbf{k}}| r^{l \leftrightarrow j}  t_{lj} e^{i \mathbf{k} \mathbf{r^{l\leftrightarrow j}}} |\psi_{m\mathbf{k}}\rangle 

\end{aligned}$$
$$

In the implementation define matrix 
$$
\partial H_i := i \mathbf{r}^{l\leftrightarrow j}_i t_{lj} e^{i \mathbf{k} \mathbf{r}^{l\leftrightarrow j}}
$$
