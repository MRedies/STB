Building on Thonhausers review (eq 21)


$$
\begin{align}
M_{orb} &= \frac{e}{2\hbar c}Im \sum_n \int \frac{d^3k}{(2\pi)^3} f_{nk} \left< \partial_k u_{nk}| \times (H_k + E_{nk} -2 \mu)| \partial_k u_{nk} \right> \\
 &= \frac{e}{2\hbar c}Im \sum_n \int \frac{d^3k}{(2\pi)^3} f_{nk} \left< \partial_k u_{nk}| \times (H_k - E_{nk})| \partial_k u_{nk} \right> \\
 &+ \frac{e}{2\hbar c}Im \sum_n \int \frac{d^3k}{(2\pi)^3} f_{nk} \left< \partial_k u_{nk}| \times (2E_{nk} - 2\mu)| \partial_k u_{nk} \right>
\end{align}
$$
then we can use: 
$$
\left|\partial u_{nk}\right> = \sum_{m \neq n} \frac{\left<u_{mk}|v|u_{nk} \right> }{E_{nk} - E_{mk}}\left|u_{mk}\right> =: \sum_{m \neq n} \frac{V_{mn}}{E_{nk} - E_{mk}}\left|u_{mk}\right>
$$
applying this to the first term:


$$
\begin{align}
	&\left< \partial_k u_{nk}| \times (H_k - E_{nk})| \partial_k u_{nk} \right>^i  \\
	=&\varepsilon_{ijk} \sum_{m \neq n} \frac{V_{nm}^j}{E_{nk} - E_{mk}}  \left< u_{mk}\right|H_k - E_{nk}  \sum_{m' \neq n} \frac{V_{mn}^k}{E_{nk} - E_{m'k}}\left|u_{m'k} \right> \\
	=&\varepsilon_{ijk}   \sum_{m \neq n} \frac{V_{nm}^j}{E_{nk} - E_{mk}} \sum_{m' \neq n} \frac{V_{m'n}^k}{E_{nk} - E_{m'k}}  \left< u_{mk}\right|H_k - E_{nk} \left|u_{m'k} \right> \\
	=&\varepsilon_{ijk}   \sum_{m \neq n} \frac{V_{nm}^j}{E_{nk} - E_{mk}} \sum_{m' \neq n} \frac{V_{m'n}^k}{E_{nk} - E_{m'k}}  \underbrace{\left< u_{mk}\right|E_{m'k} - E_{nk} \left|u_{m'k} \right>}_{-\delta_{ m m'} E_{nk} -E_{m'k}} \\
	=&\varepsilon_{ijk}   \sum_{m \neq n} \frac{V_{nm}^j}{E_{nk} - E_{mk}} \frac{V_{mn}^k}{E_{nk} - E_{mk}}  \left< u_{mk}\right|E_{mk} - E_{nk} \left|u_{mk} \right> \\
	=&-\varepsilon_{ijk}   \sum_{m \neq n} \frac{V_{nm}^j V_{mn}^k}{E_{nk} - E_{mk}}   \\
\end{align}
$$
while the other term leads to:
$$
\begin{align}
	 &\left< \partial_k u_{nk}| \times (2E_{nk} - 2\mu)| \partial_k u_{nk} \right>\\
	 =& \varepsilon_{ijk} \sum_{m \neq n} \frac{V_{nm}^j}{E_{nk} - E_{mk}}  \left<u_{mk} \right| \left( 2E_{nk} - 2\mu \right) \sum_{m' \neq n} \frac{V_{m'n}^k}{E_{nk} - E_{m'k}} \left|u_{m'k}\right> \\
	 =& \varepsilon_{ijk}  \sum_{m\neq n} \frac{V_{nm}^j  V_{mn}^k}{\left(E_{nk} - E_{mk} \right)^2} \left( 2E_{nk} - 2\mu \right)
\end{align}
$$
Leading to an overall:
$$
\begin{align}
M_{orb}^i &= \frac{e}{2\hbar c}Im \sum_n \int \frac{d^3k}{(2\pi)^3} f_{nk} \left< \partial_k u_{nk}| \times (H_k + E_{nk} -2 \mu)| \partial_k u_{nk} \right> \\
 &= \frac{e}{2\hbar c}Im \sum_n \int \frac{d^3k}{(2\pi)^3} f_{nk} \left< \partial_k u_{nk}| \times (H_k - E_{nk})| \partial_k u_{nk} \right> \\
 &+ \frac{e}{2\hbar c}Im \sum_n \int \frac{d^3k}{(2\pi)^3} f_{nk} \left< \partial_k u_{nk}| \times (2E_{nk} - 2\mu)| \partial_k u_{nk} \right> \\
 &= \varepsilon_{ijk}  \frac{e}{2\hbar c}Im \sum_n\sum_{m\neq n} \int \frac{d^3k}{(2\pi)^3} f_{nk} \left(\frac{-V_{nm}^j V_{mn}^k}{E_{nk} - E_{mk}} + 2 (E_{nk} - \mu ) \frac{V_{nm}^j V_{mn}^k}{\left(E_{nk} - E_{mk}\right)^2}\right) \\
  &= \varepsilon_{ijk}  \frac{e}{2\hbar c}Im \sum_n\sum_{m\neq n} \int \frac{d^3k}{(2\pi)^3} f_{nk} \left(\frac{V_{nm}^j V_{mn}^k}{\Delta E} - 2 (\mu -E_{nk}  ) \frac{V_{nm}^j V_{mn}^k}{\left(\Delta E\right)^2}\right) 
\end{align}
$$
where we define (notice switch in order):
$$
\Delta E := E_{mk} - E_{nk}
$$
