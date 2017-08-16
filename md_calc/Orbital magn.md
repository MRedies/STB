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
	=&\varepsilon_{ijk} \sum_{m \neq n} \frac{V_{mn^j}}{E_{nk} - E_{mk}}  \left< u_{mk}\right|H_k - E_{nk}  \sum_{m' \neq n} \frac{V_{mn^k}}{E_{nk} - E_{m'k}}\left|u_{m'k} \right> \\
	=&\varepsilon_{ijk}   \sum_{m' \neq n} \frac{V_{mn^k}}{E_{nk} - E_{m'k}} \sum_{m \neq n} \frac{V_{mn^j}}{E_{nk} - E_{mk}}  \left< u_{mk}\right|H_k - E_{nk} \left|u_{m'k} \right> \\
	=&\varepsilon_{ijk}   \sum_{m' \neq n} \frac{V_{mn^k}}{E_{nk} - E_{m'k}} \sum_{m \neq n} \frac{V_{mn^j}}{E_{nk} - E_{mk}}  \left< u_{mk}\right|E_{m'k} - E_{nk} \left|u_{m'k} \right> \\
	=&\varepsilon_{ijk}   \sum_{m \neq n} \frac{V_{mn^k}}{E_{nk} - E_{mk}} \frac{V_{mn^j}}{E_{nk} - E_{mk}}  \left< u_{mk}\right|E_{mk} - E_{nk} \left|u_{m'k} \right> \\
	=&\varepsilon_{ijk}   \sum_{m \neq n} \frac{V_{mn^k} V_{mn^j}}{E_{nk} - E_{mk}}   \\
\end{align}
$$
