Starting from Jans Hamiltonian:
$$
H = -t \sum_{<ij>\alpha} c^\dagger_{i\alpha}c_{j\alpha} + i t_{so} \sum_{<ij>\alpha \beta} \mathbf{\hat{e}}_z \cdot (\sigma \times d_{ij}) c^\dagger_{i\alpha}c_{j\beta} + \lambda \sum_{i\alpha\beta}(\mathbf{\hat{m}} \cdot \sigma) c^\dagger_{i\alpha} c_{i\beta} - \lambda_{nl} \sum_{<ij> \alpha \beta} (\mathbf{\hat{m}} \cdot \sigma) c^\dagger_{i\alpha} c_{j\beta}
$$
we neglect the last term, because it is no clear for a non-ferro magnetic case:
$$
H = -t \sum_{<ij>\alpha} c^\dagger_{i\alpha}c_{j\alpha} + i t_{so} \sum_{<ij>\alpha \beta} \mathbf{\hat{e}}_z \cdot (\sigma \times d_{ij}) c^\dagger_{i\alpha}c_{j\beta} + \lambda \sum_{i\alpha\beta}(\mathbf{\hat{m}} \cdot \sigma) c^\dagger_{i\alpha} c_{i\beta}
$$
we can evaluate:
$$
H := i t_{so} \sum_{<ij>\alpha \beta} \mathbf{\hat{e}}_z \cdot (\sigma \times d_{ij}) c^\dagger_{i\alpha}c_{j\beta}  = i t_{so} \sum_{<ij>\alpha \beta} \left( \sigma_x \cdot d^y_{ij} - \sigma_y \cdot d^x_{ij} \right) c^\dagger_{i\alpha}c_{j\beta}
$$
To see if the model is hermitian we want:
$$
\overline{h_{ij}^{\alpha \beta}} = h_{ji}^{\beta \alpha}
$$

$$
\begin{align}
\overline{h_{ij}^{\alpha \beta}} &= \overline{ i t_{so} \sum_{<ij>\alpha \beta} \left( \sigma_x^{\alpha \beta} \cdot d^y_{ij} - \sigma_y^{\alpha \beta} \cdot d^x_{ij} \right) c^\dagger_{i\alpha}c_{j\beta}} \\
&= - i t_{so} \sum_{<ij>\alpha \beta} \overline{\left( \sigma_x^{\alpha \beta} \cdot d^y_{ij} - \sigma_y^{\alpha \beta} \cdot d^x_{ij} \right)} c_{i\alpha}c^\dagger_{j\beta} \\
&= i t_{so} \sum_{<ij>\alpha \beta} \overline{\left( \sigma_x^{\alpha \beta} \cdot d^y_{ij} - \sigma_y^{\alpha \beta} \cdot d^x_{ij} \right)} c^\dagger_{j\beta}c_{i\alpha} \\
&= i t_{so} \sum_{<ij>\alpha \beta} \left( \overline{\sigma_x^{\alpha \beta} \cdot d^y_{ij}} - \overline{\sigma_y^{\alpha \beta} \cdot d^x_{ij}} \right) c^\dagger_{j\beta}c_{i\alpha} \\
&= i t_{so} \sum_{<ij>\alpha \beta} \left(\sigma_x^{\alpha \beta} \cdot d^y_{ij} + \sigma_y^{\alpha \beta} \cdot d^x_{ij} \right) c^\dagger_{j\beta}c_{i\alpha} \\
&= i t_{so} \sum_{<ij>\alpha \beta} \left(-\sigma_x^{\alpha \beta} \cdot d^y_{ji} - \sigma_y^{\alpha \beta} \cdot d^x_{ji} \right) c^\dagger_{j\beta}c_{i\alpha} \\
&= i t_{so} \sum_{<ij>\alpha \beta} \left(-\sigma_x^{\beta \alpha} \cdot d^y_{ji} + \sigma_y^{\beta \alpha} \cdot d^x_{ji} \right) c^\dagger_{j\beta}c_{i\alpha} \\
&= -i t_{so} \sum_{<ij>\alpha \beta} \left(+\sigma_x^{\beta \alpha} \cdot d^y_{ji} - \sigma_y^{\beta \alpha} \cdot d^x_{ji} \right) c^\dagger_{j\beta}c_{i\alpha} \\
&= -h_{ji}^{\beta \alpha}
\end{align}
$$
