# Third Damper Transfer Function Derivation
(1) $ \frac{F_{m_s}+m_s g}{V_{in} b_1}(s) = \frac{1}{b_1}\cdot\frac{F_{m_s}+Se_{m_sg}}{V_{in}} = \frac{1}{b_1}\cdot\Big( \frac{F_{m_s}}{V_{in}} + \frac{Se_{m_sg}}{V_{in}}\Big)$
Since $F_{m_s} = \dot{p_{m_s}}$, after a Laplace transform $F_{m_s}(s) = p_{m_s}(s)*s$, 
Substituting into (1):
$\frac{1}{b_1}\cdot\Big( \frac{p_{m_s}}{V_{in}}\cdot s + \frac{Se_{m_sg}}{V_{in}}\Big)$
Then to rewrite $\frac{Se_{m_sg}}{V_{in}}$ in terms of state-variables transfer functions, the transfer function is decomposed into $(\frac{p_{m_s}}{Se_{m_sg}})^{-1}\frac{p_{m_s}}{V_{in}} = \frac{Se_{m_sg}}{p_{m_s}}\frac{p_{m_s}}{V_{in}} = \frac{Se_{m_sg}}{V_{in}}$
Therefore the Transfer Function can be rewritten in terms of state variables with:
$\frac{F_{m_s}+m_s g}{V_{in} b_1}(s) = \frac{1}{b_1}*\Big(\frac{p_{m_s}}{V_{in}}(s)\cdot s +  (\frac{p_{m_s}}{Se_{m_sg}}(s))^{-1}\frac{p_{m_s}}{V_{in}}(s)\Big)$