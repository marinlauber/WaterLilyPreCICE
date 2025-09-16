## 3d-od Coupling of an elastic sphere to a Windkessel model

### The Windkessel model

$$
 \begin{cases}\begin{split}
    \frac{dV_{LV}^{0D}}{dt} &= Q_{MV} - Q_{AO}\\
    \frac{dP_{AO}}{dt} &= \frac{Q_{AO}}{C} - \frac{P_{AO}}{RC}
\end{split}\end{cases}
$$

$$
Q_{MV} = \begin{cases}(P_\text{fill}-P_{LV})/R_{MV}&\text{ if }P_\text{Fill} \ge P_{LV}\\
\epsilon &\text{ else }\end{cases}
$$
$$
Q_{AO} = \begin{cases}(P_{LV} - P_{AO})/R_{AV}&\text{ if } P_{LV} \ge P_{AO}\\
\epsilon &\text{ else }\end{cases}
$$

$$
P_{LV} = \frac{V_{LV}^{0D} - V_{0}}{C_{LV}}
$$