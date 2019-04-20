This is a python module for solving Fredholm integral equation of the first kind.

# Background

A Fredholm integral equation of the first kind is written as,

$$f(x)=\int_{a}^{b}K(x,s)p(s)\mathrm{d}s.$$

The problem is to find $p(s)$, given that $f(x)$ and $K(x,s)$ are known. This equation occurs quite often in many different areas. In many cases, $f(x)$ is a probability density function (PDF) and $K(x,s)$ is a conditional PDF as $K(x|s)$. Thus $p(s)$ must also be a PDF, which impose two constraints on the solution $p(s)$, that is $p(s)\geq 0$ and $\int p(s)\mathrm{d}s = 1$.
