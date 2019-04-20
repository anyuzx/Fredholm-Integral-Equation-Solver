This is a python module for solving Fredholm integral equation of the first kind.

# Background

A Fredholm integral equation of the first kind is written as,

$$f(x)=\int_{a}^{b}K(x,s)p(s)\mathrm{d}s.$$

The problem is to find $p(s)$, given that $f(x)$ and $K(x,s)$ are known. This equation occurs quite often in many different areas. In many context, both $f(x)$ and $g(s)$ are probability density function (PDF), and $K(x,s)$ is a conditional PDF as $K(x|s)$. There are two constraints on the solution $p(s)$, that is $p(s)\geq 0$ and $\int p(s)\mathrm{d}s = 1$. 

# Discretization-based Method

Here we describe a discretization-based method to solve the Fredholm integral equation. The integral equation is approximately replaced by a Riemann summation over grids,

$$f(x_i)=\sum_j \Delta_s K(x_i, s_j) p(s_j)$$

where $\Delta_s$ is the grid size along the dimension $s$ and $x_i$, $s_j$ are the grid points and $i$ and $j$ indicate their indices. When grid size $\Delta_s\to0$, the summation converges to the true integral. It is more convenient to write it in the matrix form,

$$\boldsymbol{f} = \Delta_s \boldsymbol{K} \boldsymbol{p}$$

