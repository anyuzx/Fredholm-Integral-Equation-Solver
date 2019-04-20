This is a python module for solving Fredholm integral equation of the first kind.

# Background

A Fredholm integral equation of the first kind is written as,

$$f(x)=\int_{a}^{b}K(x,s)p(s)\mathrm{d}s.$$

The problem is to find $p(s)$, given that $f(x)$ and $K(x,s)$ are known. This equation occurs quite often in many different areas. In many context, both $f(x)$ and $g(s)$ are probability density function (PDF), and $K(x,s)$ is a conditional PDF as $K(x|s)$. There are two constraints on the solution $p(s)$, that is $p(s)\geq 0$ and $\int p(s)\mathrm{d}s = 1$. 

# Discretization-based Method

Here we describe a discretization-based method to solve the Fredholm integral equation. The integral equation is replaced by a summation approximately,

$$f(\boldsymbol{x})$$
