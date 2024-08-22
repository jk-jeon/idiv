# `idiv`

This library is a collection of generic algorithms and supporting utilities for performing analysis of Warren-Granlund-Montgomery style optimized integer divisions and other related problems, including formatting/parsing of integers and floating-point numbers. The project has started when I first discovered an algorithm for finding a triple $(k,m,s)$ of integers such that

$$
  \left\lfloor nx \right\rfloor = \left\lfloor \frac{nm + s}{2^{k}} \right\rfloor
$$

holds for all $n = 1,\ \cdots\ ,n_{\max}$, for given real number $x$ and a positive integer $n_{\max}$, which I described in my [blog post](https://jk-jeon.github.io/posts/2023/08/optimal-bounds-integer-division/). Since then, the algorithm has been continuously improved and generalized, and in particular I obtained an algorithm for computing the precise set of $(\xi,\zeta)\in\mathbb{R}^{2}$ satisfying

$$
  \left\lfloor nx + y \right\rfloor = \left\lfloor n\xi + \zeta \right\rfloor
$$

for all $n = n_{\min},\ \cdots\ ,n_{\max}$, for given real numbers $x,y$ and a range $\left\\{n_{\min},\ \cdots\ ,n_{\max}\right\\}$ of consecutive integers. More detailed explanation of what this is all about and how exactly the algorithm works can be found in [this](https://github.com/jk-jeon/idiv/blob/main/docs/pdf/A%20Note%20on%20Floor%20Computation.pdf) informal paper.

# Building and installing

See the [BUILDING](BUILDING.md) document.

# Licensing

All code is licensed under either of

 * Apache License Version 2.0 with LLVM Exceptions ([LICENSE-Apache2-LLVM](LICENSE-Apache2-LLVM) or https://llvm.org/foundation/relicensing/LICENSE.txt) or
 * Boost Software License Version 1.0 ([LICENSE-Boost](LICENSE-Boost) or https://www.boost.org/LICENSE_1_0.txt).

