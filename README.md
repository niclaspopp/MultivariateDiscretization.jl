# MultivariateDiscretization

[![Build Status](https://travis-ci.com/niclaspopp/MultivariateDiscretization.jl.svg?branch=master)](https://travis-ci.com/niclaspopp/MultivariateDiscretization.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/niclaspopp/MultivariateDiscretization.jl?svg=true)](https://ci.appveyor.com/project/niclaspopp/MultivariateDiscretization-jl)
[![Coverage](https://codecov.io/gh/niclaspopp/MultivariateDiscretization.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/niclaspopp/MultivariateDiscretization.jl)
[![Stable]](https://github.com/niclaspopp/MultivariateDiscretization.jl/blob/master/MultivariateDiscretization/Doc%20v0.1.0.ipynb)

MultivariateDiscretization.jl is a Julia package that can be used for mulitvariate, interaction preserving discretization. Currently the following methods are implemented:

+ Uniform width binning
+ Bayesian Blocks using [Discretizers.jl](https://github.com/Tchanders/Discretizers.jl).
+ [Correlation preserving discretization](https://ieeexplore.ieee.org/document/1490525)
+ [Interaction preserving discretization](https://link.springer.com/article/10.1007/s10618-014-0350-5)

For testing purposes a function to evaluate the [Cumulative Jensen Shannon divergence](http://eda.mmci.uni-saarland.de/prj/cjs/) is included as well.

Documentation can be found [here](https://github.com/niclaspopp/MultivariateDiscretization.jl/blob/master/MultivariateDiscretization/Doc%20v0.1.0.ipynb).
