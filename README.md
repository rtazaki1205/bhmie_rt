# Introduction

This package provides a light scattering solution of a homogeneous sphere by Lorenz-Mie theory.
This code is based on BHMIE ([Bohren & Huffman 1983](https://ui.adsabs.harvard.edu/abs/1983asls.book.....B)) 
implemented by B. T. Draine (https://www.astro.princeton.edu/~draine/scattering.html).
However, this code adopts the continued fraction method to find the logarithmic 
derivative ([Lentz et al. 1976](https://ui.adsabs.harvard.edu/abs/1976ApOpt..15..668L)) unlike the standard BHMIE.

As a benchmark test, the results obtained by this code are compared with the ones 
obtained by Wiscombe's MIEV0 code ([Wiscombe 1980](https://ui.adsabs.harvard.edu/abs/1980ApOpt..19.1505W)).

# Terms of use

The codes are distributed under the [MITlicense](https://opensource.org/licenses/MIT) and can be used, changed
and redistributed freely. If you use this package to publish papers, please cite the relevant papers.

# History

 - Initial release
