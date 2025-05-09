############
Introduction
############

This document is the Manual of Nemoh Software in version v3.0 Released on 2nd December 2022. It serves as a guide for using, installing and running the software.

Nemoh in its original version is an open-source potential flow boundary element solver for computing first-order hydrodynamic coefficients in the frequency domain. Nemoh v3.0 includes the treatment of irregular frequencies, has an extended module to post-process the first-order hydrodynamic results, and compute complete Quadratic Transfer Functions (QTFs).
The new developments/extensions in this Nemoh v3.0 are summarised as follows.

-  GNU General Public v3 license is used instead of Apache.

-  Newly written source codes for better readability and coding practice.

-  Green function uses finer interpolation points.

-  The influence coefficient is constructed by applying Gauss-quadrature integration over a panel.

-  Irregular frequency removal method is available.

-  Two new options of the linear system solvers, LU decomposition and GMRES iterative solver, are available for enhanced computational efficiency.

-  Full difference- and sum-frequencies QTF module.
