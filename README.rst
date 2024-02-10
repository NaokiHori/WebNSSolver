#############
Web NS Solver
#############

|License|_ |LastCommit|_ |Deployment|_

.. |License| image:: https://img.shields.io/github/license/NaokiHori/WebNSSolver
.. _License: https://opensource.org/license/MIT

.. |LastCommit| image:: https://img.shields.io/github/last-commit/NaokiHori/WebNSSolver/main
.. _LastCommit: https://github.com/NaokiHori/WebNSSolver/commits/main

.. |Deployment| image:: https://github.com/NaokiHori/WebNSSolver/actions/workflows/deployment.yml/badge.svg?branch=main
.. _Deployment: https://github.com/NaokiHori/WebNSSolver/actions/workflows/deployment.yml

Navier-Stokes solver in Browsers.

.. image:: https://naokihori.github.io/WebNSSolver/thumbnail.jpg
   :target: https://naokihori.github.io/WebNSSolver/index.html
   :width: 800

*****************************************************************
`Main page <https://naokihori.github.io/WebNSSolver/index.html>`_
*****************************************************************

=============
Visualisation
=============

Click the canvas to change the scalar field to be drawn.
For now 7 options are available:

* Temperature

* X velocity

* Y velocity

* Velocity magnitude

* Vorticity

* Q-value

* Tracer particles

==========
Parameters
==========

By default the Rayleigh and the Prantdl numbers are set to ``8.5e7`` and ``4.4e0``, respectively.
Two URL parameters ``ra`` and ``pr`` are available to change them, e.g.: ``https://naokihori.github.io/WebNSSolver/index.html?ra=1.0e8&pr=0.25``.

The minimum / maximum Rayleigh numbers are ``1.0e+4`` / ``4.0e+8``.
The minimum / maximum Prantdl numbers are ``1.0e-1`` / ``1.0e+1``.

*************
Documentation
*************

Implementation details are documented `here <https://naokihori.github.io/WebNSSolver/doc/web_ns_solver/index.html>`_.

