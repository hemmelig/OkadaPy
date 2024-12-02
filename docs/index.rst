.. figure:: img/OkadaPyLogo.png
   :figwidth: 70 %
   :width: 90%
   :align: center
   :alt: OkadaPy: an efficient implementation of the Okada, 1992, equations.

OkadaPy
=======

:mod:`OkadaPy` is a Python package that provides an efficient implementation of the analytical equations for deformation in an isotropic elastic half-space, as set out in Okada, 1992.

These equations describe the displacement induced by either a finite rectangular source or an infinitesimal point source in an isotropic elastic half-space. It is a free and open-source alternative to a number of functions available in the MATLAB package Coulomb3.4, primarily:

* computing displacements induced by the aforementioned sources;
* computing the corresponding strain and stress tensors from these displacements;
* and visualising these fields.

The source code for the project is hosted on |github|.

This package is written by Conor A. Bacon, and is distributed under
the GPLv3 License, Copyright Conor A. Bacon, 2024.


.. |github| raw:: html

    <a href="https://github.com/hemmelig/OkadaPy" target="_blank">GitHub</a>

Citation
--------
If you use this package in your work, please cite the following:

OkadaPy Developers (2024). OkadaPy: v0.0.1 (v0.0.1). Zenodo. https://doi.org/10.5281/zenodo.

Toda, S., Stein, R. S., Sevilgen, V., & Lin, J. (2011). Coulomb 3.3 Graphic-rich deformation and stress-change software for earthquake, tectonic, and volcano research and teaching—user guide. US Geological Survey open-file report, 1060(2011), 63.

as well as the relevant version of the source code on `Zenodo <https://doi.org/10.5281/zenodo.>`_.

Contributing to OkadaPy
-----------------------
Contributions to OkadaPy are welcomed. The first stop should be to reach out, either directly or—preferably—via the GitHub Issues panel, to discuss the proposed changes. Next, simply fork the OkadaPy repository, install the package with the developer options enabled (that is, using `pip install okadapy[dev]`), make your changes/add your new contribution, then make a [pull request](https://help.github.com/articles/about-pull-requests/). All contributors to OkadaPy will be listed as authors on the releases.

Bug reports, suggestions for new features and enhancements, and even links to projects that have made use of OkadaPy are most welcome.

See our [contributions page](https://github.com/hemmelig/OkadaPy/blob/main/.github/CONTRIBUTING.md) for more information.

Contact
-------
Any comments/questions can be directed to:
* **Conor Bacon** - cbacon [ at ] ldeo.columbia.edu

License
-------
OkadaPy is **free** and **open source**, distributed under the GPLv3 License. Please see the [LICENSE](LICENSE) file for a complete description of the rights and freedoms that this provides the user.

Contents:
---------

.. toctree::
   :numbered:
   :maxdepth: 1

   installation
   tutorials
   sourcecode
