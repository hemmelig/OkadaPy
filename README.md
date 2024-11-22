<p align="center">
  <!-- DOI -->
  <a href="https://doi.org/10.5281/zenodo.">
    <img src="https://zenodo.org/badge/DOI/10.5281/zenodo..svg" />
  </a>
  <!-- ReadTheDocs -->
  <a href="https://okadapy.readthedocs.io/en/latest">
    <img src="https://readthedocs.org/projects/okadapy/badge/?version=latest" />
  </a>
  <!-- License -->
  <a href="https://www.gnu.org/licenses/gpl-3.0">
    <img src="https://img.shields.io/badge/License-GPLv3-blue.svg" />
  </a>
</p>

<p align="center">
  <a href="https://okadapy.readthedocs.io/en/latest/index.html">OkadaPy</a> is a cross-platform Python package for the analytical computation of displacement, strain, and stress fields.</a>.
</p>

<p align="center">
  <img src="https://raw.githubusercontent.com/hemmelig/OkadaPy/main/docs/img/OkadaPyLogo.png" width="80%" />
</p>

Key features
------------
OkadaPy provides an efficient implementation of the equations set out in Okada, 1992, for the computation of the displacement induced by either a finite rectangular source or an infinitesimal point source in an isotropic elastic halfspace. It is a free and open-source alternative to a number of functions available in the MATLAB package Coulomb, primarily:

- computing displacements induced by the aforementioned sources;
- computing the corresponding strain and stress tensors from these displacements;
- and visualising these fields.

Documentation
-------------
Documentation for OkadaPy is hosted [here](https://okadapy.readthedocs.io/en/latest/index.html).

Installation
------------
OkadaPy requires Python version 3.11 and above. The package will be made available through the Python Package Index (PyPI), at which time it will be possible to install via `pip`:

```console
pip install okadapy
```

Until then, it can be installed from source, providing a suitable C compiler is available:

```console
git clone https://github.com/hemmelig/OkadaPy
cd OkadaPy
pip install .
```

For further information regarding installation—including virtual environment management and installation from source—please consult [our documentation](https://okadapy.readthedocs.io/en/latest/installation.html).

Usage
-----
Tutorial material is in the works for the package.

I previously wrote a wrapper package for the MATLAB package Coulomb, which made use of the MATLAB engine for Python. This is a replacement for this tool, built around a core library of compiled C code.

Citation
--------
If you use OkadaPy in your work, please cite the following:

OkadaPy Developers (2024). OkadaPy: v0.0.1 (v0.0.1). Zenodo. https://doi.org/10.5281/zenodo.

Toda, S., Stein, R. S., Sevilgen, V., & Lin, J. (2011). Coulomb 3.3 Graphic-rich deformation and stress-change software for earthquake, tectonic, and volcano research and teaching—user guide. US Geological Survey open-file report, 1060(2011), 63.

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
