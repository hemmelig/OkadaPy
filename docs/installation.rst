Installation
============
:mod:`OkadaPy` is a predominantly Python package with some routines written and optimised in C. These are built and linked to OkadaPy at installation - if installing from source you will need to ensure that there is a suitable compiler available (see :ref:`C compilers`).

However, most users can bypass this step by installing OkadaPy using ``pip``.

Supported operating systems
---------------------------
OkadaPy was developed and tested on macOS 12.5.1 (M1 chip), with the intention of being "platform agnostic". As of November 2024, the package has been successfully built and run on:

- Ubuntu 16.04/18.04/20.04/22.04
- Red Hat Enterprise Linux
- Debian
- Windows 10
- macOS Monterey 12 (including M1)

Prerequisites
-------------
OkadaPy supports Python 3.11 or newer (3.11/3.12). We recommend using `uv` as a package manager and virtual environment system to isolate and install the specific dependencies of OkadaPy.

Instructions for downloading and installing `uv` can be found `here <https://github.com/astral-sh/uv>`_.

Installation via `pip`
----------------------
The simplest way to get a working copy of OkadaPy is to install it from the Python Package Index (PyPI) using ``pip`` (the Python package installer).

To do this you first need to set up an enivironment. We recommend creating a minimal environment initially:

.. code-block:: bash
    
    uv venv --python=3.12
    source .venv/bin/activate

All other dependencies will be handled during the installation of OkadaPy. After activating your environment, type the following command into terminal:

.. code-block:: bash
    
    uv pip install okadapy

This will install OkadaPy **and** its explicit dependencies!

.. note:: Installing the package this way will not provide you with the examples. These can be retrieved directly from the GitHub repository (see :ref:`testing your installation`).

The full list of dependencies is:

- matplotlib
- numpy
- pandas
- pyproj >= 2.5
- pygmt (optional plotting backend)

.. note:: We are currently not pinning the version of any dependencies. We aim to keep on top of any new syntax changes etc. as new versions of these packages are released - but please submit an issue if you come across something!

If you want to explore the example notebooks, you will also need to install `ipython` and `jupyter`. This can be done with conda (making sure your environment is still activated) as:

.. code-block:: bash

    uv pip install jupyter ipython

Other installation methods
--------------------------
From source
***********

.. note:: In order to install from source, you will need an accessible C compiler, such as `gcc` or `clang` (see :ref:`C compilers`).

`Clone the repository <https://help.github.com/en/github/creating-cloning-and-archiving-repositories/cloning-a-repository>`_ from our `GitHub <https://github.com/hemmelig/OkadaPy>`_ (note: you will need ``git`` installed on your system), or alternatively download the source code directly through the GitHub web interface. Once you have a local copy, navigate to the new ``OkadaPy`` directory.

As with the installation from PyPI, you should set up a virtual environment:

.. code-block:: bash
    
    uv venv --python=3.12
    source .venv/bin/activate

Since building from source means the C routines included with the package need to be compiled during installation, you will need an appropriate C compiler, such as `gcc` (Linux) or `clang` (macOS). If you are using Linux and `uv`, you will need to first set an environment variable specifying GCC should be used. This is necessary because the `uv` binaries are built with `clang`.

.. code-block:: bash
    export CC=gcc
    uv pip install .

You can optionally pass a ``-e`` argument to install the package in 'editable' mode.

You should now be able to import :mod:`okada` within a Python session:

.. warning:: You should try this import in any directory that is *not* the root of the git repository (i.e. ``OkadaPy/``. Here, the local ``okada`` directory will override the version of OkadaPy installed in your environment site-packages!

.. code-block:: bash
    
    cd  # Moving out of OkadaPy directory - see warning above!
    python
    >>> import okada
    >>> okada.__version__

If successful, this should return '|Version|'.

Testing your installation
-------------------------
In order to test your installation, you will need to have cloned the GitHub repository (see :ref:`installation from source <From Source>`). This will ensure you have all of the required benchmarked data (which is not included in pip/conda installs).

To run the tests, navigate to ``OkadaPy/tests`` and run the benchmark testing script. This compares the outputs of OkadaPy with those computed for an identical input file with Coulomb3.4:

.. code-block:: bash

    python run_benchmark_tests.py

If your installation is working as intended, this should execute with no failures.

C compilers
-----------
In order to install and use OkadaPy from source, you will need a C compiler.

If you already have a suitable compiler (e.g. `gcc`, `MSVC`, `clang`) at the OS level, then you can proceed with installation. If this fails, then read on for tips to overcome common issues!

Checking for a C compiler
*************************
On Linux or macOS, to check if you have a C compiler, open a terminal and type:

.. code-block:: bash
    
    which gcc
    gcc --version

If a compiler is present, the first command will return ``/usr/bin/gcc``. However, this does not guarantee it is present! The second command will confirm this.

On **Linux** the second command should output something like:

.. code-block:: console

    gcc (Ubuntu 11.3.0-1ubuntu1~22.04) 11.3.0
    Copyright (C) 2021 Free Software Foundation, Inc.
    This is free software; see the source for copying conditions. There is NO
    warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

As long as the version is relatively recent (version 9 or later), you should be good to go. To additionally confirm that you have a C++ compiler installed, type:

.. code-block:: bash
    
    which g++
    g++ --version

For which you should obtain a similar result.

On **macOS** it will be obvious if the compiler is not actually installed -- you will be faced with a prompt to install the Xcode Command Line Tools. You can go ahead and install this (press ``Install`` and wait for the process to complete). If these are already installed, the second command should output something like:

.. code-block:: console

    Configured with: --prefix=/Library/Developer/CommandLineTools/usr --with-gxx-include-dir=/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include/c++/4.2.1
    Apple clang version 12.0.5 (clang-1205.0.22.11)
    Target: x86_64-apple-darwin20.5.0
    Thread model: posix
    InstalledDir: /Library/Developer/CommandLineTools/usr/bin

.. warning:: Even if `clang` is installed, it may not have all necessary libraries included. See :ref:`OpenMP on macOS`.

Note that this indicates that the system compiler is ``clang``, and that the accompanying C++ compiler is also installed. These are all supplied as part of the Xcode Command Line Tools (see e.g. `here <https://mac.install.guide/commandlinetools/index.html>`_ for a rundown).

If you do not have a compiler, or to be sure, we provide a simple guide for :ref:`Linux`, :ref:`macOS` and :ref:`Windows` operating systems below.

Linux
*****
If you have root access, the simplest route is to install `gcc` and `gxx` at system-level. You should search for the correct way to do this for your Linux Distribution. For example, on Ubuntu you would type:

.. code-block:: console

    sudo apt-get install build-essential

This includes `gcc`, `g++` as well as `make`. The commands will differ on other distros (CentOS, Red Hat, etc.).

Alternatively, you can install `gcc` and `g++` `through conda <https://docs.conda.io/projects/conda-build/en/latest/resources/compiler-tools.html>`_. Make sure you have activated your environment, then type:

.. code-block:: bash
    
    conda install -c conda-forge gcc_linux-64 gxx_linux-64

You can test this was successful with the same procedure detailed :ref:`above<Checking for a C compiler>`. Once installed, you can proceed with the OkadaPy :ref:`installation from source <Other installation methods>`.

macOS
*****
By default, there is no C compiler included with macOS. If you have previously installed the Xcode Command Line Tools (via the web or the App Store), the `clang` compiler will be installed. However, this may not include all necessary libraries to install OkadaPy (see :ref:`OpenMP on macOS`).

Whether you already have Xcode installed or not, there are two options to install what is required: the user can either install all dependencies :ref:`through conda <conda>` - noting that they will only be available in that specific environment - or using `HomeBrew <https://brew.sh/>`_. We generally recommend using conda, unless the user is already familiar with brew (in which case, see :ref:`brew`).

OpenMP on macOS
+++++++++++++++
The default C compiler on macOS does not include support for OpenMP. This will result in the following error during installation from source:

.. code-block:: console

    ld: library not found for -lomp
    clang: error: linker command failed with exit code 1 (use -v to see invocation)
    error: command '/usr/bin/clang' failed with exit code 1

As above, this can either be solved with :ref:`conda` or :ref:`brew`.

conda
+++++
First create and/or activate your environment:

.. code-block:: bash

    conda create -n okada python=3.12  # if not already created
    conda activate okada  # replace with alternative environment name if desired

Then use conda to install the compiler (along with the OpenMP libraries). **Note the syntax is different if your machine is running on an Apple Silicon (M1, M2, etc.) chip**:

.. code-block:: bash

    conda install -c conda-forge llvm-openmp clang_osx-64 clangxx_osx-64  # Intel chip
    conda install -c conda-forge llvm-openmp clang_osx-arm64 clangxx_osx-arm64  # Apple Silicon chip (M1, M2 etc.)

.. note:: If you did not already have Xcode Command Line Tools installed, you will be prompted to install them now. Click ``Install`` and wait for installation to complete.

You should now open a fresh terminal, and activate your environment. To test the installation was successful, type:

.. code-block:: bash

    echo $CC
    $CC --version

This should return something like:

.. code-block:: console

    echo $CC
    x86_64-apple-darwin13.4.0-clang
    $CC --version
    clang version 14.0.6
    Target: x86_64-apple-darwin13.4.0
    Thread model: posix
    InstalledDir: /Users/user/miniconda3/envs/okada/bin

You can proceed with the OkadaPy :ref:`installation from source <Other installation methods>`.

brew
++++
If brew is not already installed (check with ``which brew``), follow the instructions on the `HomeBrew frontpage <https://brew.sh/>`_. This will offer to install the Xcode Command Line Tools if they are not already present (press 'ENTER' or 'Y' to accept this suggestion).

You can then proceed to install the OpenMP libraries with brew:

.. code-block:: bash
    
    brew install libomp

You can safely ignore the warning about explicitly adding the relevant LDFLAGS etc. - this is already handled in the OkadaPy ``setup.py`` script.

You can proceed with the OkadaPy :ref:`installation from source <Other installation methods>`.

*Legacy*: brew gcc
++++++++++++++++++
Alternatively, you can use the `gcc` compiler to install OkadaPy. As with `clang`, we recommend installing GCC through ``Homebrew``. First, check if you already have `gcc` installed, with:

.. code-block:: bash

    which gcc

If this doesn't return anything, continue to installing `gcc`. If this returns the path to a gcc executable (e.g. `/usr/bin/gcc`), then you should check the version, with:

.. code-block:: bash

    gcc --version

If the version string includes `Apple clang`, or is a version number lower than 9, you should proceed to install with ``Homebrew``:

.. code-block:: bash
    
    brew install gcc
    brew link gcc

Note that the ``brew link`` command should add ``gcc`` to your path, but might not succeed if a previous ``gcc`` install was present. To test this, type:

.. code-block:: bash

    which gcc
    gcc --version

If the linking was successful, this should point to a new gcc executable, and the version string should contain ``gcc (Homebrew GCC 9.4.0) 9.4.0`` or similar. If not, you will need to manually link the new ``gcc`` executable. To do this, find the path to your new ``gcc``` installation with:

.. code-block:: bash

    brew --prefix gcc

Then create a symlink to this executable:

.. code-block:: bash

    ln -s /usr/local/bin/gcc /path/to/brew/gcc

Where ``/path/to/brew/gcc`` is the path returned by the ``brew --prefix`` command.

Finally, test this has worked by repeating the check from above:

.. code-block:: bash

    which gcc
    gcc --version

This should now return the ``Homebrew`` ``gcc`` version string. If not, please get in touch and we will try to help if we can...

Windows
*******
Compilation and linking of the C extensions has been successful using the Microsoft Visual C++ (MSVC) build tools. 

We strongly recommend that you download and install these tools in order to use OkadaPy. You can either install Visual Studio in its entirety, or just the Build Tools - `available here <https://visualstudio.microsoft.com/downloads/>`_.

You will need to restart your computer once the installation process has completed. We recommend using the anaconda command line interface (unix shell-like) to install OkadaPy over command prompt.

Once installed, you can proceed with the OkadaPy :ref:`installation from source <Other installation methods>`.
