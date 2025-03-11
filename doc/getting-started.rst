
###############
Getting started
###############

************
Installation
************

Pre-built binaries
==================

Executables for Linux and Windows are provided in the `Releases page <https://gitlab.com/lheea/Nemoh/-/releases>`__. They can be used directly without the need to do the compilation procedure described in the next subsection. Simply download the executables and consider :ref:`getting-started:Adding Nemoh to the system ``PATH```.

Compilation
===========

Compilation should only be required if you are in one of the following situations:

- No pre-built binaries are available for your platform on the `Releases page <https://gitlab.com/lheea/Nemoh/-/releases>`__, or

- The pre-built binaries for your platform do not work, or

- You modified the code

Otherwise, we recommend using the :ref:`getting-started:pre-built binaries`.

Dependencies
------------

The following tools are required to build Nemoh:

-  A Fortran compiler. The code has been tested using:

   -  GNU's `gfortran <https://gcc.gnu.org/wiki/GFortran>`__,

   -  Intel's `ifort <https://www.intel.com/content/www/us/en/developer/tools/oneapi/fortran-compiler.html#gs.jik1s6>`__.

-  `CMake <https://cmake.org/>`__, a cross-platform tool for building and testing the software package.


Additionally, the following external libraries are used by Nemoh:

-  BLAS, https://netlib.org/blas/

-  LAPACK, https://netlib.org/lapack/

Any implementation of BLAS/LAPACK should work, but only Netlib's implementation has been tested.

Building the executables
------------------------

Compile all Nemoh executables using CMake (from the root of the repository):

.. code:: bash

   cmake -S. -Bbuild
   cmake --build build

The resulting executables are in the ``bin/`` directory. To compile only one of the executables, use the ``--target`` option of CMake. The available targets are:

-  for Nemoh: ``mesh``, ``preProc``, ``hydrosCal``, ``solver``, ``postProc``

-  for Nemoh QTF: ``QTFpreProc``, ``QTFsolver``, ``QTFpostProc``

The choice of the compiler is left to CMake, but can be overridden by setting the ``CMAKE_Fortran_COMPILER`` at the configuration step, e.g.:

.. code:: bash

   cmake -S. -Bbuild -DCMAKE_Fortran_COMPILER=gfortran

Checking the compilation
------------------------

After building, the quick tests can be run from the ``build/`` directory:

.. code:: bash

   ctest -V -j <N_concurrent>

Where ``<N_ concurrent>`` is the number of simultaneous workers (processes). The tests can be restricted using their labels and the ``-L`` option of ``ctest``:

.. code:: bash

   ctest -V -j <N_concurrent> -L <label>

Where label is one of the following:

-  ``Nemoh1``: only the non-QTF test cases

-  ``PREPROC``: only the pre-processing operations

-  ``SOLVER``: only the solving operations (depend on the pre-processing tests)

-  ``POSTPROC``: only the post-processing operations (depend on the pre-processing and solving tests)

-  ``Nemoh2``: only the QTF test cases

-  ``QTF``: only the computation of the QTF (depend on the prior non-QTF Nemoh computation)

Tests with unsatisfied requirements will fail.

Adding Nemoh to the system ``PATH``
===================================

Nemoh does not come with a proper package or installer for system. Instead, the executables are portable and may be placed and executed anywhere.
The most straightforward way of running Nemoh is by refering to its executables by their full path on the system, but that is not very convenient.

To be able to call the Nemoh programs from anywhere on your system without providing their full path, one needs to add their location to the system's ``PATH``.
The procedure depends on the operating system. Examples are provided here for Linux and Windows, but keep in mind that it may differ depending on your version or distribution.

In the following, ``<nemoh-programs-path>`` refers to the location where you keep the Nemoh executables.
If you just compiled Nemoh and want to keep it as is, that would be the ``bin/`` directory where you compiled Nemoh.
Alternatively, e.g. if you just downloaded the :ref:`getting-started:Pre-built binaries`, you can copy the executables anywhere that suits you before setting the system's ``PATH``.

Linux
-----

The following assumes you use ``bash`` as your terminal.

For the current terminal session only
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

#. Add Nemoh programs' location to the ``PATH`` environment variable:

   .. code:: bash

      export PATH=$PATH:<nemoh-programs-path>

#. Check that you system finds the programs:

   .. code:: bash

      which solver

   Should return the full location of the ``solver`` executable.

For all future terminal sessions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

#. Set your ``bash`` profile:

   .. code:: bash

      echo 'export PATH=$PATH:<nemoh-programs-path>' >> ~/.bashrc

#. Restart the terminal

#. Check that you system finds the programs:

   .. code:: bash

      which solver

   Should return the full location of the ``solver`` executable.

Windows
-------

The following assumes you use ``cmd`` (standard command prompt) as your terminal.

For the current terminal session only
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

#. Add Nemoh programs' location to the ``path`` environment variable:

   .. code:: winbatch

      set path=%path%;C:<nemoh-programs-path>

#. Check that you system finds the programs:

   .. code:: winbatch

      where solver.exe

   Should return the full location of the ``solver.exe`` executable.

For all future terminal sessions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

#. Permanently add Nemoh programs' location to the ``path`` environment variable:

   .. code:: winbatch

      setx path "%path%;<nemoh-programs-path>"

#. Restart the terminal

#.  Check that you system finds the programs in a new terminal:

   .. code:: winbatch

      where solver.exe

   Should return the full location of the ``solver.exe`` executable.

.. note:: If you're using PowerShell, use ``which`` instead of ``where``.


*****
Usage
*****

Executable files
================

Nemoh is made of the following executables files (programs):

-  Nemoh1: ``mesh``, ``preProc``, ``hydrosCal``, ``solver``, ``postProc``,

-  Nemoh2: ``QTFpreProc``, ``QTFsolver``, ``QTFpostProc``.


The programs of Nemoh1 and Nemoh2 normally have to be executed following this order.
See :ref:`description:Software features` for details about each program's purpose, and :ref:`input-output:Input/Output` for what files are used and outputted by each program.

Any Nemoh program (``<program>``) can be executed directly in the project directory as:

.. code:: bash

   <program>

Alternatively, one can specify the project directory (``<project-dir>``) as an argument:

.. code:: bash

   <program> <project-dir>

The project directory should contain the input files described in :ref:`input-output:Input/Output` (``Nemoh.cal``, ``mesh/``, ``Mechanics/``, etc...).

.. note:: A Matlab wrapper is provided to use the programs from a Matlab environment. More details are provided in :ref:`getting-started:supporting matlab files`.


Running the test cases
======================

To simplify the procedure for Linux platforms (or Windows with Make), a ``Makefile`` is provided in the ``TestCases/`` directory. It is then possible to run the Nemoh1 test cases by executing the following commands in a Terminal (each line being a test case):

.. code:: bash

   make run_1_cylinder
   make run_2_2Bodies
   make run_3_nonsymmetrical
   make run_4_postprocessing
   make run_5_quicktest
   make run_6_box_coarsemesh
   make run_7_Solvers_Check_OC3
   make run_8a_Cylinder_irregfreq

For the QTF test cases, the following commands can be used:

.. code:: bash

   make run_8b_QTF_Cylinder
   make run_9_oc4_semisub
   make run_10a_softwind
   make run_10b_softwind_FS
   make run_11_QTF_OC3_Hywind.

Commands to clean the test cases are also available to clean all the output files. They can apply either to a specific tests case, *e.g.*

.. code:: bash

   make clean_1_cylinder

Or to remove a range of test cases

.. code:: bash

   make clean_all_testsNemoh1
   make clean_all_testsNemoh2
   make clean_all_tests

The description and the benchmark results of those test cases are described in :ref:`test-cases:test cases`.


Supporting Matlab files
=======================

.. warning::

   The standard way of running Nemoh is via the command line.
   The Nemoh Matlab wrapper is only provided for convenience and quick post-processing.
   It is not actively developed, and may or may not be updated following future evolutions in Nemoh's main code.
   As-is, it only provides access to a subset of features.
   In particular, do not rely on the Maltab wrapper for mesh generation on advanced geometries.
   **Issues opened in GitLab regarding the Matlab wrapper will be closed with a link to this section.**

The following directories, containing a set of Matlab functions, are provided in ``matlabRoutines/``:

-  ``NemohWrapper``: This is for running Nemoh executables in MATLAB environment.

-  ``GMSHconverter``: There are two codes, first, for converting body mesh file output from GMSH to Nemoh, DIODORE and HydroSTAR formats and second, for converting free-surface mesh file output from GMSH to Nemoh and HydroSTAR formats.

-  ``postproc_testcases``: There are two main codes for plotting results from Nemoh and HydroSTAR. First, for plotting hydrodynamic coefficients results and second for plotting QTF results. This code can be executed after all data in one specific test cases are obtained.
