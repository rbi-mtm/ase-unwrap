ase-unwrap
===========================

A simple tool for unwrapping `ASE <https://wiki.fysik.dtu.dk/ase/>`_-readable atomic structures.

Installation
----------------------------------

The following installation procedure was tested:

.. code:: shell

        conda create -n ase-unwrap python=3.9
        conda activate ase-unwrap
        python -m pip install -e .

Usage
----------------------------------

Command-line usage:

.. code:: shell

        ase-unwrap <ASE_readable_file>

The code relies on `ASE neighborlists <https://wiki.fysik.dtu.dk/ase/_modules/ase/neighborlist.html>`_ to determine the interatomic connectivity. Bonds are found automatically using ``natural_cutoffs`` (covalent radii). The default radii are usually a reasonable estimation, but a multiplier for all cutoffs may be introduced as follows:

.. code:: shell

        ase-unwrap <ASE_readable_file> <cutoff_multiplier>
        ase-unwrap <ASE_readable_file> 0.9
        ase-unwrap <ASE_readable_file> 1.1

The output file is automatically generated. If you wish to use a non-default output name:

.. code:: shell

        ase-unwrap <ASE_readable_file> <cutoff_multiplier> <output_filename>
        ase-unwrap <ASE_readable_file> 1.0 my_filename.cif

The ``src/ase_unwrap/unwrap.py`` module may also be run as a Python script:

.. code:: shell

        python unwrap.py <ASE_readable_file> <cutoff_multiplier> <output_filename>

or used as a library:

.. code-block:: python

        import ase.io
        from ase_unwrap.unwrap import unwrap

        atoms = read('unwrapped_structure.traj')

        unwrap(atoms)
