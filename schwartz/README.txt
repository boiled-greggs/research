# Making a PythTB input file

Will need:
    - lattice parameters
    - orbitals in crystal coordinates (supplied by .vasp file, hopefully)
    - hopping parameters

Later:
    - kpath information (need to learn more about this)
    - making model and calculating bands, then displaying
        + I think I just need to copy from the example files to get the right info

Notes:

begin by creating the input file, then write the first few lines
    - these lines can be copied from existing ones

then, need to parse .vasp file into lattice coordinates and orbitals:
    - first line says "PDB file" and next line has some number, not sure what (*should learn this*)
        + **Third Line** begins the lattice vectors
        + there are three lattice vectors
        + parse and write to file
    - next line has atomic information (type of atoms)
    - next line has number of unique atoms
    - next line says "Direct"
    - THEN, we see the orbital positions in crystal coordinates,
        if they were exported in the most helpful way
        + parse this and then write to input file

next, hopping needs to be defined
    - inside the cell, this is relatively straightforward
    - edge cases:
        + will get some hops within cell from existing algorithm
        + must get more hopping info from shifting cell and seeing if any
            shifted atoms are within hopping distance


From here, I need to find out how to make a path for the calculation to be done.
I really don't know how to do this so I'll ask Aldo what he thinks.
