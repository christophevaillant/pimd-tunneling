# pimd-tunneling

Path Integral Molecular Dynamics calculations of tunneling splittings
of degenerate rearrangements of molecules, with parallel and serial
implementations, as well as instanton approximation code.

Uses either an Andersen thermostat or a Langevin thermostat with
modified Verlet algorithm to propagate MD trajectories during a
thermodynamic integration. For more details, see Matyus et al, J Chem
Phys 144, 114108 (2016). To extend to multiwell systems, see Matyus
and Althorpe, J Chem Phys 144, 114109 (2016). Also includes Ring
Polymer Instanton calculations; see Richardson and Althorpe, J Chem
Phys 134, 054109 (2011).

Various potential energy surfaces (PES) are implemented; see the
bottom of this for references. To add a new PES, create a new
mcmod_PESNAME.f90 file (the mcmod_malon.f90 is a good example to use),
and modify the makefile accordingly. Please also add the reference for
the PES to this README file.

Uses MKL libraries, as well as some Fortran routines from Numerical
Recipes for Fortran (Press et al). Finally, uses the L-BFGS
implementation by Byrd et al:
[1] R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, ``A limited
memory algorithm for bound constrained optimization'',
SIAM J. Scientific Computing 16 (1995), no. 5, pp. 1190--1208.

[2] C. Zhu, R.H. Byrd, P. Lu, J. Nocedal, ``L-BFGS-B: a
limited memory FORTRAN code for solving bound constrained
optimization problems'', Tech. Report, NAM-11, EECS Department,
Northwestern University, 1994.

----------------------------------
PES References
----------------------------------

1. Malonaldehyde (mcmod_malon.f90)
Mizukami et al, J Chem Phys 141, 114108 (2016)

2. Malonaldehyde (mcmod_wmalon.f90)
Wang et al, J Chem Phys 128, 224314 (2008)

3. Water Dimer (mcmod_waterdimer.f90); MB-pol potential
Babin et al, J Chem Theory Comput 10, 1599 (2014)

4. Clathrate (mcmod_clathrate.f90)
Babin et al, J Chem Theory Comput 10, 1599 (2014)
Akin-Ojo and Szalewicz, J Chem Phyus 123, 134311 (2005)