In addition to the original grid files used as input, this directory contains:

- Data files containing just the centers, to facilitate comparison with the test output; these were generated with commands like `ncks -v x,y -d nxp,1,,2 -d nyp,1,,2 C96_grid.tile6.nc C96_grid.tile6.centers.nc` (which extracts the centers from the arrays that interleave corners and centers).
