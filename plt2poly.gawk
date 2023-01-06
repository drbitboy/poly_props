#!/bin/awk -f

### Convert triangular plate shape model
### to polyhedron shape model suitable for volInt

### Input triangular plate shape model format (no blank lines):

###   nvertices nplates
###   X[0]           Y[0]           Z[0]
###   X[1]           Y[1]           Z[1]
###   X[2]           Y[2]           Z[2]
###   ...
###   X[nvertices-1] Y[nvertices-1] Z[nvertices-1]
###   IV[0][0]         IV[0][1]         IV[0][2]
###   IV[1][0]         IV[1][1]         IV[1][2]
###   IV[2][0]         IV[2][1]         IV[2][2]
###   ...
###   IV[nplates-1][0] IV[nplates-1][1] IV[nplates-1][2]


### Output polyhedron shape model format (see volInt.c):

###   nvertices
###   X[0]    Y[0]    Z[0]
###   X[1]    Y[1]    Z[1]
###   X[2]    Y[2]    Z[2]
###   ...
###   X[nv-1] Y[nv-1] Z[nv-1]
###   nplates
###   3 IV[0][0]         IV[0][1]         IV[0][2]
###   3 IV[1][0]         IV[1][1]         IV[1][2]
###   3 IV[2][0]         IV[2][1]         IV[2][2]
###   ...
###   3 IV[nplates-1][0] IV[nplates-1][1] IV[nplates-1][2]

BEGIN { 
  zero=0.0+0.0
  space = " "
}

{
  ### Read up to three elements as floats from triangular plate file
  x=$1+zero
  y=$2+zero
  z=$3+zero

  ### Change values near zero to be zeros
  if ( x < 1e-9 && x > -1e-9) x = 0.0
  if ( y < 1e-9 && y > -1e-9) y = 0.0
  if ( z < 1e-9 && z > -1e-9) z = 0.0
}

### For the first line of the triangular plate file, ...
NR==1 {
  ### ... parse counts of vertices and of plates
  nverts = x
  nplates = y
  ### Print vertex count only; plate count follows vertices
  print nverts
  ### Skip to the next input line (i.e. to first vertex)
  next
}

### While there are vertices remaining to be printed out ...
nverts > 0 {
  ### ... print each vertex [X Y Z]
  print x space y space z
  ### Decrement remaining vertex count;
  ### when last vertex is printed, print the triangular plate count
  if (--nverts == 0) { print nplates }
  ### Skip to the next input line
  next
}

### After all vertices are read & printed, ...
nverts == 0 {
  ### ... print each triangular plate's 0-based vertices' indices,
  ### with a numeric prefix inidicating the number of vertices i.e. 3
  print "3 " $0
}
