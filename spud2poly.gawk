#!/bin/awk -f

### Convert SPUD shape model
### to polyhedron shape model suitable for volInt

### Input SPUD shape model format (no blank lines):

###   -90.0           0.0           RADIUS
###   -90.0           deltalon      RADIUS
###   ...
###   -90.0           360.0         RADIUS
###   -90.0+deltalat  0.0           RADIUS
###   ...
###   -90.0+deltalat  deltalon      RADIUS
###   -90.0+deltalat  360.0         RADIUS
###   ...
###   90.0            360.0         RADIUS

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

### Initialization
BEGIN { 
  ### vertex and plate numbers
  ivert=0-1
  iplat=0-1

  ### Dictionary of previous parallel's vertex numbers; key is longitude
  split("", iverts)

  ### String of all vertices' Cartesian coordinates, with newlines
  vertstring="\n"
  ### String of all plates' zero-based vertices, with newlines
  platstring=""
  ### Constants
  rpd = atan2(1.0,1.0) / 45.0
  zero=0.0+0.0
  space = " "
  newline = "\n"
}

### For each [Lat Lon Radius] line ...
{
  ### ... increment vertex number, parse vertex as floats from SPUD file
  ++ivert
  lat=$1+zero
  wlon=$2+zero
  rad=$3+zero

  ### Convert spherical coordinates to Cartesian coordinates
  ### N.B. Longitude is West
  z = rad * sin(lat * rpd)
  rcoslat = rad * cos(lat * rpd)
  x = cos(wlon * rpd) * rcoslat
  y = 0.0 - (sin(wlon * rpd) * rcoslat)

  ### Change values near zero to be zeros
  if ( x < 1e-9 && x > -1e-9) x = 0.0
  if ( y < 1e-9 && y > -1e-9) y = 0.0
  if ( z < 1e-9 && z > -1e-9) z = 0.0

  ### Append vertex to vertex string
  vertstring = vertstring x space y space z newline
}

### For the first line of the SPUD file, ...
NR==1 {
  ### ... initialize previous/last values
  prevlat=lat
  lastlat=lat
  lastwlon=wlon
  lastivert=ivert
}

### When new parallel starts, ...
lat > lastlat {
  ### ... update previous parallel's lat, so 1s longitude does nothing
  prevlat=lastlat
}

### Once there are two latitudes on the current parallel ...
prevlat < lastlat && lastlat==lat {
  ### Increment the plate number
  ++iplat
  ### Print triangles at the poles, quadrilateral elswhere
  ### - get previous parallel's vertices from iverts dictionary
  if ( prevlat == -90.0) {
    addstr = 3 space lastivert space ivert        space iverts[wlon]
  } else if ( lat == 90.0) {
    addstr = 3 space ivert     space iverts[wlon] space iverts[lastwlon]
  } else {
    addstr = 4 space lastivert space ivert        space iverts[wlon]      space iverts[lastwlon]
  }
  ### Extend plate string
  platstring = platstring addstr newline
}

{
  ### Update vertex number at last lon for next parallel
  iverts[lastwlon]=lastivert
  ### Keep track last of lat, lon, vertex number
  lastlat=lat
  lastwlon=wlon
  lastivert=ivert
}

### Final action after all lines have been read ...
END { 
  ### ... print vertex count, vertices, plate count, plates
  print ivert+1
  print vertstring
  print newline iplat+1
  print platstring
}
