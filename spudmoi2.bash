#!/bin/bash

me="$0"
medir="`dirname \"$me\"`"

one="$1"

[ -r "$one" ] || ( echo "Usage:  $0 file.shape; exiting " ; false ) || exit 1

spud2poly.gawk "$one" | volInt - \
| gawk '
BEGIN {
  split("",lats)
  split("",lons)
  split("",radii)
  split("",eigvecs)
  spudCount=0
  dLat=-1
  dLon=-1
  lastLat=-90
  lastLon=0.0
  stderr="/dev/stderr"
  eigvecN=-1
  cm=","
}
{
  inShape=0
}
FILENAME!="-" {
  inShape=1
  lat=0.0+$1
  lon=0.0+$2
  radius=0.0+$3
  dt=lat-lastLat
  dg=(360.0+lon-lastLon)%360
}
function addLatLon() {
  lats[spudCount]=lat
  lons[spudCount]=lon
  radii[spudCount++]=radius
  lastLat=lat
  lastLon=lon
  if ( lat==-90 && lon==360) nLon=spudCount
  if ( lat==90 && lon==360) nLat=spudCount/nLon
}
inShape==1 && FNR==1 {
  addLatLon()
  next
}
inShape==1 && dg>0 && dt==0 {
  if (dLon<0) dLon=dg
  if ( dg!=dLon ) {
    print "### Longitude error in shape file " FILENAME " line " FNR >> stderr
    print "### " $0 >> stderr
    exit 1
  }
  addLatLon()
  next
}
inShape==1 && dt>0 && dg==0 {
  if (dLat<0) dLat=dt
  if ( dt!=dLat ) {
    print "Latitude error in shape file " FILENAME " line " FNR >> stderr
    print "### " $0 >> stderr
    exit 2
  }
  addLatLon()
  next
}
inShape==1 {
  print "### Error in shape file " FILENAME " line " FNR >> stderr
  print "### " $0 >> stderr
  print "### " lastLat, lat, dt >> stderr
  print "### " lastLon, lon, dg >> stderr
  exit 3
}
########################################################################
### Parse the volInt output

{
  print "### " $0 >> stderr
}

$0~/^Eigenvectors [(]in columns[)]: *$/ && eigvecN==-1 {
  eigvecN=0+0
  ###print "### " $0
  next
}
function evIdx( f, s) {
  return f cm s
}
eigvecN>(0-1) && eigvecN<3 {
  ###print "### -> " eigvecN, $0
  eigvecs[evIdx(eigvecN,0)]=0.0+$1
  eigvecs[evIdx(eigvecN,1)]=0.0+$2
  eigvecs[evIdx(eigvecN,2)]=0.0+$3
  ++eigvecN
  next
}
function rInterp( latArg, lonArg) {

  lclLat=(latArg+90)
  if ( lclLat > 180) lclLat=180
  if ( lclLat < 0) lclLat=0
  fLat=lclLat%dLat
  iLat0=(lclLat-fLat)/dLat
  iLat1=iLat0+1
  if ( iLat1 >= nLat) iLat1=nLat-1

  lclLon=lonArg
  if ( lclLon > 360) lclLon=360
  if ( lclLon < 0) lclLon=0
  fLon=lclLon%dLon
  iLon0=(lclLon-fLon)/dLon
  if ( iLon0 >= nLon) iLon0=nLon
  iLon1=iLon0+1
  if ( iLon1 >= nLon) iLon1=nLon-1

  ll=radii[iLon0+(iLat0*nLon)]
  lr=radii[iLon1+(iLat0*nLon)]
  ul=radii[iLon0+(iLat1*nLon)]
  ur=radii[iLon1+(iLat1*nLon)]

  lrad = ((fLon*lr) + ((dLon-fLon)*ll)) / dLon
  urad = ((fLon*ur) + ((dLon-fLon)*ul)) / dLon

  ###print "### ll,lr,ul,ur = " ll,lr,ul,ur
  ###print "### latArg,lonArg,lclLat,lclLon = " latArg,lonArg,lclLat,lclLon
  ###print "### iLat0,iLat1,iLon0,iLon1,fLat,fLon = " iLat0,iLat1,iLon0,iLon1, fLat, fLon
  return ((fLat*urad) + ((dLat-fLat)*lrad)) / dLat
}
END {
  ###print "### " dLat, dLon, spudCount, nLat, nLon
  if (1) for ( i=0; i<3; ++i ) {
    for ( j=0; j<3; ++j ) {
       idx=evIdx(i,j)
       ###printf " [%s]%s", idx, eigvecs[idx]
    }
    ###print ""
  }

  rpd=atan2(1.0,1.0)/45.0
  for ( iRad=0; iRad<spudCount; ++iRad) {
    lat=lats[iRad]
    latr=lat*rpd
    lon=lons[iRad]
    lonr=lon*rpd
    clat=cos(latr)
    slat=sin(latr)
    clon=cos(lonr)
    slon=sin(lonr)
    split("",xyz)
    xyz[0]=clon*clat
    xyz[1]=slon*clat
    xyz[2]=slat

    split("",oldXyz)
    for ( i=0; i<3; ++i ) { oldXyz[1]=0.+0.0 }

    for ( i=0; i<3; ++i ) {
      for ( j=0; j<3; ++j ) {
        oldXyz[i] += (xyz[j] * eigvecs[evIdx(i,j)])
      }
    }
    ###print "### x,y,z = " xyz[0],xyz[1],xyz[2]
    ###print "### old x,y,z = " oldXyz[0],oldXyz[1],oldXyz[2]

    oldLon = (360.0 + atan2(oldXyz[1],oldXyz[0])/rpd) % 360

    oldLat = atan2(oldXyz[2],sqrt(oldXyz[0]*oldXyz[0]+oldXyz[1]*oldXyz[1]))/rpd

    ###print "### lat,lon,oldLat,oldLon = " lat,lon,oldLat,oldLon

    printf "%10.4lf%10.4lf%10.4lf\n", lat, lon, rInterp( oldLat, oldLon)
  }
}
' - "$one"

exit

