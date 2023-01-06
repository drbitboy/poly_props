#!/bin/bash

me="$0"
medir="`dirname \"$me\"`"

one="$1"

[ -r "$one" ] || ( echo "Usage:  $0 file.shape; exiting " ; false ) || exit 1

spud2poly.gawk "$one" | volInt -
