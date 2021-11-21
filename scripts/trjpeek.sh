#!/bin/bash

if [[ $# -lt 4 ]]; then
    echo trjpeek.sh uname node jobdirname destpath
    exit
fi

uname="$1"
node="$2"
jobdirname="$3"
destpath="$4"

rsync --progress -th -r $uname@xlence.disfeb.unimi.it:/scratches/$node/$uname/$jobdirname/\*{-in.cms,_trj} $destpath
