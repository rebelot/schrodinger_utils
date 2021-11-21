#!/bin/bash


USAGE="mkhosts.sh [-h|--help] username version destination
example: ./mkhosts.sh marty.mcfly.stud 2015-3 ~/.schrodinger/schrodinger.hosts_2015-3"

[[ $1 =~ (-h|--help) ]] && { echo "$USAGE"; exit; }
[[ $# -ne 3 ]] && { echo "$USAGE"; exit; }

DIR="$(dirname "$0")/../templates/"

sed -e "s/USERNAME/$1/" -e "s/VERSION/$2/" "$DIR"/schrodinger.hosts.template > "$3"
