#!/bin/bash

params=""
while (( "$#" )); do
  case "$1" in
    -h|--help)
      echo 'usage: trj_eac cms trj [-e|--extract ASL args] [-a|--align ASL args] [-c|--center ASL args] [-o|--output filename]'
      exit 1
      ;;
    -e|--extract)
      easl=$2
      [[ ! $3 =~ ^- ]] && { earg=$3; shift; }
      shift 2
      ;;
    -a|--align)
      aasl=$2
      [[ ! $3 =~ ^- ]] && { aarg=$3; shift; }
      shift 2
      ;;
    -c|--center)
      casl=$2
      [[ ! $3 =~ ^- ]] && { carg=$3; shift; }
      shift 2
      ;;
    -o|--output)
      out="$2-"
      shift 2
      ;;
    --) # end argument parsing
      shift
      break
      ;;
    -*|--*=) # unsupported flags
      echo "Error: Unsupported flag $1" >&2
      exit 1
      ;;
    *) # preserve positional arguments
      params="$params $1"
      shift
      ;;
  esac
done

params=($params)
cms="${params[0]}"
trj="${params[1]}"

if [[ ! -z $easl ]]; then
    echo "extracting \"$easl\" subsystem from $cms..."
    [[ -z $step ]] && step=e || step="$step""e"
    $SCHRODINGER/run trj_extract_subsystem.py "$cms" "$out$step" -t "$trj" -asl "$easl" $earg
    cms="$out$step""-out.cms"
    trj="$out$step""_trj"
    echo "written $cms, $trj"
fi

if [[ ! -z $aasl ]]; then
    echo "aligning \"$aasl\" subsystem from $cms..."
    [[ -z $step ]] && step=a || step="$step""a"
    $SCHRODINGER/run trj_align.py  "$cms" "$trj" "$out$step" -asl "$aasl" $aarg
    cms="$out$step""-out.cms"
    trj="$out$step""_trj"
    echo "written $cms, $trj"
fi

if [[ ! -z $casl ]]; then
    echo "centering \"$casl\" subsystem from $cms..."
    [[ -z $step ]] && step=c || step="$step""c"
    $SCHRODINGER/run trj_center.py "$cms" "$out$step" -t $trj -asl "$casl" $carg
    cms="$out$step""-out.cms"
    trj="$out$step""_trj"
    echo "written $cms, $trj"
fi
