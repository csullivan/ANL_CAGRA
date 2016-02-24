#!/bin/bash

#source ~/.cshrc_oldroot
MERGDIR="/mnt/analysis/rcnp-e441/cagragr"
ROOTDIR="/mnt/analysis/rcnp-e441/analysis/rootfiles"
CHATDIR="$PWD/bin"
cd $CHATDIR
 
if [ $# -ne 2 ] 
  then
   echo "Specify only the run number to sort"
  exit 1
fi

RUN=$1
RFILE=$2

echo "GEBSort started sorting run $RUN at `date`"
echo "$CHATDIR/GEBSort -input disk $MERGDIR/GEBMerged_run$RUN.gtd_000 -rootfile $ROOTDIR/$RFILE.root RECREATE -chat $CHATDIR/GEBSort.chat"
$CHATDIR/GEBSort -input disk $MERGDIR/GEBMerged_run$RUN.gtd_000 -rootfile $ROOTDIR/$RFILE.root RECREATE -chat $CHATDIR/GEBSort.chat
echo "$CHATDIR/GEBSort -input disk $MERGDIR/GEBMerged_run$RUN.gtd_000 -rootfile $ROOTDIR/$RFILE.root RECREATE -chat $CHATDIR/GEBSort.chat"
echo "GEBSort DONE at `date`"

#exit
