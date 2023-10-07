#!/usr/bin/bash


# 
# Setup arguments
# 

jobId=${1:-FULL} # default FAST
nEvents=${2:-10}
sim=${3:-true}
pid=${4:-5} # default muon


# 
# Print job parameters
# 
echo "Usage:\n rbfull <jobId> <nEvents> <run_sim> <pid>"
echo "jobId=${jobId}, nEvents=${nEvents}, run_sim=${sim}, pid=${pid}"
# source env.sh

# 
# Run simulation if needed
# 
if [ "$sim" = true ] ; then
    strongrandom=`od -vAn -N3 -tu4 < /dev/urandom | tr -d '[:space:]'`
    echo "strong random ${strongrandom}"
    echo root4star -b -q -l 'sim/gen.C( '"${nEvents}"','"${strongrandom}"')'
    time root4star -b -q -l 'sim/gen.C( '"${nEvents}"','"${strongrandom}"')'
else
    echo "Skipping simulation - use existing sim.fzd"
fi


# 
# Run tracking
# 
echo root4star -b -q -l 'sim/fwd_tracking.C( '"${nEvents}"', "sim.fzd", "sim/fst_track.xml" )' 
time root4star -b -q -l 'sim/fwd_tracking.C( '"${nEvents}"', "sim.fzd", "sim/fst_track.xml" )'


# 
# Copy output files into JOB prefixed names
# 
if [ -f "fwdtree.root" ]; then
    mv fwdtree.root ${jobId}_fwdtree.root
else 
    echo "WARNING: no fwdtree.root produced"
fi 

if [ -f "fst_track.root" ]; then
    mv fst_track.root ${jobId}_trackingQA.root
else 
    echo "WARNING: no fst_track.root produced (tracking QA)"
fi 

if [ -f "fcstrk.root" ]; then
    mv fcstrk.root ${jobId}_fcstrk.root
else 
    echo "No fcstrk.root output found, double check if FCS is included"
fi


if [ -f "FwdAna.root" ]; then
    mv FwdAna.root ${jobId}_FwdAna.root
else 
    echo "No FwdAna.root output found, double check if FCS is included"
fi


