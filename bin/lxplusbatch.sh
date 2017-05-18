# Lxplus Batch Job Script
set INPUT_FILE=“input_gluino_decay_cc.DAT”

source /afs/cern.ch/sw/lcg/external/gcc/4.9.3/x86_64-slc6/setup.sh
source /afs/cern.ch/sw/lcg/app/releases/ROOT/6.06.00/x86_64-slc6-gcc49-opt/root/bin/thisroot.sh
export LHAPATH=/afs/cern.ch/user/s/ssilvado/local/share/LHAPDF/
cd /afs/cern.ch/user/s/ssilvado/superchicv2.04
./afs/cern.ch/user/s/ssilvado/superchicv2.04/bin/init < $INPUT_FILE
./afs/cern.ch/user/s/ssilvado/superchicv2.04/bin/superchic < $INPUT_FILE
