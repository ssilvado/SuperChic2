# SuperChic2
# By default the program makes use of the LHAPDF library

please add the following lines to the shell (bash) login script:
export LHAPDFSYS=/yourpath/LHAPDF-X.Y.Z
export PATH=${PATH}:${LHAPDFSYS}/bin
LD LIBRARY PATH=${LD LIBRARY PATH}:${LHAPDFSYS}/lib

$ make
$ cd bin/
$ ./init < input.DAT
$ ./superchic < input.DAT
