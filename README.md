# SuperChic2

By default the program makes use of the LHAPDF library, so please add the following lines to the shell (bash) login script:

export LHAPDFSYS=/yourpath/LHAPDF-X.Y.Z

export PATH=${PATH}:${LHAPDFSYS}/bin

LD LIBRARY PATH=${LD LIBRARY PATH}:${LHAPDFSYS}/lib

To compile and to run:

$ make

$ cd bin/

$ ./init < input.DAT

$ ./superchic < input.DAT

