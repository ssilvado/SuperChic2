# SuperChic2

By default the program makes use of the LHAPDF library, to build LHAPDF:

  	wget http://www.hepforge.org/archive/lhapdf/LHAPDF-6.X.Y.tar.gz -O- | tar xz
  	cd LHAPDF-6.X.Y
  	./configure --prefix=$PWD/../local
  	make -j2 && make install
  	cd ..
	
After, add the following lines to the shell (bash) login script:

  	export LHAPDFSYS=/yourpath/LHAPDF-X.Y.Z
	export PATH=${PATH}:${LHAPDFSYS}/bin
  	LD LIBRARY PATH=${LD LIBRARY PATH}:${LHAPDFSYS}/lib

To compile and to run SuperChic2:

  	make
  	cd bin/
  	./init < input.DAT
  	./superchic < input.DAT
