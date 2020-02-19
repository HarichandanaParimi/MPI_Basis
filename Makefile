

FILESTOINCLUDE=approx.cpp assignment-mpi-basic.pdf params.sh \
		hello_world/hello.cpp hello_world/Makefile hello_world/run.sh \
		matvec/Makefile matvec/mpi_matmul.cpp matvec/plot.sh matvec/queue.sh matvec/run_matmul.sh matvec/test.sh \
		num_int/cases.txt  num_int/libfunctions.a  num_int/Makefile  num_int/mpi_num_int.cpp  num_int/plot.sh  num_int/queue.sh  num_int/run_numint.sh  num_int/test.sh



assignment-mpi-basic.pdf: assignment-mpi-basic.tex
	pdflatex assignment-mpi-basic.tex

assignment-mpi-basic.tgz: $(FILESTOINCLUDE)
	tar zcvf $@ $(FILESTOINCLUDE)
