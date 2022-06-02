filename=lee


o:
	gfortran ${filename}.f90 -o ${filename}


clean:
	rm -f ${filename}
