## Uncomment the appropriate fortran 90 compiler.
F90=pgf90
#F90=gfortran
#F90=ifort

molanal.new : molanal.new.c
	pgcc -O3 molanal.new.c -lm -o molanal.new

v2 : molanal2.c
	pgcc -O3 molanal2.c -lm -o v2.x

metal : molanal.metal.c
	pgcc -g molanal.metal.c -lm -o metal.x

h2o : molanal_h2o.c
	pgcc -O3 molanal_h2o.c -lm -o h2o.x

gr_calc : gr_calc.f90
	$(F90) gr_calc.f90 -o gr_calc

clean:
	rm -f *.o gr_calc molanal.new
