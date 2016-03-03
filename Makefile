# ARPACK++ v1.2 2/18/2000
# c++ interface to ARPACK code.
# examples/product/simple directory makefile.

# including other makefiles.
CC=g++
HERE = /home/sgeraedt/basis_change
MYDIR = /home/sgeraedt/myClibrary/
CFLAGS=-Wall -I$(HERE) -I$(MYDIR)
LIBS= $(MYDIR)utils.o  -lgfortran 

a.out: weir3.o basis.o
	$(CC) -O3 $(CFLAGS) -o a.out weir3.o basis.o $(LIBS)
	
clean:
	rm -f *~ *.o a.out

%.o:	%.cpp
	$(CC) -O3 $(CFLAGS) -c $<

%.o:	%.f90
	gfortran -O3 $(CFLAGS) -c $<
