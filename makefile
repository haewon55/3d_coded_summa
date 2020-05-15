EXECS= 3d_coded_summa replication_3d_summa sys_3d_coded_summa
MPICC?=mpicc
CC = gcc 

IDIR = $MKLROOT/include
LDIR = $MKLROOT/lib/intel64_lin 
CFLAGS=-I$(IDIR) -L$(LDIR)  

LIBS = -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm

all: ${EXECS}

$(EXECS): %: %.c
	${MPICC} -o $@ $< $(CFLAGS) $(LIBS) 

clean:
	rm ${EXECS}
