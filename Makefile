# if there is no Makefile.in then use the template
# --------------------------------------------------
ifeq ($(shell test -e Makefile.in && echo 1), 1)
MAKEFILE_IN = $(PWD)/Makefile.in
else
MAKEFILE_IN = $(PWD)/Makefile.in.template
endif
include $(MAKEFILE_IN)



ifeq ($(HAVE_HDF5), 1)
HDF5_I = -I$(HDF5_HOME)/include
HDF5_L = -L$(HDF5_HOME)/lib -lhdf5
endif

ifeq ($(HAVE_FFTW), 1)
FFTW_I = -I$(FFTW_HOME)/include
FFTW_L = -L$(FFTW_HOME)/lib -lfftw3
endif

ifeq ($(HAVE_RNPL), 1)
RNPL_I = -I$(RNPL_HOME)/include
RNPL_L = -L$(RNPL_HOME)/lib -lrnpl
endif


COW_LIB = cow/libcow.a

SRC = $(wildcard *.c)
OBJ = $(SRC:.c=.o)


default : ffe


%.o : %.c ffe.h
	$(CC) $(CFLAGS) $< -c $(HDF5_I) $(FFTW_I) $(RNPL_I)

ffe : $(OBJ) $(COW_LIB)
	$(CC) $(CFLAGS) $(CLIBS) $^ $(HDF5_L) $(FFTW_L) $(RNPL_L) -o $@

$(COW_LIB) : .FORCE
	$(MAKE) -C cow MAKEFILE_IN=$(MAKEFILE_IN)

clean :
	$(MAKE) -C cow clean MAKEFILE_IN=$(MAKEFILE_IN)
	$(RM) ffe $(OBJ)

.FORCE :
