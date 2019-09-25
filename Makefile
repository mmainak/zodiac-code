#
# This is the makefile for zodiac.
# To compile the code, just type make.  Such an approach makes
# recompilation of the code easy, recompiling only as necessary
# to account for recent changes to the code.
#
COMPILER = mpif90 -O3 #-mcmodel=medium -openmp
COMPOPTS = -I/usr/local/include -I/opt/intel/fce/10.1.008/lib/ #-i-dynamic 



#LINKOPTS = -L/usr/local/lib  -lrfftw -lfftw
# Location where the optional netcdf include file (netcdf.inc) is installed
INCLUDEDIR = /usr/local/include



PARALLEL = TRUE
LES      = FALSE

# Option to compile with the NetCDF libraries
NETCDF = TRUE




ifeq ($(LES),TRUE)
LES_CHAN = les_chan_mpi.o les_chan_th_mpi.o
else
LES_CHAN =  no_les.o 
endif


# Use the parameters to set flags
 ifeq ($(NETCDF),TRUE)
 LINKOPTS = -L/usr/local/lib  -lrfftw -lfftw -lnetcdf -lnetcdff
 NETCDF_o = netcdf.o
 else
 LINKOPTS = -L/usr/local/lib  -lrfftw -lfftw 
 NETCDF_o = no_netcdf.o
 endif



ifeq ($(PARALLEL),TRUE)
zodiac: zodiac.f periodic.o channel.o $(LES_CHAN) $(NETCDF_o)\
	duct.o cavity.o fft.o boundary.o flow_statistic.o mpi_chan.o wall_model.o\
	header grid_def 
	$(COMPILER) $(COMPOPTS) zodiac.f -o zodiac \
	periodic.o channel.o $(LES_CHAN) $(NETCDF_o) \
	duct.o cavity.o fft.o boundary.o flow_statistic.o mpi_chan.o wall_model.o $(LINKOPTS)
else
zodiac: zodiac.f periodic.o channe.o $(LES_CHAN) $(NETCDF_o) \
	duct.o cavity.o fft.o mpi_chan_serial.o wall_model.o\
	header grid_def 
	$(COMPILER) $(COMPOPTS) zodiac.f -o zodiac \
        periodic.o channel.o $(LES_CHAN) \
	duct.o cavity.o fft.o mpi_chan_serial.o wall_model.o $(LINKOPTS)
endif

periodic.o: periodic.f fft.o header grid_def
	$(COMPILER) $(COMPOPTS) -c periodic.f

ifeq ($(PARALLEL),TRUE)
channel.o: channel.f fft.o mpi_chan.o boundary.o flow_statistic.o wall_model.o header grid_def
	$(COMPILER) $(COMPOPTS) -c channel.f
else
channel.o: channel.f fft.o mpi_chan_serial.o boundary.o flow_statistic.o wall_model.o header grid_def
	$(COMPILER) $(COMPOPTS) -c channel.f
endif

ifeq ($(LES),TRUE) 
les_chan_mpi.o: les_chan_mpi.f fft.o header header_les grid_def
	$(COMPILER) $(COMPOPTS) -c les_chan_mpi.f

les_chan_th_mpi.o: les_chan_th_mpi.f fft.o header header_les grid_def
	$(COMPILER) $(COMPOPTS) -c les_chan_th_mpi.f
else
no_les.o: no_les.f
	$(COMPILER) $(COMPOPTS) -c no_les.f
endif

ifeq ($(NETCDF),TRUE)
netcdf.o: netcdf.f header grid_def
	$(COMPILER) $(COMPOPTS) -c netcdf.f
else
no_netcdf.o: no_netcdf.f 
	$(COMPILER) $(COMPOPTS) -c no_netcdf.f
endif

ifeq ($(PARALLEL),TRUE)
mpi_chan.o: mpi_chan.f header header_mpi grid_def
	$(COMPILER) $(COMPOPTS) -c mpi_chan.f
else
mpi_chan_serial.o: mpi_chan_serial.f header header_mpi grid_def
	$(COMPILER) $(COMPOPTS) -c mpi_chan_serial.f
endif

duct.o: duct.f header grid_def
	$(COMPILER) $(COMPOPTS) -c duct.f

cavity.o: cavity.f header grid_def
	$(COMPILER) $(COMPOPTS) -c cavity.f

fft.o:  fft.f header grid_def
	$(COMPILER) $(COMPOPTS) -c fft.f

wall_model.o:	wall_model.f header grid_def
	$(COMPILER) $(COMPOPTS) -c wall_model.f

flow_statistic.o: flow_statistic.f header grid_def
	$(COMPILER) $(COMPOPTS) -c flow_statistic.f

boundary.o: boundary.f header grid_def
	$(COMPILER) $(COMPOPTS) -c boundary.f

clean:
	rm -f *.o *.o* output.txt grids/sponge_* fort.* *~ zodiac core

# Compiler specific notes:
#
# Compilation with Absoft Linux Fortran 77 appears to be impossible, as it
# cannot handle the INTEGER*8 option required by FFTW.  If someone finds
# a way around this, please let me know.
# 
# Compilation with Absoft Linux Fortran 90 is possible, but the option
# -YEXT_NAMES=LCS must be used as one of the link options so the compiler
# can find the lowercase external library function names.
#
# Compilation with Lahey Fortran 95 (lf95) is possible, but there is an
# underscore incompatability with the FFTW libraries, which are compiled
# with g77.  To get around this, you need to go into fft.f and add 
# trailing underscores to the name of every fftw function where they
# appear throughout the code.

