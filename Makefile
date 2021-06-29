F90	   = gfortran
F90LINKER = $(F90)
OPTFLAGS  = -O3 -llapack -lblas
FFLAGS    = $(OPTFLAGS)
DATE      = $(shell date '+%y%m%d')

TARGET = SKS-1Dx
default: $(TARGET)
OBJS = \
  output_runout.o \
  hll_flux_H.o \
  hllc_flux_H.o \
  mnewt_FC.o \
  cons_to_prim_FC.o \
  godunov_flux_FC.o \
  hll_flux_L.o \
  hllc_flux_L.o \
  cons_to_prim.o \
  fractional_step.o \
  cfl_condition.o \
  output_points.o \
  output_mass.o \
  output_front.o \
  output.o \
  set_initial.o \
  main.o
CNST = set_parameter.o variable.o

$(TARGET):  $(CNST) $(OBJS) 
	$(F90LINKER) $(OPTFLAGS) -o $(TARGET) $(CNST) $(OBJS) ;\
	/bin/rm -f *~
clean:
	/bin/rm -f *.o $(TARGET) *~ *.mod $(DATE).f90
.f90.o:
	$(F90) $(FFLAGS) -c  $< ;\
	/bin/cat $< >> $(DATE).f90

# resetting default suffixes rule
.SUFFIXES:
.SUFFIXES: .f90 .o

#depencency
main.o: set_initial.o output.o output_front.o output_mass.o output_points.o cfl_condition.o fractional_step.o output_runout.o
fractional_step.o: cons_to_prim.o hll_flux_L.o hllc_flux_L.o godunov_flux_FC.o cons_to_prim_FC.o hll_flux_H.o hllc_flux_H.o
cons_to_prim_FC.o: mnewt_FC.o
$(OBJS): $(CNST)

