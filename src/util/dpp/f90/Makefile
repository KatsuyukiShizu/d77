FORT = gfortran
FFLAGS = -Wall -fbounds-check
 
TARGET = dpp.exe
OBJECTS = global_constants.o \
read_txt.o fchk2inputs.o \
gen_cicoef.o \
gen_g16_cicoef_cisd.o \
gen_g16_cicoef_td.o \
gen_g16_cicoef_td_unres.o \
main.o
 
.SUFFIXES : .o .f90
.f90.o:
	${FORT} ${LDFLAGS} ${FFLAGS} -c $<
 
${TARGET} : ${OBJECTS}
	${FORT} -o $@ ${OBJECTS} ${LDFLAGS} ${FFLAGS}
