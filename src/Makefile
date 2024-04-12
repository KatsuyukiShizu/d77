FORT = gfortran
FFLAGS = -Wall -fbounds-check

#FORT = ifort
#FFLAGS = -g -heap-arrays -parallel -O2 -mkl -Warn all

TARGET = d77.exe
OBJECTS = ifmod.o \
global_constants.o \
global_read_input.o \
func_int_pgf.o \
global_read_data.o \
main.o \
write_messages.o \
calc_int_cgf.o \
dipole_int_pgf.o \
vc_int_pgf.o \
soc_int_pgf.o \
func_int_slater.o \
calc_def.o \
calc_rlmn_0.o \
calc_opr_dmat.o \
calc_dmat_cgf.o \
calc_dmat_soc_cgf.o \
num2char.o \
write_g16_log_format.o \
write_mat.o \
 
.SUFFIXES : .o .f90
.f90.o:
	${FORT} ${LDFLAGS} ${FFLAGS} -c $<
 
${TARGET} : ${OBJECTS}
	${FORT} -o $@ ${OBJECTS} ${LDFLAGS} ${FFLAGS}
