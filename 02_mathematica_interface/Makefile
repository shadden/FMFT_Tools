# This makefile can be used to build all or some of the sample
# programs.  To build all of them, use the command
# 'make all'.  To build one, say addtwo, use the command
# 'make addtwo'.

# Portions of this makefile require the use of GNU make.
# see http://www.gnu.org/software/make for more information.

VERSION=10.0
MLINKDIR = /Applications/Mathematica.app/SystemFiles/Links/MathLink/DeveloperKit/MacOSX-x86-64
SYS = MacOSX-x86-64
CADDSDIR = ${MLINKDIR}/CompilerAdditions
EXTRA_CFLAGS=

INCDIR = ${CADDSDIR}
LIBDIR = ${CADDSDIR}/AlternativeLibraries

MPREP = ${CADDSDIR}/mprep
RM = rm

CC = /usr/bin/clang
CXX = /usr/bin/clang++

BINARIES = ML_fmft

all : $(BINARIES)

ML_fmft : ML_fmfttm.o ML_fmft.o fmft.o nrutil.o
	${CXX} ${EXTRA_CFLAGS} -I${INCDIR} $^  -L${LIBDIR} -lMLi4 -lstdc++ -framework Foundation -o $@
fmft.o: fmft.c nrutil.o
	${CC} -c ${EXTRA_CFLAGS} $^
nrutil.o: nrutil.c
	${CC} -c ${EXTRA_CFLAGS} $^
.c.o :
	${CC} -c ${EXTRA_CFLAGS} -I${INCDIR} $<

addtwotm.c : addtwo.tm
	${MPREP} $? -o $@
ML_fmfttm.c : ML_fmft.tm
	${MPREP} $? -o $@
clean :
	@ ${RM} -rf *.o *tm.c $(BINARIES)
