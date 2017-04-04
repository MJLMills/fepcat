FC = gfortran
FCFLAGS = -ffree-line-length-none -fimplicit-none -fmax-errors=0 -std=f95 -Wall -Werror -pedantic-errors

PROGRAMS = Fepcat

all: $(PROGRAMS)

StatisticalFunctions.o: StatisticalFunctions.f90
Output.o:               Output.f90
FileIO.o:               FileIO.f90
Matrix.o:               Matrix.f90
ArrayUtil.o:            ArrayUtil.f90
Log.o:                  FileIO.f90
DCDHeader.o:            DCDHeader.f90
InternalCoords.o:       InternalCoords.f90
DCDFiles.o:             DCDFiles.f90   DCDHeader.o
Input.o:                Input.f90      FileIO.o Log.o DCDFiles.o
Data.o:                 Data.f90       Input.o Matrix.o StatisticalFunctions.o InternalCoords.o
Util.o:                 Util.f90       Log.f90 Input.f90 Data.f90
Analysis.o:             Analysis.f90   Output.o StatisticalFunctions.o Input.o Data.o
FreeEnergy.o:           FreeEnergy.f90 Data.o Input.o StatisticalFunctions.o ArrayUtil.o
Qfep.o:                 Qfep.f90       Data.o StatisticalFunctions.o FileIO.o FreeEnergy.o
Fepcat.o:               Fepcat.f90     Data.o Input.o Log.o Qfep.o Analysis.o Util.o

Fepcat: Fepcat.o Input.o Output.o Data.o FileIO.o Analysis.o StatisticalFunctions.o FreeEnergy.o Qfep.o Log.o ArrayUtil.o Matrix.o DCDHeader.o DCDFiles.o InternalCoords.o Util.o

%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.f90*
	$(FC) $(FCFLAGS) -c $<

clean:
	rm -f *.o *.mod *.MOD

install:
	mv Fep ../bin
