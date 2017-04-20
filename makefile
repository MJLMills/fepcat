FC = gfortran
FCFLAGS = -ffree-line-length-none -fimplicit-none -fmax-errors=0 -std=f95 -Wall -Werror -pedantic-errors

PROGRAMS = Fepcat Qfep FepMovie Fep2D

all: $(PROGRAMS)

StatisticalFunctions.o: StatisticalFunctions.f90
Output.o:               Output.f90
FileIO.o:               FileIO.f90
Matrix.o:               Matrix.f90
ArrayUtil.o:            ArrayUtil.f90
Log.o:                  FileIO.f90
DCDHeader.o:            DCDHeader.f90
InternalCoords.o:       InternalCoords.f90
Input2D.o:              Input2D.f90 FileIO.o 
Movies.o:               Movies.f90          FileIO.o Output.o
DCDFiles.o:             DCDFiles.f90        DCDHeader.o
Input.o:                Input.f90           FileIO.o Log.o DCDFiles.o
Data.o:                 Data.f90            Input.o Input2D.o Matrix.o StatisticalFunctions.o InternalCoords.o
FepMovie.o:             FepMovie.f90        Util.o Data.o Movies.o Log.o
Util.o:                 Util.f90            Log.f90 Input.o Data.o
Analysis.o:             Analysis.f90        Output.o StatisticalFunctions.o Input.o Data.o FreeEnergy.o
Fep2D.o:                Fep2D.f90           Analysis.o Util.o Input2D.o Input.o
FreeEnergy.o:           FreeEnergy.f90      Data.o Input.o StatisticalFunctions.o ArrayUtil.o
Qfep.o:                 Qfep.f90            Data.o StatisticalFunctions.o FileIO.o FreeEnergy.o Util.o
Fepcat.o:               Fepcat.f90          Data.o Input.o Log.o Analysis.o Util.o

Fepcat:   Fepcat.o Input.o Output.o Data.o FileIO.o Analysis.o StatisticalFunctions.o FreeEnergy.o Log.o ArrayUtil.o Matrix.o DCDHeader.o DCDFiles.o InternalCoords.o Util.o Movies.o Input2d.o
Qfep:     Qfep.o Data.o Input.o Matrix.o StatisticalFunctions.o InternalCoords.o Fileio.o FreeEnergy.o ArrayUtil.o Log.o DCDFiles.o DCDHeader.o Util.o Input2d.o
FepMovie: FepMovie.o Util.o Data.o Movies.o Log.o Input.o Matrix.o StatisticalFunctions.o InternalCoords.o FileIO.o Output.o DCDFiles.o DCDHeader.o FreeEnergy.o ArrayUtil.o Input2d.o
Fep2D:    Fep2D.o Analysis.o Output.o StatisticalFunctions.o Input.o Data.o FreeEnergy.o FileIO.o Log.o DCDFiles.o Matrix.o InternalCoords.o DCDHeader.o ArrayUtil.o Util.o Input2D.o

%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.f90*
	$(FC) $(FCFLAGS) -c $<

clean:
	rm -f *.o *.mod *.MOD

install:
	mv Fepcat ../bin
	mv Qfep ../bin
