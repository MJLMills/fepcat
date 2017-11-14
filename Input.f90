!Subprograms for doing task-specific input (read namelists, data)
MODULE Input

  ! #DES: Subprograms for doing task-specific input (read namelist, input files)
  ! #DES: Also stores the input data and provides access to it for other modules

  IMPLICIT NONE

!  PRIVATE
!  PUBLIC :: ReadInput, DeallocateInputArrays

  ! #DES: The namelist settings.nml needs to contain the variables below, and must match the data in the energy files.
  ! Nstates is the number of states in the EVB potential. 
  ! The particular states used to generate the mapping potential are stateA and stateB.
  ! NFepSteps is the number of steps taken from stateA to stateB in the FEP simulation.
  ! Ntimesteps is the number of timesteps in each file.
  ! Nskip tells the program to ignore the first Nskip timesteps in each file.
  ! The final state for LRA study is stepTS.
  ! nBins sets the number of bins to use in the FEP/US histogram.
  ! minPop is the minimum number of points that must be in a FEP/US histogram bin to be considered usable.
  ! fileBase is the part of the filename preceding the FEP index (in I3 format)
  ! temperature is the thermodynamic temperature at which the simulation was run in Kelvin
  ! rcCoeffA and rcCoeffB are used to construct a reaction coordinate from the energy difference
  ! RC = (rcCoeffA * E_A) + (rcCoeffB * E_B)

  NAMELIST /settings/ nStates, stateA, stateB, nFepSteps, nSkip, stepTS, nBins, minPop, fileBase, trajecBase, &
                      temperature, rcCoeffA, rcCoeffB, useEnergyGapCoupling, dGTS, dGPS, outDir
  INTEGER :: nFepSteps, nStates, nSkip, stateA, stateB, stepTS, nBins, minPop
  REAL(8) :: temperature, rcCoeffA, rcCoeffB, dGTS, dGPS
  LOGICAL :: useEnergyGapCoupling
  CHARACTER(100) :: fileBase, trajecBase, outDir

  ! #DES: The stateEnergy and coeffs arrays hold all of the simulation data (irrespective of nSkip, nTimesteps)
  REAL(8), ALLOCATABLE :: stateEnergy(:,:,:,:), coeffs(:,:,:)
  ! The trajectory array holds all the Cartesian coordinates of the set of FEP simulations - Single Precision cos of the dcd format.
  REAL(4), ALLOCATABLE :: trajectory(:,:,:,:)
  ! #DES: These arrays hold the EVB parameters
  REAL(8), ALLOCATABLE :: alpha(:), couplingConstant(:,:), couplingGaussExpFactor(:,:), couplingExpExpFactor(:,:)
  REAL(8) :: beta
  LOGICAL :: targetsPresent = .FALSE.

  ! This data is used to determine what data in stateEnergy and coeffs is to be used in the FEP/US procedure
  ! The block mask(fep,skip(fep)+1:nTimesteps(fep)) must be set to .TRUE. initially, the rest to .FALSE.

  INTEGER, ALLOCATABLE :: nTimesteps(:) ! nFepSteps - number of actual timesteps in each FEP simulation
  INTEGER, ALLOCATABLE :: skip(:)       ! nFepSteps - How many steps to ignore at the start of each step; allows a per-simulation equilibration period
  LOGICAL, ALLOCATABLE :: mask(:,:)     ! nFepSteps x nTimesteps - pointwise mask for simulation data; turn points on/off
  INTEGER :: maxTimesteps               ! Length of the longest FEP simulation (largest of nTimesteps(:)); used to allocate arrays

  ! These will eventually have to be inputs unless a logical switch for Q data is TRUE
  CHARACTER(20), ALLOCATABLE :: energyNames(:)
  INTEGER :: nEnergyTypes

  CONTAINS

    SUBROUTINE AllocateInputArrays

      ! #DES: Allocate and zero the input arrays based on the supplied namelist inputs

      IMPLICIT NONE

      ALLOCATE(stateEnergy(maxTimesteps,nFepSteps,nStates,nEnergyTypes)); stateEnergy = 0.0d0
      ALLOCATE(skip(nFepSteps));                                          skip(:) = 0
      ALLOCATE(mask(nFepSteps,maxTimesteps));                             mask(:,:) = .TRUE.
      ALLOCATE(coeffs(maxTimesteps,nFepSteps,nStates));                   coeffs(:,:,:) = 0.0d0

      ALLOCATE(alpha(nStates));                                           alpha(:)                    = 0.0d0
      ALLOCATE(couplingConstant(nStates,nStates));                        couplingConstant(:,:)       = 0.0d0
      ALLOCATE(couplingGaussExpFactor(nStates,nStates));                  couplingGaussExpFactor(:,:) = 0.0d0
      ALLOCATE(couplingExpExpFactor(nStates,nStates));                    couplingExpExpFactor(:,:)   = 0.0d0

    END SUBROUTINE AllocateInputArrays

!*

    SUBROUTINE DeallocateInputArrays

      ! #DES: Free the memory associated with the input arrays

      IMPLICIT NONE

      IF (ALLOCATED(stateEnergy))             DEALLOCATE(stateEnergy)
      IF (ALLOCATED(coeffs))                  DEALLOCATE(coeffs)
      IF (ALLOCATED(alpha))                   DEALLOCATE(alpha)
      IF (ALLOCATED(couplingConstant))        DEALLOCATE(couplingConstant)
      IF (ALLOCATED(couplingGaussExpFactor))  DEALLOCATE(couplingGaussExpFactor)
      IF (ALLOCATED(couplingExpExpFactor))    DEALLOCATE(couplingExpExpFactor)
      IF (ALLOCATED(skip))                    DEALLOCATE(skip)
      IF (ALLOCATED(nTimesteps))              DEALLOCATE(nTimesteps)
      IF (ALLOCATED(mask))                    DEALLOCATE(mask)
      IF (ALLOCATED(trajectory))              DEALLOCATE(trajectory)

    END SUBROUTINE DeallocateInputArrays

!*

    SUBROUTINE ReadInput(readCoords)

      ! #DES: Master routine for reading all needed input from files

      USE Log, ONLY : logUnit
      IMPLICIT NONE
      LOGICAL, INTENT(IN) :: readCoords

      CALL ProcessNamelist(logUnit)
      CALL ReadEnergyTypes(logUnit)
      CALL ReadSimulationLengths()
      CALL WriteSimulationLengths(logUnit)
      CALL AllocateInputArrays()
      CALL ReadEVBParameters(logUnit)
      CALL PrintEVBParameters(logUnit)
      CALL ReadConvergenceParameters(logUnit)
      CALL WriteConvergenceParameters(logUnit)
      CALL ReadFormattedEnergyFiles(logUnit)
      IF (readCoords .EQV. .TRUE.) CALL ReadCoordinateFiles(logUnit)

    END SUBROUTINE ReadInput

!*

    SUBROUTINE ReadSimulationLengths

      ! #DES: Determine the number of timesteps in each simulation

      USE FileIO, ONLY : fileLengths
      IMPLICIT NONE
      CHARACTER(500) :: fileNames(nFepSteps)
      INTEGER :: step

      ALLOCATE(nTimeSteps(nFepSteps))
      CALL CreateFileNames('en',fileBase,fileNames)
      nTimesteps(:) = fileLengths(fileNames)
      nTimeSteps(:) = nTimesteps(:) / nStates
      DO step = 1, nFepSteps
        IF (ALLOCATED(mask)) mask(step,nTimesteps(step)+1:maxTimesteps) = .FALSE.
      ENDDO
      maxTimesteps = MAXVAL(nTimesteps)

      IF (nSkip < 0 .OR. nSkip >= maxTimesteps-1) THEN
        WRITE(*,*) "ReadSimulationLength - Illegal value of nSkip: ", nSkip
        IF (nSkip < 0) STOP "nSkip < 0"
        IF (nSkip >= maxTimesteps-1) STOP "nSkip >= maxTimesteps - 1"
      ENDIF
      IF (minPop < 1 .OR. minPop > maxTimesteps)  STOP ": ReadSimulationLength: Illegal value of minPop"

    END SUBROUTINE ReadSimulationLengths

!*

    SUBROUTINE WriteSimulationLengths(logUnit)

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: logUnit
      CHARACTER(500) :: fileNames(nFepSteps)
      INTEGER :: step

      CALL CreateFileNames('en',fileBase,fileNames)
      WRITE(logUnit,'(A)') "* Simulation Lengths *"; WRITE(logUnit,*)
      DO step = 1, nFepSteps
        WRITE(logUnit,'(A,I8)') TRIM(ADJUSTL(fileNames(step))), nTimesteps(step)
      ENDDO
      WRITE(logUnit,*)
      WRITE(logUnit,'(A,I8)') "Longest Simulation Length: ", maxTimesteps
      WRITE(logUnit,*)

    END SUBROUTINE WriteSimulationLengths

!*

    SUBROUTINE WriteConvergenceParameters(logUnit)

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: logUnit
      CHARACTER(500) :: fileNames(nFepSteps)
      INTEGER :: step

      WRITE(logUnit,'(A)') "Convergence parameters:"
      WRITE(logUnit,*)
      WRITE(logUnit,'(A)') "Energy File  Fep Step   Skip"
      CALL createFileNames('en',fileBase,fileNames)
      DO step = 1, nFepSteps
        WRITE(logUnit,'(A,2I8)')  TRIM(ADJUSTL(fileNames(step))), step, skip(step)
      ENDDO
      WRITE(logUnit,*)

    END SUBROUTINE WriteConvergenceParameters

!*

    SUBROUTINE ReadConvergenceParameters(logUnit)

      USE FileIO, ONLY : OpenFile, CloseFile

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: logUnit
      INTEGER, PARAMETER :: prmUnit = 16
      INTEGER :: step
      LOGICAL :: fileExists
      CHARACTER(100) :: discard

      WRITE(logUnit,'(A)', ADVANCE='NO') "Determining FEP Simulation Convergence Parameters... "
      CALL OpenFile(prmUnit,"skip.prm","read",success=fileExists)
      IF (fileExists) THEN
        DO step = 1, nFepSteps
          READ(prmUnit,*) discard, skip(step)
          IF (skip(step) > nTimesteps(step)) STOP "ERROR: nSkip > nTimesteps (ReadConvergenceParameters)"
          mask(step,1:skip(step)) = .FALSE.
        ENDDO
        CALL CloseFile(prmUnit)
        WRITE(logUnit,'(A)') "data read from skip.prm"
        WRITE(logUnit,*)
      ELSE
        skip(:) = nSkip
        mask(:,1:nSkip) = .FALSE.
        CALL WriteConvergenceParamFile()
        WRITE(logUnit,'(A)') "no convergence parameter file provided; template written to skip.prm"
        WRITE(logUnit,*)
      ENDIF

    END SUBROUTINE ReadConvergenceParameters

!*

    SUBROUTINE WriteConvergenceParamFile()

      ! #DES: Write an empty convergence parameter file to skip.prm for the number of FEP simulations

      USE FileIO, ONLY : OpenFile, CloseFile

      IMPLICIT NONE
      INTEGER, PARAMETER :: prmUnit = 16
      CHARACTER(500) :: fileNames(nFepSteps)
      INTEGER :: step

      CALL OpenFile(prmUnit,"skip.prm","write")
      CALL CreateFileNames('en',fileBase,fileNames)

      DO step = 1, nFepSteps
        WRITE(prmUnit,'(A,I8)') TRIM(ADJUSTL(fileNames(step))), skip(step)
      ENDDO

      CALL CloseFile(prmUnit)

    END SUBROUTINE WriteConvergenceParamFile

!*

    SUBROUTINE ReadEVBParameters(logUnit)

      ! #DES: Read a text file containing diagonal constants and of-diagonal function parameters
      ! #DES: If this cannot be found, run the calculation with all zeroes, dump a .prm file and note in the log file

      USE FileIO, ONLY : OpenFile, CloseFile

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: logUnit
      INTEGER, PARAMETER :: prmUnit = 20
      INTEGER :: state, i, j, dummy
      LOGICAL :: fileExists

      WRITE(logUnit,'(A34)',ADVANCE='NO') "Determining EVB Parameter Data... "

      CALL OpenFile(prmUnit,"evb.prm","read",success=fileExists)

      IF (fileExists) THEN

        DO state = 1, nStates
          READ(prmUnit,*) alpha(state)
        ENDDO

        DO i = 1, nStates
          DO j = i+1, nstates
            READ(prmUnit,*) dummy, dummy, couplingConstant(i,j), couplingExpExpFactor(i,j), couplingGaussExpFactor(i,j)
            couplingConstant(j,i)     = couplingConstant(i,j)
            couplingExpExpFactor(j,i) = couplingExpExpFactor(i,j)
            couplingGaussExpFactor(j,i) = couplingGaussExpFactor(i,j)
          ENDDO
        ENDDO

        CALL CloseFile(prmUnit)
        WRITE(logUnit,'(A22)') "data read from evb.prm"
        WRITE(logUnit,*)

      ELSE

        IF (targetsPresent .EQV. .TRUE.) THEN
          WRITE(logUnit,'(A)') "no EVB parameter file provided: using targets to guess EVB parameters"
          WRITE(logUnit,*)
        ELSE
          CALL WriteEmptyEVBParameterFile()
          WRITE(logUnit,'(A)') "no EVB parameter file or free energy targets provided; template written to evb.prm"
          WRITE(logUnit,*)
        ENDIF

      ENDIF

    END SUBROUTINE ReadEVBParameters

!*

    SUBROUTINE PrintEVBParameters(logUnit)

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: logUnit
      INTEGER :: state, i, j

      WRITE(logUnit,'(A15)') "EVB parameters:"
      WRITE(logUnit,*)

      WRITE(logUnit,'(A22)') "* Diagonal Constants *"
      WRITE(logUnit,'(A13)') "State   Alpha"
      DO state = 1, nStates
        WRITE(logUnit,'(I3,3X,F7.3)') state, alpha(state)
      ENDDO
      WRITE(logUnit,*)

      WRITE(logUnit,'(A27)') "* Off-Diagonal Parameters *"
      WRITE(logUnit,'(A)') "i   j   A   mu   eta"
      DO i = 1, nStates
        DO j = i+1, nStates
          WRITE(logUnit,'(2I4,3F7.3)') i, j, couplingConstant(i,j), couplingExpExpFactor(i,j), couplingGaussExpFactor(i,j)
        ENDDO
      ENDDO
      WRITE(logUnit,*)

      WRITE(logUnit,*) "* Coupling Coordinate *"
      IF (useEnergyGapCoupling) WRITE(logUnit,'(A)') "Using energy gap as coupling coordinate"

    END SUBROUTINE PrintEVBParameters

!*

    SUBROUTINE WriteEmptyEVBParameterFile

      ! #DES: Write an empty EVB parameters file to evb.prm for the number of states in settings.nml

      USE FileIO, ONLY : OpenFile, CloseFile

      IMPLICIT NONE
      INTEGER, PARAMETER :: prmUnit = 20
      INTEGER :: i, j

      CALL OpenFile(prmUnit,"evb.prm","write")

      DO i = 1, nStates
        WRITE(prmUnit,'(F7.3)') 0.0d0
      ENDDO

      DO i = 1, nStates
        DO j = i+1 , nStates
          WRITE(prmUnit,'(2I4,3F7.3)') i, j, 0.0d0, 0.0d0, 0.0d0
        ENDDO 
      ENDDO

      CALL CloseFile(prmUnit)

    END SUBROUTINE WriteEmptyEVBParameterFile

!*

    SUBROUTINE ProcessNamelist(logUnit)

      ! #DES: Read the namelist settings.nml, check the inputs, determine dependent variables and print to log file

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: logUnit
      CHARACTER(100) :: command

      CALL ReadNameList()
      CALL CheckNameList()

      !Compute dependent information
      beta = 1.0d0 / (0.0019872041d0 * temperature)
      IF (dGTS /= HUGE(0.0d0) .AND. dGPS /= HUGE(0.0d0)) targetsPresent = .TRUE.
      command = "mkdir "//TRIM(ADJUSTL(outdir))
      !CALL SYSTEM(command) ! Write functional directory creation routine

      CALL PrintNameList(logUnit)
      
    END SUBROUTINE ProcessNamelist

!*

    SUBROUTINE ReadNamelist

      ! #DES: Open, read and close the namelist file settings.nml

      USE FileIO, ONLY : OpenFile, CloseFile
      IMPLICIT NONE
      INTEGER, PARAMETER :: nmlUnit = 20

      CALL OpenFile(nmlUnit,"settings.nml","read")

      nStates     = 0
      stateA      = 0
      stateB      = 0
      nFepSteps   = 0
      stepTS      = 0
      temperature = -1.0d0

      ! useable defaults
      nSkip      =  0
      rcCoeffA   =  1.0d0
      rcCoeffB   = -1.0d0
      nBins      =  0
      minPop     =  1
      fileBase   =  "EnergyFile"
      outDir     =  "output"
      useEnergyGapCoupling = .TRUE.
      dGTS = HUGE(0.0d0)
      dGPS = HUGE(0.0d0)

      READ(nmlUnit,NML=settings)
      CALL CloseFile(nmlUnit)

  END SUBROUTINE ReadNameList

!*

  SUBROUTINE CheckNameList

    ! #DES: Confirm that all input values are in the correct bounds

    IMPLICIT NONE

    IF (nStates < 2)                            STOP ": ReadNameList: ILLEGAL NUMBER OF EVB STATES REQUESTED"
    IF (stateA < 1 .OR. stateA > nStates)       STOP ": ReadNameList: Illegal value of stateA"
    IF (stateB < 1 .OR. stateB > nStates)       STOP ": ReadNameList: Illegal value of stateB"
    IF (nFepSteps  < 2)                         STOP ": ReadNameList: ILLEGAL NUMBER OF FEP STEPS REQUESTED"
    IF (stepTS < 1 .OR. stepTS > nFepSteps)     STOP ": ReadNameList: Illegal value of StepTS"
    IF (nBins < 1)                              STOP ": ReadNameList: Illegal value of nBins"
    IF (fileBase == "")                         STOP ": ReadNameList: Illegal value of fileBase"
    IF (temperature < 0.0d0)                    STOP ": ReadNameList: Illegal value of temperature"

  END SUBROUTINE CheckNameList

!*

  SUBROUTINE PrintNameList(logUnit)

    ! #DES: Reprint all namelist input to the log file

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: logUnit

    WRITE(logUnit,'(A)') "FEP/US parameters read from settings.nml:"
    WRITE(logUnit,*)
    WRITE(logUnit,'(A)') "* EVB Settings *"
    WRITE(logUnit,'(A,I6)')   "Number of EVB States         : ", nStates
    WRITE(logUnit,'(A,I6)')   "EVB State a                  : ", stateA
    WRITE(logUnit,'(A,I6)')   "EVB State b                  : ", stateB
    WRITE(logUnit,*) 
    WRITE(logUnit,'(A)') "* MD/FEP Settings *"
    WRITE(logUnit,'(A,I6)')   "Number of Mapping Simulations: ", nFepSteps
    WRITE(logUnit,'(A,I6)')   "Simulation Steps to Skip     : ", nSkip
    WRITE(logUnit,'(A,F6.2)') "Simulation Temperature       : ", temperature
    WRITE(logUnit,'(A,A)')    "Energy File Prefix           : ", TRIM(ADJUSTL(fileBase))
    WRITE(logUnit,*)
    WRITE(logUnit,'(A)') "* FEP/US Settings *"
    WRITE(logUnit,'(A,I6)')   "Step for LRA calculations    : ", stepTS
    WRITE(logUnit,'(A,I6)')   "Number of FEP/US Bins        : ", nBins
    WRITE(logUnit,'(A,I6)')   "Population Threshold of Bins : ", minPop
    WRITE(logUnit,'(A,L6)')   "Use Energy Gap for Coupling  : ", useEnergyGapCoupling
    WRITE(logUnit,'(A,A)')    "Output Directory             : ", outDir

    IF (targetsPresent .EQV. .TRUE.) THEN
      WRITE(logUnit,'(A,F6.2)') "Target Free Energy RS -> TS  : ", dGTS
      WRITE(logUnit,'(A,F6.2)') "Target Free Energy RS -> PS  : ", dGPS
      WRITE(logUnit,*)
    ENDIF

    WRITE(logUnit,*)

  END SUBROUTINE PrintNameList

!*

  SUBROUTINE CreateFileNames(ext,prefix,fileNames)

    IMPLICIT NONE
    CHARACTER(*), INTENT(IN) :: ext, prefix
    CHARACTER(500), INTENT(OUT) :: fileNames(nFepSteps)
    CHARACTER(3) :: stepString
    INTEGER :: step

    DO step = 1, nFepSteps
      WRITE(stepString,'(I0.3)') step - 1
      fileNames(step) = TRIM(ADJUSTL(prefix))//stepString//'.'//TRIM(ADJUSTL(ext))
    ENDDO

  END SUBROUTINE CreateFileNames

!*

  SUBROUTINE ReadCoordinateFiles(logUnit)

    ! #DES: Read the trajectory from the .dcd files corresponding to each .en file

    USE DCDFiles, ONLY : ReadDCD, CountAtoms

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: logUnit
    CHARACTER(500) :: fileNames(nFepSteps)
    INTEGER :: step, nAtoms, maxAtoms

    WRITE(logUnit,'(A)') "Reading Simulation Trajectory Data"
    WRITE(logUnit,*)
    CALL CreateFileNames('dcd',trajecBase,fileNames)

    maxAtoms = 0
    DO step = 1, NFepSteps
      CALL CountAtoms(fileNames(step),nAtoms)
      IF (nAtoms > maxAtoms) maxAtoms = nAtoms
    ENDDO

    ALLOCATE(trajectory(nFepSteps,3,maxAtoms,maxTimesteps))

    DO step = 1, nFepSteps
      WRITE(logUnit,*) "Reading trajectory from ", TRIM(ADJUSTL(fileNames(step)))
      CALL ReadDCD(fileNames(step),trajectory(step,:,:,:))
    ENDDO

    WRITE(logUnit,*)
    WRITE(logUnit,'(A)') "Finished Reading Simulation Trajectory Data"
    WRITE(logUnit,*)

  END SUBROUTINE ReadCoordinateFiles

!*

  SUBROUTINE ReadFormattedEnergyFiles(logUnit)

    ! #DES: Read the nFepsSteps energy files into the arrays stateEnergy and coeffs
    ! Need to edit this such that it does not need to know beforehand how many timesteps are in each fep file
    ! and if a fep file is missing/incomplete, instead of dying just have however many points were read OK in stateEnergy(fepstep,:)
    ! This saves knowing the length of simulations, allows analysis on partially completed simulations.

    USE FileIO, ONLY : OpenFile, CloseFile, fileLengths

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: logUnit
    INTEGER, PARAMETER :: eneUnit = 10
    CHARACTER(500) :: fileNames(nFepSteps)
    INTEGER :: step, state, timestep, ios, type
    LOGICAL :: openSuccess

    WRITE(logUnit,*) ""
    WRITE(logUnit,'(A30)') "Reading Simulation Energy Data"
    WRITE(logUnit,*) ""

    CALL CreateFileNames('en',fileBase,fileNames)

    ReadSteps: DO step = 1, nFepSteps

      CALL OpenFile(eneUnit,fileNames(step),"read",success=openSuccess)
      IF (nTimesteps(step) == 0) CYCLE ReadSteps

      ! This would be the place to allocate space for this fep simulation's energy data
      DO timestep = 1, nTimesteps(step)
        DO state = 1, nStates

          READ(eneUnit,'(F15.8)',ADVANCE='NO',IOSTAT=ios) coeffs(timestep,step,state)
          IF (ios /= 0) THEN
            WRITE(logUnit,'(A,I0.2,A,I0.6,A,I0.4)') "Failed to read expected coefficient for state ", state, " at timestep ", timestep, " of FEP simulation ", step 
            WRITE(logUNit,'(A)') "Moving to next FEP simulation"
            CYCLE ReadSteps
          ENDIF

          DO type = 1, nEnergyTypes
            READ(eneUnit,'(F15.8)',ADVANCE='NO',IOSTAT=ios) stateEnergy(timestep,step,state,type)
          ENDDO; READ(eneUnit,*)

        ENDDO
      ENDDO
      CALL CloseFile(eneUnit)

      WRITE(logUnit,'(A4,1X,I0.6,1X,A14,1X,A)') "Read", timestep-1, "timesteps from", TRIM(ADJUSTL(fileNames(step)))

    ENDDO ReadSteps

    WRITE(logUnit,*)
    WRITE(logUnit,'(A39)') "Finished Reading Simulation Energy Data"
    WRITE(logUnit,*)

  END SUBROUTINE ReadFormattedEnergyFiles

!*

  SUBROUTINE ReadEnergyTypes(logUnit)

    USE FileIO, ONLY : OpenFile, CloseFile, FileLength
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: logUnit
    INTEGER, PARAMETER :: eneUnit = 22
    LOGICAL :: fileExists = .FALSE.
    INTEGER :: line

    WRITE(logUnit,'(A)') "* Determine Energy Types *"; WRITE(logUnit,*)
    CALL OpenFile(eneUnit,"EnergyNames.dat","read",success=fileExists)

    IF (fileExists .EQV. .TRUE.) THEN

      WRITE(logUnit,'(A)') "Reading energy names from EnergyNames.dat"
      ! This is a fudge to get around the way FileLength is written - add a 2nd version that just takes a unit #, counts, then rewinds
      CALL CloseFile(eneUnit)
      nEnergyTypes = FileLength("EnergyNames.dat")
      CALL OpenFile(eneUnit,"EnergyNames.dat","read",success=fileExists)

      ALLOCATE(energyNames(nEnergyTypes))
      DO line = 1, nEnergyTypes
        READ(eneUnit,'(A)') EnergyNames(line)
        IF (TRIM(ADJUSTL(energyNames(line))) == "") STOP "Error: Remove blank spaces in EnergyNames.dat"
      ENDDO

    ELSE

      WRITE(logUnit,'(A)') "No file EnergyNames.dat - using Q Defaults"
      ! No file so assume Q types
      CALL DefineQTypes()

    ENDIF

    WRITE(logUnit,'(A,1X,I4)') "Number of Types: ", nEnergyTypes
    DO line = 1, nEnergyTypes
      WRITE(logUnit,'(A)') trim(adjustl(energyNames(line)))
    ENDDO
    WRITE(logUnit,*)

  END SUBROUTINE ReadEnergyTypes
!*

  SUBROUTINE DefineQTypes

    ! #DES: Defines the energy types as found in Q output files for convenience

    IMPLICIT NONE

    nEnergyTypes = 14
    ALLOCATE(energyNames(nEnergyTypes))

    energyNames(1)  = "    Q_TOTAL"; energyNames(2)  = "     Q_BOND";
    energyNames(3)  = "    Q_ANGLE"; energyNames(4)  = "  Q_TORSION";
    energyNames(5)  = " Q_IMPROPER"; energyNames(6)  = "     Q_ELEC";
    energyNames(7)  = "      Q_VDW"; energyNames(8)  = "   Q-Q_ELEC";
    energyNames(9)  = "    Q-Q_VDW"; energyNames(10) = "   Q-P_ELEC";
    energyNames(11) = "    Q-P_VDW"; energyNames(12) = "   Q-W_ELEC";
    energyNames(13) = "    Q-W_VDW"; energyNames(14) = "    Q_RESTR";

  END SUBROUTINE DefineQTypes

!*

  SUBROUTINE WriteEnergyTypes(fileUnit)

    ! #DES: Write the names of the energy types to fileUnit in .csv format

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: fileUnit
    INTEGER :: type

    DO type = 1, SIZE(energyNames)
      WRITE(fileUnit,'(A11,A1)',ADVANCE='no') adjustr(energyNames(type)), ","
    ENDDO

  END SUBROUTINE WriteEnergyTypes

END MODULE Input
