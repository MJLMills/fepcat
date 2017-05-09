PROGRAM FepMovie

  USE Util, ONLY : Startup, Cleanup
  USE Data, ONLY : ComputeDerivedData
  USE Log,  ONLY : logUnit

  IMPLICIT NONE

    CALL Startup()
    CALL ComputeDerivedData(logUnit,doTiming=.FALSE.,readCoords=.FALSE.,doFEPUS=.TRUE.)
    CALL Driver()
    CALL CleanUp()  

  CONTAINS

!*

  SUBROUTINE Driver()

    USE Data,   ONLY       : mappingEnergies, energyGap, groundStateEnergy, lambda
    USE Input,  ONLY       : mask, Nbins, minPop
    USE FileIO, ONLY       : OpenFile, CloseFile
    USE Movies, ONLY       : MakeFepMovie, MakeFepUsMovie
    USE MovieOptions, ONLY : ProcessNameList, movieOutputDir, plotShellScript

    IMPLICIT NONE
    INTEGER, PARAMETER :: shUnit = 87

    CALL ProcessNameList()
    CALL OpenFile(shUnit,TRIM(ADJUSTL(movieOutputDir))//"/"//TRIM(ADJUSTL(plotShellScript)),"write")
    WRITE(shUnit,'(A)') "rm list.dat"

    CALL MakeFepMovie(mappingEnergies(:,:,:,1),mask(:,:),lambda(:),shUnit)
    CALL MakeFepusMovie(energyGap(:,:),groundStateEnergy(:,:),mappingEnergies(:,:,:,1),mask,Nbins,minPop,shUnit)

    CALL CloseFile(shUnit)

  END SUBROUTINE Driver

END PROGRAM FepMovie
