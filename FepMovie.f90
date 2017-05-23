PROGRAM FepMovie

  USE Util, ONLY : Startup, Cleanup
  USE Data, ONLY : ComputeDerivedData
  USE Log,  ONLY : logUnit

  IMPLICIT NONE

    CALL Startup(readCoords=.FALSE.)
    CALL ComputeDerivedData(logUnit,doTiming=.FALSE.,readCoords=.FALSE.,doFEPUS=.TRUE.)
    CALL Driver()
    CALL CleanUp()

  CONTAINS

!*

  SUBROUTINE Driver()

    USE Data,   ONLY       : mappingEnergies, energyGap, groundStateEnergy, lambda
    USE Input,  ONLY       : mask, Nbins, minPop
    USE FileIO, ONLY       : OpenFile, CloseFile
    USE Movies, ONLY       : SetupMovies, TeardownMovies, MakeFepMovie, MakeFepUsMovie

    IMPLICIT NONE

    CALL SetupMovies()
    CALL MakeFepMovie(mappingEnergies(:,:,:,1),mask(:,:),lambda(:))
    CALL MakeFepusMovie(energyGap(:,:),groundStateEnergy(:,:),mappingEnergies(:,:,:,1),mask(:,:),Nbins,minPop)
    CALL TeardownMovies()

  END SUBROUTINE Driver

END PROGRAM FepMovie

