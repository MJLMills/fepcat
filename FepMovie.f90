PROGRAM FepMovie

  USE Util, ONLY : Startup, Cleanup
  USE Data, ONLY : ComputeDerivedData
  USE Log,  ONLY : logUnit

  IMPLICIT NONE

    CALL Startup()
    CALL ComputeDerivedData(logUnit,doTiming=.FALSE.,readCoords=.FALSE.,doFEPUS=.FALSE.)
    CALL Driver()
    CALL CleanUp()  

  CONTAINS

!*

  SUBROUTINE Driver

    USE Data,   ONLY : mappingEnergies, lambda
    USE Input,  ONLY : mask
    USE Movies, ONLY : MakeFepMovie
    IMPLICIT NONE

    CALL MakeFepMovie(mappingEnergies(:,:,:,1),mask(:,:),lambda(:),skip=100)

  END SUBROUTINE Driver

END PROGRAM FepMovie
