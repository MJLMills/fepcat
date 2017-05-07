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

  SUBROUTINE Driver

    USE Data,   ONLY : mappingEnergies, energyGap, groundStateEnergy !, lambda
    USE Input,  ONLY : mask, Nbins, minPop
    USE Movies, ONLY : MakeFepusMovie
    USE MovieOptions, ONLY : ProcessNameList

    IMPLICIT NONE

    CALL ProcessNameList
!    CALL MakeFepMovie(mappingEnergies(:,:,:,1),mask(:,:),lambda(:),skip=1)
    CALL MakeFepusMovie(energyGap,groundStateEnergy,mappingEnergies(:,:,:,1),mask,Nbins,minPop,skip=1)

  END SUBROUTINE Driver

END PROGRAM FepMovie
