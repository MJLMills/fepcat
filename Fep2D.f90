PROGRAM Fep2D

  USE Util, ONLY : Startup, Cleanup
  USE Data, ONLY : ComputeDerivedData
  USE Log, ONLY  : logUnit

  IMPLICIT NONE

  CALL Startup()
  CALL ComputeDerivedData(logUnit,doTiming=.FALSE.)
  CALL Run2D()
  CALL CleanUp()  

  CONTAINS

    SUBROUTINE Run2D()

      USE Analysis, ONLY : Test2D
      IMPLICIT NONE

      CALL Test2D()

    END SUBROUTINE Run2D

END PROGRAM Fep2D
