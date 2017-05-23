MODULE Output

  ! #DES: Subprograms for doing general output

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: WriteCSV2D, WriteXYZ

  CONTAINS

!    SUBROUTINE WriteEVBParameterFile(alpha,A,mu)

      ! #DES: Write an empty EVB parameters file to evb.prm for the number of states in settings.nml

!      USE FileIO, ONLY : OpenFile, CloseFile

!      IMPLICIT NONE
!      REAL(8), INTENT(IN) :: alpha(:), A(:,:), mu(:,:)
!      INTEGER, PARAMETER :: prmUnit = 20
!      INTEGER :: i, j

!      CALL OpenFile(prmUnit,"evb.prm","write")

!      DO i = 1, SIZE(alpha)
!        WRITE(prmUnit,'(F7.3)') alpha(i)
!      ENDDO

!      DO i = 1, SIZE(A,1)
!        DO j = i+1 , SIZE(A,1)
!          WRITE(prmUnit,'(2I4,2F7.3)') i, j, A(i,j), mu(i,j)
!        ENDDO
!      ENDDO

!      CALL CloseFile(prmUnit)

!    END SUBROUTINE WriteEVBParameterFile


    SUBROUTINE WriteCSV2D(head,data,outUnit)

      ! #DES: Write the 2D input 'data' in csv format to unit 'outUnit' as .csv with header 'head'
      ! #DES: 'data' is indexed (row,column)

      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: head(:) !columns
      REAL(8), INTENT(IN) :: data(:,:) !rows, columns
      INTEGER, INTENT(IN) :: outUnit
      INTEGER :: i, j

      DO i = 1, SIZE(head)
        WRITE(outUnit,'(A)',ADVANCE='NO') TRIM(ADJUSTL(head(i)))
        IF (i < SIZE(head)) WRITE(outUnit,'(A1)',ADVANCE='NO') ","
      ENDDO
      WRITE(outUnit,*) ""

      DO i = 1, SIZE(data,1)
        DO j = 1, SIZE(data,2)
          WRITE(outUnit,'(F20.10)',ADVANCE='NO') data(i,j)
          IF(j < SIZE(data,2)) WRITE(outUnit,'(A1)',ADVANCE='NO') ","
        ENDDO
        WRITE(outUnit,*) ""
      ENDDO

    ENDSUBROUTINE

!*

    SUBROUTINE WriteXYZ(title,elements,coordinates,outUnit)

      IMPLICIT NONE
      CHARACTER(*) :: title, elements(:)
      REAL(8), INTENT(IN) :: coordinates(:,:)
      INTEGER, INTENT(IN) :: outUnit
      INTEGER :: atom

      WRITE(outUnit,'(I8)') SIZE(elements)
      WRITE(outUnit,'(A)')  TRIM(ADJUSTL(title))

      DO atom = 1, SIZE(elements)
        WRITE(outUnit,'(A,1X,3F12.7)') TRIM(ADJUSTL(elements(atom))), coordinates(:,atom)
      ENDDO

    END SUBROUTINE WriteXYZ

END MODULE Output
