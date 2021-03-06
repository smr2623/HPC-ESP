 MODULE FPSM2D_OUTPUT
!PURPOSE: SAVE SIMULATION STATUS AT GIVEN TIME
!BY COLLECTING THE DATA FROM ALL PROCESSES
!AND WRITING TO A SINGLE FILE USING THE ROOT PROCESS

 USE FPSM2D_PINPUT
 USE FPSM2D_PARALLELISM
 USE FPSM2D_FILEUNITS
 
 IMPLICIT NONE
 
 TYPE FPSM2D_OUTPUT_TYPE
  PRIVATE
  LOGICAL :: LROOT
  INTEGER :: NBLOCKSIZE
  REAL, DIMENSION(:,:), POINTER :: PLANE
 END TYPE FPSM2D_OUTPUT_TYPE
 
 CONTAINS
!===============================================

 SUBROUTINE FPSM2D_OUTPUT_NEW(THIS,INPUT,OK)
  TYPE(FPSM2D_OUTPUT_TYPE), INTENT(INOUT) :: THIS
  TYPE(FPSM2D_PINPUT_TYPE), INTENT(IN) :: INPUT
  LOGICAL, INTENT(OUT) :: OK
  !.................................
  INTEGER :: ME,IERR
  INTEGER, DIMENSION(2) :: N,ND
  
   CALL FPSM2D_PINPUT_GETSPACE(INPUT,N,ND)
   THIS%NBLOCKSIZE=N(1)*ND(2)
   
   !LROOT  
   CALL MPI_COMM_RANK(MPI_COMM_WORLD,ME,IERR) 
   THIS%LROOT=ME.EQ.0
  
   IF(THIS%LROOT) THEN
    ALLOCATE(THIS%PLANE(N(1),N(2)))
    OPEN(FPSM2D_FILEUNITS_OUTPUT,FILE=FPSM2D_FILEUNITS_OUTPUT_FNAME,&
	     ACCESS='STREAM',FORM='UNFORMATTED')
    print*,'output file is open'
   ENDIF
   
 END SUBROUTINE FPSM2D_OUTPUT_NEW
!===============================================

 SUBROUTINE FPSM2D_OUTPUT_SAVE(THIS,UNOW)
  TYPE(FPSM2D_OUTPUT_TYPE), INTENT(IN) :: THIS
  REAL, DIMENSION(:,:,:), INTENT(IN) :: UNOW
  !----------------------------------------
  INTEGER :: IC,IERR
  
  DO IC=1,2
   CALL MPI_GATHER(UNOW(:,:,IC),THIS%NBLOCKSIZE,MPI_REAL,&
                   THIS%PLANE,THIS%NBLOCKSIZE,MPI_REAL,0,MPI_COMM_WORLD,IERR)
   IF(THIS%LROOT) WRITE(FPSM2D_FILEUNITS_OUTPUT) THIS%PLANE				   
  ENDDO
  
   IF(THIS%LROOT) print*,'partial results sent to output file'
  
 
 END SUBROUTINE FPSM2D_OUTPUT_SAVE
!===============================================

 SUBROUTINE FPSM2D_OUTPUT_KILL(THIS)
  TYPE(FPSM2D_OUTPUT_TYPE), INTENT(IN) :: THIS
  
  IF(THIS%LROOT) THEN
   CLOSE(FPSM2D_FILEUNITS_OUTPUT)
   DEALLOCATE(THIS%PLANE)
  ENDIF
  
 
 END SUBROUTINE FPSM2D_OUTPUT_KILL
!===============================================

END MODULE FPSM2D_OUTPUT
 