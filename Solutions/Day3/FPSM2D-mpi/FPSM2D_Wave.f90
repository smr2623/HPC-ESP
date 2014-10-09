MODULE FPSM2D_WAVE
!PURPOSE:
!HOLD THE TEMPORARY STATE OF THE SIMULATION
!UPDATE IT AT EACH TIME STEP USING 2ND ORDER FD
!AND SAVE IT WHEN REQUIRED

 USE FPSM2D_PINPUT
 USE FPSM2D_OPERATOR
 USE FPSM2D_OUTPUT
 
 IMPLICIT NONE
 
 TYPE FPSM2D_WAVE_TYPE
  PRIVATE
   INTEGER :: INEW,IOLD,INOW
   TYPE(FPSM2D_OPERATOR_TYPE) :: L
   TYPE(FPSM2D_OUTPUT_TYPE) :: SAFE 
   REAL, DIMENSION(:,:,:,:), POINTER :: U => NULL()
    REAL, DIMENSION(:,:,:), POINTER :: DU
 END TYPE FPSM2D_WAVE_TYPE

 CONTAINS
!===========================================
 
 SUBROUTINE FPSM2D_WAVE_NEW(THIS,INPUT,OK)
  TYPE(FPSM2D_WAVE_TYPE), INTENT(INOUT) :: THIS
  TYPE(FPSM2D_PINPUT_TYPE), INTENT(IN) :: INPUT
  LOGICAL, INTENT(OUT) :: OK
  !------------------------------
  INTEGER, DIMENSION(2) :: N,ND
  
  CALL FPSM2D_PINPUT_GETSPACE(INPUT,N,ND)
  
  ALLOCATE(THIS%U(N(1),ND(2),2,3))
  ALLOCATE(THIS%DU(N(1),ND(2),2))
  
  CALL FPSM2D_OPERATOR_NEW(THIS%L,INPUT,OK)
  
  IF(OK) CALL FPSM2D_OUTPUT_NEW(THIS%SAFE,INPUT,OK)
  
  IF(OK) THEN
	THIS%INEW=1
	THIS%IOLD=2
	THIS%INOW=3
  
	CALL FPSM2D_PINPUT_GETINITIAL(INPUT,&
                                  THIS%U(:,:,:,THIS%IOLD),&
								  THIS%U(:,:,:,THIS%INOW))
  ENDIF
  
  IF(.NOT.OK) PRINT*, 'ERROR IN FPSM2D_WAVE_NEW' 
 
 END SUBROUTINE FPSM2D_WAVE_NEW
!===========================================

 SUBROUTINE FPSM2D_WAVE_KILL(THIS)
  TYPE(FPSM2D_WAVE_TYPE), INTENT(INOUT) :: THIS
  
  CALL FPSM2D_OUTPUT_KILL(THIS%SAFE)
  CALL FPSM2D_OPERATOR_KILL(THIS%L)
  
  DEALLOCATE(THIS%U)
  DEALLOCATE(THIS%DU)
 
 END SUBROUTINE FPSM2D_WAVE_KILL
!===========================================

  SUBROUTINE FPSM2D_WAVE_UPDATE(THIS)
  !CALLED AT EACH TIME STEP
  !2ND ORDER FD ADVANCING
   TYPE(FPSM2D_WAVE_TYPE), INTENT(INOUT) :: THIS
   !----------------------
   INTEGER :: IBUF
   
   CALL FPSM2D_OPERATOR_APPLY(THIS%L,THIS%U(:,:,:,THIS%INOW),THIS%DU)
     	  
   !THIS%DU IS ALREADY MULTIPLIED WITH DT*DT
   THIS%U(:,:,:,THIS%INEW)=THIS%U(:,:,:,THIS%INOW)+THIS%U(:,:,:,THIS%INOW)-&
                         THIS%U(:,:,:,THIS%IOLD)+THIS%DU
  
   IBUF=THIS%IOLD
   THIS%IOLD=THIS%INOW
   THIS%INOW=THIS%INEW
   THIS%INEW=IBUF 
   
  END SUBROUTINE FPSM2D_WAVE_UPDATE
!===========================================

  SUBROUTINE FPSM2D_WAVE_OUTPUT(THIS)
   TYPE(FPSM2D_WAVE_TYPE), INTENT(INOUT) :: THIS
   
    CALL FPSM2D_OUTPUT_SAVE(THIS%SAFE,THIS%U(:,:,:,THIS%INOW))
  
  END SUBROUTINE FPSM2D_WAVE_OUTPUT
!===========================================
 
END MODULE FPSM2D_WAVE