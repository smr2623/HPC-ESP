 MODULE DFOURIER
  
  IMPLICIT NONE 
  
  INCLUDE 'fftw3.f'
 
  INTEGER, PARAMETER, PRIVATE:: DP=kind(1.0d+00) !double precision kind 
  
  TYPE DFOURIER_TYPE
   PRIVATE
     REAL :: DXINV
     INTEGER :: NLENGTH,NWIDTH
     INTEGER*8, DIMENSION(2) :: I8PLAN
     REAL, DIMENSION(:), POINTER :: WS
     COMPLEX, DIMENSION(:), POINTER :: COUT=>NULL()
     COMPLEX, DIMENSION(:), POINTER :: CIK=>NULL()
  END TYPE DFOURIER_TYPE
 
 PRIVATE :: KW

 CONTAINS
!==================================================

 SUBROUTINE DFOURIER_NEW(THIS,N,NW,D,OK)
  TYPE(DFOURIER_TYPE), INTENT(INOUT) :: THIS
  INTEGER, INTENT(IN) :: N,NW
  REAL, INTENT(IN) :: D
  LOGICAL, INTENT(OUT) :: OK
  !-------------------------
  INTEGER :: NHP1
  
  OK=((N.GT.1).AND.(MOD(N,2).EQ.0))
  IF(OK) THEN
    THIS%NLENGTH=N
	THIS%NWIDTH=NW
	THIS%DXINV=1.0/D
    NHP1=1+N/2
	ALLOCATE(THIS%WS(N))
    ALLOCATE(THIS%COUT(NHP1))
    ALLOCATE(THIS%CIK(NHP1))
    CALL KW(NHP1,THIS%CIK)
    CALL sfftw_plan_dft_r2c_1d(THIS%I8PLAN(1),N,THIS%WS,THIS%COUT,FFTW_MEASURE)  
    CALL sfftw_plan_dft_c2r_1d(THIS%I8PLAN(2),N,THIS%COUT,THIS%WS,FFTW_MEASURE)
  ELSE
   PRINT*, 'ERROR IN DFOURIER_NEW' 
  ENDIF
 
 
 END SUBROUTINE DFOURIER_NEW
!===========================================

 SUBROUTINE DFOURIER_KILL(THIS)
  TYPE(DFOURIER_TYPE), INTENT(INOUT) :: THIS
  
   call sfftw_destroy_plan(THIS%I8PLAN(1))
   call sfftw_destroy_plan(THIS%I8PLAN(2))
   DEALLOCATE(THIS%COUT)
   DEALLOCATE(THIS%CIK)
   DEALLOCATE(THIS%WS)
 
 END SUBROUTINE DFOURIER_KILL
 !===========================================
 
 SUBROUTINE DFOURIER_APPLY(THIS,U,DU)
  TYPE(DFOURIER_TYPE), INTENT(INOUT) :: THIS
  REAL, DIMENSION(:,:), INTENT(IN) :: U
  REAL, DIMENSION(:,:), INTENT(OUT) :: DU
  !-----------------------------------------
  INTEGER :: IW
  
  DO IW=1,THIS%NWIDTH
			
   THIS%WS=U(1:THIS%NLENGTH,IW)		
   call sfftw_execute(THIS%I8PLAN(1))  !fft
   THIS%COUT=THIS%CIK*THIS%COUT !multiply by ik
   call sfftw_execute(THIS%I8PLAN(2))  !backward fft
   DU(1:THIS%NLENGTH,IW)=THIS%WS
   
  ENDDO
 
  !IMPORTANT
  !du is NOW in terms of unit length equal to space sampling step
 
  DU=DU*THIS%DXINV !adjust for the unit length
	
 
 END SUBROUTINE DFOURIER_APPLY
 !===========================================

   SUBROUTINE KW(NHP1,CIK)
   integer, intent(in) :: NHP1
   complex, dimension(NHP1), intent(out) :: CIK
   !-----------------------------------
   real(kind=DP) :: DDK
   integer :: K,N
      
   !width of the space domain is assumed to be N (i.e. DX=1.0)
   !therefore the wavenumber sampling step is 
   ! DK=8.0_DP*atan(1.0_DP)*DINVN 
   
   N=NHP1-1
   N=N+N
   DDK=8.0_DP*atan(1.0_DP)/real(N*N,DP)  !=2PI/(N^2)
  
   CIK(1)=cmplx(0.0,0.0)
  do K=2,NHP1 	                  
      CIK(K)=cmplx(0.0,real(K-1,DP)*DDK)    ! ik     
  enddo
  
  !division by N required for FFT & FFT^{-1}
  !already included in DDK
  
  END SUBROUTINE KW
 !====================================================== 
 
 END MODULE DFOURIER