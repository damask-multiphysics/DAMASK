!##############################################################
 MODULE CPFEM
!##############################################################
!    *** CPFEM engine ***
!
 use prec, only: pReal,pInt
 implicit none
!
! ****************************************************************
! *** General variables for the material behaviour calculation ***
! ****************************************************************
 real(pReal), dimension (:,:,:),        allocatable :: CPFEM_cs               ! Cauchy stress
 real(pReal), dimension (:,:,:,:),      allocatable :: CPFEM_dcsdE            ! Cauchy stress tangent
 real(pReal), dimension (:,:,:,:),      allocatable :: CPFEM_dcsdE_knownGood  ! known good tangent

 logical :: CPFEM_init_done    = .false.       ! remember whether init has been done already
 logical :: CPFEM_calc_done    = .false.       ! remember whether first IP has already calced the results

 real(pReal), parameter :: CPFEM_odd_stress = 1e15_pReal, CPFEM_odd_jacobian = 1e50_pReal
!
 CONTAINS
!
!*********************************************************
!***    allocate the arrays defined in module CPFEM    ***
!***    and initialize them                            ***
!*********************************************************
 SUBROUTINE CPFEM_init()

 use prec,           only: pInt,pReal
 use FEsolving,      only: parallelExecution,symmetricSolver,FEsolving_execElem,FEsolving_execIP
 use mesh,           only: mesh_element,mesh_NcpElems,mesh_maxNips,FE_Nips
 use material,       only: homogenization_maxNgrains
 use constitutive,   only: constitutive_maxSizePostResults
 use crystallite,    only: crystallite_Nresults
 use homogenization, only: homogenization_maxSizePostResults

 implicit none
 integer(pInt) e,i,g

 allocate(CPFEM_cs(6,mesh_maxNips,mesh_NcpElems)) ;                CPFEM_cs              = 0.0_pReal
 allocate(CPFEM_dcsdE(6,6,mesh_maxNips,mesh_NcpElems)) ;           CPFEM_dcsde           = 0.0_pReal
 allocate(CPFEM_dcsdE_knownGood(6,6,mesh_maxNips,mesh_NcpElems)) ; CPFEM_dcsde_knownGood = 0.0_pReal


!    *** Output to MARC output file ***
!$OMP CRITICAL (write2out)
 write(6,*)
 write(6,*) '<<<+-  cpfem init  -+>>>'
 write(6,*)
 write(6,'(a32,x,6(i5,x))') 'CPFEM_cs:              ', shape(CPFEM_cs)
 write(6,'(a32,x,6(i5,x))') 'CPFEM_dcsde:           ', shape(CPFEM_dcsde)
 write(6,'(a32,x,6(i5,x))') 'CPFEM_dcsde_knownGood: ', shape(CPFEM_dcsde_knownGood)
 write(6,*)
 write(6,*) 'parallelExecution:    ', parallelExecution
 write(6,*) 'symmetricSolver:      ', symmetricSolver
 call flush(6)
!$OMP END CRITICAL (write2out)
 return
!
 END SUBROUTINE
!
!
!***********************************************************************
!***    perform initialization at first call, update variables and   ***
!***    call the actual material model                               ***
!
!     CPFEM_mode             computation mode (regular, collection, recycle)
!     ffn                    deformation gradient for t=t0
!     ffn1                   deformation gradient for t=t1
!     Temperature            temperature
!     CPFEM_dt               time increment
!     CPFEM_en               element number
!     CPFEM_in               intergration point number
!     CPFEM_stress           stress vector in Mandel notation
!     CPFEM_updateJaco       flag to initiate computation of Jacobian
!     CPFEM_jaco             jacobian in Mandel notation
!     CPFEM_ngens            size of stress strain law
!***********************************************************************
subroutine CPFEM_general(CPFEM_mode, ffn, ffn1, Temperature, CPFEM_dt,&
                         CPFEM_en, CPFEM_in, CPFEM_stress, CPFEM_updateJaco, CPFEM_jaco, CPFEM_ngens)
! note: CPFEM_stress = Cauchy stress cs(6) and CPFEM_jaco = Consistent tangent dcs/de
!
 use prec,         only: pReal,pInt
 use FEsolving
 use debug
 use math
 use mesh,         only: mesh_init,&
                         mesh_FEasCP,mesh_element,mesh_NcpElems,mesh_maxNips,FE_Nips
 use lattice,      only: lattice_init
 use material,     only: material_init, homogenization_maxNgrains
 use constitutive, only: constitutive_init,&
                         constitutive_state0,constitutive_state
 use crystallite
 use homogenization 
 implicit none

 integer(pInt) CPFEM_en, CPFEM_in, cp_en, CPFEM_ngens, i,j,k,l,m,n
 real(pReal), dimension (3,3)        :: ffn,ffn1,Kirchhoff
 real(pReal), dimension (3,3,3,3)    :: H, H_sym
 real(pReal), dimension(CPFEM_ngens) :: CPFEM_stress
 real(pReal), dimension(CPFEM_ngens,CPFEM_ngens) :: CPFEM_jaco
 real(pReal) Temperature,CPFEM_dt,J_inverse
 integer(pInt) CPFEM_mode               ! 1: regular computation with aged results&
                                        ! 2: regular computation&
                                        ! 3: collection of FEM data&
                                        ! 4: recycling of former results (MARC speciality)&
                                        ! 5: record tangent from former converged inc&
                                        ! 6: restore tangent from former converged inc
 integer(pInt) e
 logical CPFEM_updateJaco

 if (.not. CPFEM_init_done) then        ! initialization step (three dimensional stress state check missing?)
   call math_init()
   call FE_init()
   call mesh_init()

   FEsolving_execElem = (/1,mesh_NcpElems/)
   allocate(FEsolving_execIP(2,mesh_NcpElems)); FEsolving_execIP = 1_pInt
   forall (e = 1:mesh_NcpElems) FEsolving_execIP(2,e) = FE_Nips(mesh_element(2,e))

   call lattice_init()
   call material_init()
   call constitutive_init()
   call crystallite_init()
   call homogenization_init()
   call CPFEM_init()
   CPFEM_init_done = .true.
 endif

 cp_en = mesh_FEasCP('elem',CPFEM_en)
 if (cp_en == 1 .and. CPFEM_in == 1) then
    write(6,*) '#####################################'
    write(6,'(a10,1x,f8.4,1x,a10,1x,i4,1x,a10,1x,i3,1x,a10,1x,i2,x,a10,1x,i2)') &
    'theTime',theTime,'theInc',theInc,'theCycle',theCycle,'theLovl',theLovl,&
    'mode',CPFEM_mode
    write(6,*) '#####################################'
 endif

 select case (CPFEM_mode)
    case (1,2)     ! regular computation (with aging of results if mode == 1)
       if (CPFEM_mode == 1) then                  ! age results at start of new increment
         crystallite_F0  = crystallite_partionedF ! crystallite deformation (_subF is perturbed...)
         crystallite_Fp0 = crystallite_Fp         ! crystallite plastic deformation
         crystallite_Lp0 = crystallite_Lp         ! crystallite plastic velocity
         forall (i = 1:homogenization_maxNgrains,&
                 j = 1:mesh_maxNips, &
                 k = 1:mesh_NcpElems) &
           constitutive_state0(i,j,k)%p = constitutive_state(i,j,k)%p   ! microstructure of crystallites
         write(6,'(a10,/,4(3(f10.3,x),/))') 'aged state',constitutive_state(1,1,1)%p/1e6
         do j = 1,mesh_maxNips
           do k = 1,mesh_NcpElems
             if (homogenization_sizeState(j,k) > 0_pInt) &
               homogenization_state0(j,k)%p = homogenization_state(j,k)%p   ! internal state of homogenization scheme
           enddo
         enddo
       endif

       if (outdatedFFN1 .or. any(abs(ffn1 - materialpoint_F(:,:,CPFEM_in,cp_en)) > relevantStrain)) then
         if (.not. outdatedFFN1) write(6,'(a11,x,i5,x,i2,x,a10,/,3(3(f10.3,x),/))') 'outdated at',cp_en,CPFEM_in,'FFN1 now:',ffn1(:,1),ffn1(:,2),ffn1(:,3)
         outdatedFFN1 = .true.
         CPFEM_cs(1:CPFEM_ngens,CPFEM_in,cp_en)                  = CPFEM_odd_stress
         CPFEM_dcsde(1:CPFEM_ngens,1:CPFEM_ngens,CPFEM_in,cp_en) = CPFEM_odd_jacobian*math_identity2nd(CPFEM_ngens)
       else
         if (.not. parallelExecution) then
           FEsolving_execElem(1) = cp_en
           FEsolving_execElem(2) = cp_en
           FEsolving_execIP(1,cp_en) = CPFEM_in
           FEsolving_execIP(2,cp_en) = CPFEM_in
           call materialpoint_stressAndItsTangent(CPFEM_updateJaco, CPFEM_dt)
           call materialpoint_postResults(CPFEM_dt)
         elseif (.not. CPFEM_calc_done) then
           call materialpoint_stressAndItsTangent(CPFEM_updateJaco, CPFEM_dt) ! parallel execution inside
           call materialpoint_postResults(CPFEM_dt)
           CPFEM_calc_done = .true.
         endif

!  translate from P and dP/dF to CS and dCS/dE
         Kirchhoff = math_mul33x33(materialpoint_P(:,:,CPFEM_in, cp_en),transpose(materialpoint_F(:,:,CPFEM_in, cp_en)))
         J_inverse  = 1.0_pReal/math_det3x3(materialpoint_F(:,:,CPFEM_in, cp_en))
         CPFEM_cs(1:CPFEM_ngens,CPFEM_in,cp_en) = math_Mandel33to6(J_inverse*Kirchhoff)

         H = 0.0_pReal
         forall(i=1:3,j=1:3,k=1:3,l=1:3,m=1:3,n=1:3) &
           H(i,j,k,l) = H(i,j,k,l) + &
                          materialpoint_F(j,m,CPFEM_in,cp_en) * &
                          materialpoint_F(l,n,CPFEM_in,cp_en) * &
                          materialpoint_dPdF(i,m,k,n,CPFEM_in,cp_en) - &
                          math_I3(j,l)*materialpoint_F(i,m,CPFEM_in,cp_en)*materialpoint_P(k,m,CPFEM_in,cp_en) + &
                          0.5_pReal*(math_I3(i,k)*Kirchhoff(j,l) + math_I3(j,l)*Kirchhoff(i,k) + &
                                     math_I3(i,l)*Kirchhoff(j,k) + math_I3(j,k)*Kirchhoff(i,l))
         forall(i=1:3,j=1:3,k=1:3,l=1:3) &
            H_sym(i,j,k,l)= 0.25_pReal*(H(i,j,k,l)+H(j,i,k,l)+H(i,j,l,k)+H(j,i,l,k))          ! where to use the symmetric version??
         CPFEM_dcsde(1:CPFEM_ngens,1:CPFEM_ngens,CPFEM_in,cp_en) = math_Mandel3333to66(J_inverse*H)
      endif
    case (3)    ! collect and return odd result
       materialpoint_Temperature(CPFEM_in,cp_en) = Temperature
       materialpoint_F0(:,:,CPFEM_in,cp_en)      = ffn
       materialpoint_F(:,:,CPFEM_in,cp_en)       = ffn1
       CPFEM_cs(1:CPFEM_ngens,CPFEM_in,cp_en)                  = CPFEM_odd_stress
       CPFEM_dcsde(1:CPFEM_ngens,1:CPFEM_ngens,CPFEM_in,cp_en) = CPFEM_odd_jacobian*math_identity2nd(CPFEM_ngens)
       CPFEM_calc_done = .false.

    case (4)    ! do nothing since we can recycle the former results (MARC specialty)
    case (5)    ! record consistent tangent at beginning of new FE increment (while recycling)
       CPFEM_dcsde_knownGood = CPFEM_dcsde
    case (6)    ! restore consistent tangent after FE cutback
       CPFEM_dcsde = CPFEM_dcsde_knownGood
 end select

! return the local stress and the jacobian from storage
 CPFEM_stress(1:CPFEM_ngens) = CPFEM_cs(1:CPFEM_ngens,CPFEM_in,cp_en)
 CPFEM_jaco(1:CPFEM_ngens,1:CPFEM_ngens) = CPFEM_dcsdE(1:CPFEM_ngens,1:CPFEM_ngens,CPFEM_in,cp_en)

 return

end subroutine



 END MODULE
!##############################################################