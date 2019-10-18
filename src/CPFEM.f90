!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief CPFEM engine
!--------------------------------------------------------------------------------------------------
module CPFEM
  use prec
  use numerics
  use debug
  use FEsolving
  use math
  use rotations
  use mesh
  use material
  use config
  use crystallite
  use homogenization
  use IO
  use discretization
  use DAMASK_interface
  use numerics
  use HDF5_utilities
  use results
  use lattice
  use constitutive

  implicit none
  private
  
  real(pReal),                      parameter,   private :: &
    CPFEM_odd_stress    = 1e15_pReal, &                                                             !< return value for stress in case of ping pong dummy cycle
    CPFEM_odd_jacobian  = 1e50_pReal                                                                !< return value for jacobian in case of ping pong dummy cycle
  real(pReal), dimension (:,:,:),   allocatable, private :: &
    CPFEM_cs                                                                                        !< Cauchy stress
  real(pReal), dimension (:,:,:,:), allocatable, private :: &
    CPFEM_dcsdE                                                                                     !< Cauchy stress tangent
  real(pReal), dimension (:,:,:,:), allocatable, private :: &
    CPFEM_dcsdE_knownGood                                                                           !< known good tangent
  integer(pInt),                                 public :: &
    cycleCounter =  0_pInt, &                                                                       !< needs description
    theInc       = -1_pInt, &                                                                       !< needs description
    lastLovl     =  0_pInt, &                                                                       !< lovl in previous call to marc hypela2
    lastStep     =  0_pInt                                                                          !< kstep in previous call to abaqus umat
  real(pReal),                                   public :: &
    theTime      = 0.0_pReal, &                                                                     !< needs description
    theDelta     = 0.0_pReal    
  logical,                                       public :: & 
    outdatedFFN1      = .false., &                                                                  !< needs description
    lastIncConverged  = .false., &                                                                  !< needs description
    outdatedByNewInc  = .false.                                                                     !< needs description

  logical,                                       public, protected :: &
    CPFEM_init_done       = .false.                                                                 !< remember whether init has been done already
  logical,                                       private :: &
    CPFEM_calc_done       = .false.                                                                 !< remember whether first ip has already calced the results

  integer(pInt), parameter,                      public :: &
    CPFEM_COLLECT         = 2_pInt**0_pInt, &
    CPFEM_CALCRESULTS     = 2_pInt**1_pInt, &
    CPFEM_AGERESULTS      = 2_pInt**2_pInt, &
    CPFEM_BACKUPJACOBIAN  = 2_pInt**3_pInt, &
    CPFEM_RESTOREJACOBIAN = 2_pInt**4_pInt

  public :: &
    CPFEM_general, &
    CPFEM_initAll, &
    CPFEM_results

contains


!--------------------------------------------------------------------------------------------------
!> @brief call (thread safe) all module initializations
!--------------------------------------------------------------------------------------------------
subroutine CPFEM_initAll(el,ip)
 integer(pInt), intent(in) ::                        el, &                                          !< FE el number
                                                     ip                                             !< FE integration point number

 !$OMP CRITICAL (init)
   if (.not. CPFEM_init_done) then
     call DAMASK_interface_init                                                                    ! Spectral and FEM interface to commandline
     call prec_init
     call IO_init
     call numerics_init
     call debug_init
     call config_init
     call math_init
     call rotations_init
     call FE_init
#ifdef DAMASK_HDF5
     call HDF5_utilities_init
     call results_init
#endif
     call mesh_init(ip, el)
     call lattice_init
     call material_init
     call constitutive_init
     call crystallite_init
     call homogenization_init
     call CPFEM_init
     CPFEM_init_done = .true.
   endif
 !$OMP END CRITICAL (init)

end subroutine CPFEM_initAll


!--------------------------------------------------------------------------------------------------
!> @brief allocate the arrays defined in module CPFEM and initialize them
!--------------------------------------------------------------------------------------------------
subroutine CPFEM_init

  write(6,'(/,a)')   ' <<<+-  CPFEM init  -+>>>'
  flush(6)

  allocate(CPFEM_cs(               6,discretization_nIP,discretization_nElem), source= 0.0_pReal)
  allocate(CPFEM_dcsdE(          6,6,discretization_nIP,discretization_nElem), source= 0.0_pReal)
  allocate(CPFEM_dcsdE_knownGood(6,6,discretization_nIP,discretization_nElem), source= 0.0_pReal)

  if (iand(debug_level(debug_CPFEM), debug_levelBasic) /= 0) then
    write(6,'(a32,1x,6(i8,1x))')   'CPFEM_cs:              ', shape(CPFEM_cs)
    write(6,'(a32,1x,6(i8,1x))')   'CPFEM_dcsdE:           ', shape(CPFEM_dcsdE)
    write(6,'(a32,1x,6(i8,1x),/)') 'CPFEM_dcsdE_knownGood: ', shape(CPFEM_dcsdE_knownGood)
    write(6,'(a32,l1)')            'symmetricSolver:       ', symmetricSolver
    flush(6)
  endif

end subroutine CPFEM_init


!--------------------------------------------------------------------------------------------------
!> @brief perform initialization at first call, update variables and call the actual material model
!--------------------------------------------------------------------------------------------------
subroutine CPFEM_general(mode, parallelExecution, ffn, ffn1, temperature_inp, dt, elFE, ip, cauchyStress, jacobian)

 integer(pInt), intent(in) ::                        elFE, &                                        !< FE element number
                                                     ip                                             !< integration point number
 real(pReal), intent(in) ::                          dt                                             !< time increment
 real(pReal), dimension (3,3), intent(in) ::         ffn, &                                         !< deformation gradient for t=t0
                                                     ffn1                                           !< deformation gradient for t=t1
 integer(pInt), intent(in) ::                        mode                                           !< computation mode  1: regular computation plus aging of results
 real(pReal), intent(in) ::                          temperature_inp                                !< temperature
 logical, intent(in) ::                              parallelExecution                              !< flag indicating parallel computation of requested IPs
 real(pReal), dimension(6), intent(out) ::           cauchyStress                                   !< stress as 6 vector
 real(pReal), dimension(6,6), intent(out) ::         jacobian                                       !< jacobian as 66 tensor (Consistent tangent dcs/dE)

 real(pReal)                                         J_inverse, &                                   ! inverse of Jacobian
                                                     rnd
 real(pReal), dimension (3,3) ::                     Kirchhoff, &                                   ! Piola-Kirchhoff stress in Matrix notation
                                                     cauchyStress33                                 ! stress vector in Matrix notation
 real(pReal), dimension (3,3,3,3) ::                 H_sym, &
                                                     H, &
                                                     jacobian3333                                   ! jacobian in Matrix notation

 integer(pInt)                                       elCP, &                                        ! crystal plasticity element number
                                                     i, j, k, l, m, n, ph, homog, mySource
 logical                                             updateJaco                                     ! flag indicating if JAcobian has to be updated

 elCP = mesh_FEasCP('elem',elFE)

 if (iand(debug_level(debug_CPFEM), debug_levelBasic) /= 0_pInt &
     .and. elCP == debug_e .and. ip == debug_i) then
   write(6,'(/,a)') '#############################################'
   write(6,'(a1,a22,1x,i8,a13)')   '#','element',        elCP,         '#'
   write(6,'(a1,a22,1x,i8,a13)')   '#','ip',             ip,           '#'
   write(6,'(a1,a22,1x,f15.7,a6)') '#','theTime',        theTime,      '#'
   write(6,'(a1,a22,1x,f15.7,a6)') '#','theDelta',       theDelta,     '#'
   write(6,'(a1,a22,1x,i8,a13)')   '#','theInc',         theInc,       '#'
   write(6,'(a1,a22,1x,i8,a13)')   '#','cycleCounter',   cycleCounter, '#'
   write(6,'(a1,a22,1x,i8,a13)')   '#','computationMode',mode,         '#'
   if (terminallyIll) &
   write(6,'(a,/)') '#           --- terminallyIll ---           #'
   write(6,'(a,/)') '#############################################'; flush (6)
 endif

 if (iand(mode, CPFEM_BACKUPJACOBIAN) /= 0_pInt) &
   CPFEM_dcsde_knownGood = CPFEM_dcsde
 if (iand(mode, CPFEM_RESTOREJACOBIAN) /= 0_pInt) &
   CPFEM_dcsde = CPFEM_dcsde_knownGood

 !*** age results and write restart data if requested
 if (iand(mode, CPFEM_AGERESULTS) /= 0_pInt) then
   crystallite_F0  = crystallite_partionedF                                                         ! crystallite deformation (_subF is perturbed...)
   crystallite_Fp0 = crystallite_Fp                                                                 ! crystallite plastic deformation
   crystallite_Lp0 = crystallite_Lp                                                                 ! crystallite plastic velocity
   crystallite_Fi0 = crystallite_Fi                                                                 ! crystallite intermediate deformation
   crystallite_Li0 = crystallite_Li                                                                 ! crystallite intermediate velocity
   crystallite_S0  = crystallite_S                                                                  ! crystallite 2nd Piola Kirchhoff stress

   forall ( i = 1:size(plasticState    )) plasticState(i)%state0     = plasticState(i)%state   ! copy state in this lenghty way because: A component cannot be an array if the encompassing structure is an array
   do i = 1, size(sourceState)
     do mySource = 1,phase_Nsources(i)
       sourceState(i)%p(mySource)%state0 = sourceState(i)%p(mySource)%state                    ! copy state in this lenghty way because: A component cannot be an array if the encompassing structure is an array
   enddo; enddo
   if (iand(debug_level(debug_CPFEM), debug_levelBasic) /= 0_pInt) then
     write(6,'(a)') '<< CPFEM >> aging states'
     if (debug_e <= discretization_nElem .and. debug_i <=discretization_nIP) then
       write(6,'(a,1x,i8,1x,i2,1x,i4,/,(12x,6(e20.8,1x)),/)') &
             '<< CPFEM >> aged state of elFE ip grain',debug_e, debug_i, 1, &
              plasticState(material_phaseAt(1,debug_e))%state(:,material_phasememberAt(1,debug_i,debug_e))
       endif
   endif

   do homog = 1_pInt, material_Nhomogenization
     homogState       (homog)%state0 =  homogState       (homog)%state
     thermalState     (homog)%state0 =  thermalState     (homog)%state
     damageState      (homog)%state0 =  damageState      (homog)%state
   enddo
 endif



 !*** collection of FEM input with returning of randomize odd stress and jacobian
 !*   If no parallel execution is required, there is no need to collect FEM input

 if (.not. parallelExecution) then
   chosenThermal1: select case (thermal_type(material_homogenizationAt(elCP)))
     case (THERMAL_conduction_ID) chosenThermal1
       temperature(material_homogenizationAt(elCP))%p(thermalMapping(material_homogenizationAt(elCP))%p(ip,elCP)) = &
         temperature_inp
     end select chosenThermal1
   materialpoint_F0(1:3,1:3,ip,elCP) = ffn
   materialpoint_F(1:3,1:3,ip,elCP) = ffn1

 elseif (iand(mode, CPFEM_COLLECT) /= 0_pInt) then
   call random_number(rnd)
   if (rnd < 0.5_pReal) rnd = rnd - 1.0_pReal
   CPFEM_cs(1:6,ip,elCP) = rnd * CPFEM_odd_stress
   CPFEM_dcsde(1:6,1:6,ip,elCP) = CPFEM_odd_jacobian * math_identity2nd(6)
   chosenThermal2: select case (thermal_type(material_homogenizationAt(elCP)))
     case (THERMAL_conduction_ID) chosenThermal2
       temperature(material_homogenizationAt(elCP))%p(thermalMapping(material_homogenizationAt(elCP))%p(ip,elCP)) = &
         temperature_inp
     end select chosenThermal2
   materialpoint_F0(1:3,1:3,ip,elCP) = ffn
   materialpoint_F(1:3,1:3,ip,elCP) = ffn1
   CPFEM_calc_done = .false.
 endif                                                              ! collection



 !*** calculation of stress and jacobian

 if (iand(mode, CPFEM_CALCRESULTS) /= 0_pInt) then

   !*** deformation gradient outdated or any actual deformation gradient differs more than relevantStrain from the stored one
   validCalculation: if (terminallyIll &
                    .or. outdatedFFN1 &
                    .or. any(abs(ffn1 - materialpoint_F(1:3,1:3,ip,elCP)) > defgradTolerance)) then
     if (any(abs(ffn1 - materialpoint_F(1:3,1:3,ip,elCP)) > defgradTolerance)) then
       if (iand(debug_level(debug_CPFEM), debug_levelBasic) /=  0_pInt) then
           write(6,'(a,1x,i8,1x,i2)') '<< CPFEM >> OUTDATED at elFE ip',elFE,ip
           write(6,'(a,/,3(12x,3(f10.6,1x),/))') '<< CPFEM >> FFN1 old:',&
                                             transpose(materialpoint_F(1:3,1:3,ip,elCP))
           write(6,'(a,/,3(12x,3(f10.6,1x),/))') '<< CPFEM >> FFN1 now:',transpose(ffn1)
       endif
       outdatedFFN1 = .true.
     endif
     call random_number(rnd)
     if (rnd < 0.5_pReal) rnd = rnd - 1.0_pReal
     CPFEM_cs(1:6,ip,elCP) = rnd*CPFEM_odd_stress
     CPFEM_dcsde(1:6,1:6,ip,elCP) = CPFEM_odd_jacobian*math_identity2nd(6)

   !*** deformation gradient is not outdated

   else validCalculation
     updateJaco = mod(cycleCounter,iJacoStiffness) == 0
     !* no parallel computation, so we use just one single elFE and ip for computation

     if (.not. parallelExecution) then
       FEsolving_execElem(1)     = elCP
       FEsolving_execElem(2)     = elCP
       FEsolving_execIP(1,elCP) = ip
       FEsolving_execIP(2,elCP) = ip
       if (iand(debug_level(debug_CPFEM), debug_levelExtensive) /=  0_pInt) &
         write(6,'(a,i8,1x,i2)') '<< CPFEM >> calculation for elFE ip ',elFE,ip
       call materialpoint_stressAndItsTangent(updateJaco, dt)                                     ! calculate stress and its tangent
       call materialpoint_postResults()

     !* parallel computation and calulation not yet done

     elseif (.not. CPFEM_calc_done) then
       if (iand(debug_level(debug_CPFEM), debug_levelExtensive) /= 0_pInt) &
         write(6,'(a,i8,a,i8)') '<< CPFEM >> calculation for elements ',FEsolving_execElem(1),&
                                                                 ' to ',FEsolving_execElem(2)
       call materialpoint_stressAndItsTangent(updateJaco, dt)                                       ! calculate stress and its tangent (parallel execution inside)
       call materialpoint_postResults()
       CPFEM_calc_done = .true.
     endif

     !* map stress and stiffness (or return odd values if terminally ill)
     terminalIllness: if ( terminallyIll ) then

       call random_number(rnd)
       if (rnd < 0.5_pReal) rnd = rnd - 1.0_pReal
       CPFEM_cs(1:6,ip,elCP) = rnd * CPFEM_odd_stress
       CPFEM_dcsde(1:6,1:6,ip,elCP) = CPFEM_odd_jacobian * math_identity2nd(6)

     else terminalIllness


       ! translate from P to CS
       Kirchhoff = matmul(materialpoint_P(1:3,1:3,ip,elCP), transpose(materialpoint_F(1:3,1:3,ip,elCP)))
       J_inverse  = 1.0_pReal / math_det33(materialpoint_F(1:3,1:3,ip,elCP))
       CPFEM_cs(1:6,ip,elCP) = math_sym33to6(J_inverse * Kirchhoff,weighted=.false.)

       !  translate from dP/dF to dCS/dE
       H = 0.0_pReal
       do i=1,3; do j=1,3; do k=1,3; do l=1,3; do m=1,3; do n=1,3
         H(i,j,k,l) = H(i,j,k,l) &
                    +  materialpoint_F(j,m,ip,elCP) * materialpoint_F(l,n,ip,elCP) &
                                                    * materialpoint_dPdF(i,m,k,n,ip,elCP) &
                    -  math_delta(j,l) * materialpoint_F(i,m,ip,elCP) * materialpoint_P(k,m,ip,elCP) &
                    +  0.5_pReal * (  Kirchhoff(j,l)*math_delta(i,k) + Kirchhoff(i,k)*math_delta(j,l) &
                                    + Kirchhoff(j,k)*math_delta(i,l) + Kirchhoff(i,l)*math_delta(j,k))
       enddo; enddo; enddo; enddo; enddo; enddo

       forall(i=1:3, j=1:3,k=1:3,l=1:3) &
         H_sym(i,j,k,l) = 0.25_pReal * (H(i,j,k,l) + H(j,i,k,l) + H(i,j,l,k) + H(j,i,l,k))

       CPFEM_dcsde(1:6,1:6,ip,elCP) = math_sym3333to66(J_inverse * H_sym,weighted=.false.)

     endif terminalIllness
   endif validCalculation

   !* report stress and stiffness
   if ((iand(debug_level(debug_CPFEM), debug_levelExtensive) /= 0_pInt) &
        .and. ((debug_e == elCP .and. debug_i == ip) &
               .or. .not. iand(debug_level(debug_CPFEM), debug_levelSelective) /= 0_pInt)) then
       write(6,'(a,i8,1x,i2,/,12x,6(f10.3,1x)/)') &
         '<< CPFEM >> stress/MPa at elFE ip ',   elFE, ip, CPFEM_cs(1:6,ip,elCP)*1.0e-6_pReal
       write(6,'(a,i8,1x,i2,/,6(12x,6(f10.3,1x)/))') &
         '<< CPFEM >> Jacobian/GPa at elFE ip ', elFE, ip, transpose(CPFEM_dcsdE(1:6,1:6,ip,elCP))*1.0e-9_pReal
       flush(6)
   endif

 endif

 !*** warn if stiffness close to zero
 if (all(abs(CPFEM_dcsdE(1:6,1:6,ip,elCP)) < 1e-10_pReal)) call IO_warning(601,elCP,ip)

 !*** copy to output if using commercial FEM solver
 cauchyStress = CPFEM_cs   (1:6,    ip,elCP)
 jacobian     = CPFEM_dcsdE(1:6,1:6,ip,elCP)


 !*** remember extreme values of stress ...
 cauchyStress33 = math_6toSym33(CPFEM_cs(1:6,ip,elCP),weighted=.false.)
 if (maxval(cauchyStress33) > debug_stressMax) then
   debug_stressMaxLocation = [elCP, ip]
   debug_stressMax = maxval(cauchyStress33)
 endif
 if (minval(cauchyStress33) < debug_stressMin) then
   debug_stressMinLocation = [elCP, ip]
   debug_stressMin = minval(cauchyStress33)
 endif
 !*** ... and Jacobian
 jacobian3333 = math_66toSym3333(CPFEM_dcsdE(1:6,1:6,ip,elCP),weighted=.false.)
 if (maxval(jacobian3333) > debug_jacobianMax) then
   debug_jacobianMaxLocation = [elCP, ip]
   debug_jacobianMax = maxval(jacobian3333)
 endif
 if (minval(jacobian3333) < debug_jacobianMin) then
   debug_jacobianMinLocation = [elCP, ip]
   debug_jacobianMin = minval(jacobian3333)
 endif

end subroutine CPFEM_general


!--------------------------------------------------------------------------------------------------
!> @brief triggers writing of the results
!--------------------------------------------------------------------------------------------------
subroutine CPFEM_results(inc,time)

  integer(pInt), intent(in) :: inc
  real(pReal),   intent(in) :: time

#ifdef DAMASK_HDF5
  call results_openJobFile
  call results_addIncrement(inc,time)
  call constitutive_results
  call crystallite_results
  call homogenization_results
  call results_removeLink('current') ! ToDo: put this into closeJobFile
  call results_closeJobFile
#endif

end subroutine CPFEM_results

end module CPFEM
