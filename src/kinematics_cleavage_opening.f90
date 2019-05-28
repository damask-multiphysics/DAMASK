!--------------------------------------------------------------------------------------------------
!> @author Luv Sharma, Max-Planck-Institut für Eisenforschung GmbH
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine incorporating kinematics resulting from opening of cleavage planes
!> @details to be done
!--------------------------------------------------------------------------------------------------
module kinematics_cleavage_opening
  use prec
  use IO
  use config
  use debug
  use math
  use lattice
  use material

  implicit none
  private

  integer, dimension(:), allocatable :: kinematics_cleavage_opening_instance

  type :: tParameters                                                                               !< container type for internal constitutive parameters
    integer :: &
      totalNcleavage
    integer, dimension(:),   allocatable :: &
      Ncleavage                                                                                     !< active number of cleavage systems per family
    real(pReal) :: &
      sdot0, &
      n
    real(pReal),   dimension(:),   allocatable :: &
      critDisp, &
      critLoad
  end type

! Begin Deprecated
  integer,                       dimension(:),           allocatable :: &
    kinematics_cleavage_opening_totalNcleavage                                                      !< total number of cleavage systems
    
  integer,                       dimension(:,:),         allocatable :: &
    kinematics_cleavage_opening_Ncleavage                                                           !< number of cleavage systems per family
    
  real(pReal),                   dimension(:),           allocatable :: &
    kinematics_cleavage_opening_sdot_0, &
    kinematics_cleavage_opening_N

  real(pReal),                   dimension(:,:),         allocatable :: &
    kinematics_cleavage_opening_critDisp, &
    kinematics_cleavage_opening_critLoad
! End Deprecated

  public :: &
    kinematics_cleavage_opening_init, &
    kinematics_cleavage_opening_LiAndItsTangent

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine kinematics_cleavage_opening_init

  integer, allocatable, dimension(:) :: tempInt
  real(pReal), allocatable, dimension(:) :: tempFloat

  integer :: maxNinstance,p,instance

  write(6,'(/,a)')   ' <<<+-  kinematics_'//KINEMATICS_cleavage_opening_LABEL//' init  -+>>>'

  maxNinstance = count(phase_kinematics == KINEMATICS_cleavage_opening_ID)
  if (maxNinstance == 0) return
  
  if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0) &
    write(6,'(a16,1x,i5,/)') '# instances:',maxNinstance
  
  allocate(kinematics_cleavage_opening_instance(size(config_phase)), source=0)
  do p = 1, size(config_phase)
    kinematics_cleavage_opening_instance(p) = count(phase_kinematics(:,1:p) == kinematics_cleavage_opening_ID) ! ToDo: count correct?
  enddo
    
  allocate(kinematics_cleavage_opening_critDisp(lattice_maxNcleavageFamily,maxNinstance),  source=0.0_pReal) 
  allocate(kinematics_cleavage_opening_critLoad(lattice_maxNcleavageFamily,maxNinstance),  source=0.0_pReal) 
  allocate(kinematics_cleavage_opening_Ncleavage(lattice_maxNcleavageFamily,maxNinstance), source=0)
  allocate(kinematics_cleavage_opening_totalNcleavage(maxNinstance),                       source=0)
  allocate(kinematics_cleavage_opening_sdot_0(maxNinstance),                               source=0.0_pReal) 
  allocate(kinematics_cleavage_opening_N(maxNinstance),                                    source=0.0_pReal) 

  do p = 1, size(config_phase)
    if (all(phase_kinematics(:,p) /= KINEMATICS_cleavage_opening_ID)) cycle
    instance = kinematics_cleavage_opening_instance(p)
    kinematics_cleavage_opening_sdot_0(instance) = config_phase(p)%getFloat('anisobrittle_sdot0')
    kinematics_cleavage_opening_N(instance)      = config_phase(p)%getFloat('anisobrittle_ratesensitivity')
    tempInt = config_phase(p)%getInts('ncleavage')     
    kinematics_cleavage_opening_Ncleavage(1:size(tempInt),instance) = tempInt

    tempFloat = config_phase(p)%getFloats('anisobrittle_criticaldisplacement',requiredSize=size(tempInt))
    kinematics_cleavage_opening_critDisp(1:size(tempInt),instance) = tempFloat

    tempFloat = config_phase(p)%getFloats('anisobrittle_criticalload',requiredSize=size(tempInt))
    kinematics_cleavage_opening_critLoad(1:size(tempInt),instance) = tempFloat

    kinematics_cleavage_opening_Ncleavage(1:lattice_maxNcleavageFamily,instance) = &
        min(lattice_NcleavageSystem(1:lattice_maxNcleavageFamily,p),&                            ! limit active cleavage systems per family to min of available and requested
            kinematics_cleavage_opening_Ncleavage(1:lattice_maxNcleavageFamily,instance))
      kinematics_cleavage_opening_totalNcleavage(instance)  = sum(kinematics_cleavage_opening_Ncleavage(:,instance)) ! how many cleavage systems altogether
      if (kinematics_cleavage_opening_sdot_0(instance) <= 0.0_pReal) &
        call IO_error(211,el=instance,ext_msg='sdot_0 ('//KINEMATICS_cleavage_opening_LABEL//')')
      if (any(kinematics_cleavage_opening_critDisp(1:size(tempInt),instance) < 0.0_pReal)) &
        call IO_error(211,el=instance,ext_msg='critical_displacement ('//KINEMATICS_cleavage_opening_LABEL//')')
      if (any(kinematics_cleavage_opening_critLoad(1:size(tempInt),instance) < 0.0_pReal)) &
        call IO_error(211,el=instance,ext_msg='critical_load ('//KINEMATICS_cleavage_opening_LABEL//')')
      if (kinematics_cleavage_opening_N(instance) <= 0.0_pReal) &
        call IO_error(211,el=instance,ext_msg='rate_sensitivity ('//KINEMATICS_cleavage_opening_LABEL//')')
  enddo
 
end subroutine kinematics_cleavage_opening_init
 
!--------------------------------------------------------------------------------------------------
!> @brief  contains the constitutive equation for calculating the velocity gradient  
!--------------------------------------------------------------------------------------------------
subroutine kinematics_cleavage_opening_LiAndItsTangent(Ld, dLd_dTstar, S, ipc, ip, el)
 
  integer, intent(in) :: &
    ipc, &                                                                                          !< grain number
    ip, &                                                                                           !< integration point number
    el                                                                                              !< element number
  real(pReal),   intent(in),  dimension(3,3) :: &
    S
  real(pReal),   intent(out), dimension(3,3) :: &
    Ld                                                                                              !< damage velocity gradient
  real(pReal),   intent(out), dimension(3,3,3,3) :: &
    dLd_dTstar                                                                                      !< derivative of Ld with respect to Tstar (4th-order tensor)
  integer :: &
    instance, phase, &
    homog, damageOffset, &
    f, i, index_myFamily, k, l, m, n
  real(pReal) :: &
    traction_d, traction_t, traction_n, traction_crit, &
    udotd, dudotd_dt, udott, dudott_dt, udotn, dudotn_dt

  phase = material_phase(ipc,ip,el)
  instance = kinematics_cleavage_opening_instance(phase)
  homog = material_homogenizationAt(el)
  damageOffset = damageMapping(homog)%p(ip,el)
  
  Ld = 0.0_pReal
  dLd_dTstar = 0.0_pReal
  do f = 1,lattice_maxNcleavageFamily
    index_myFamily = sum(lattice_NcleavageSystem(1:f-1,phase))                                      ! at which index starts my family
    do i = 1,kinematics_cleavage_opening_Ncleavage(f,instance)                                      ! process each (active) cleavage system in family
      traction_d    = math_mul33xx33(S,lattice_Scleavage(1:3,1:3,1,index_myFamily+i,phase))
      traction_t    = math_mul33xx33(S,lattice_Scleavage(1:3,1:3,2,index_myFamily+i,phase))
      traction_n    = math_mul33xx33(S,lattice_Scleavage(1:3,1:3,3,index_myFamily+i,phase))
      traction_crit = kinematics_cleavage_opening_critLoad(f,instance)* &
                      damage(homog)%p(damageOffset)*damage(homog)%p(damageOffset)
      udotd = &
        sign(1.0_pReal,traction_d)* &
        kinematics_cleavage_opening_sdot_0(instance)* &
        (max(0.0_pReal, abs(traction_d) - traction_crit)/traction_crit)**kinematics_cleavage_opening_N(instance)
      if (abs(udotd) > tol_math_check) then
        Ld = Ld + udotd*lattice_Scleavage(1:3,1:3,1,index_myFamily+i,phase)
        dudotd_dt = sign(1.0_pReal,traction_d)*udotd*kinematics_cleavage_opening_N(instance)/ &
                    max(0.0_pReal, abs(traction_d) - traction_crit)
        forall (k=1:3,l=1:3,m=1:3,n=1:3) &
          dLd_dTstar(k,l,m,n) = dLd_dTstar(k,l,m,n) + &
            dudotd_dt*lattice_Scleavage(k,l,1,index_myFamily+i,phase)* &
                      lattice_Scleavage(m,n,1,index_myFamily+i,phase)
      endif                

      udott = &
        sign(1.0_pReal,traction_t)* &
        kinematics_cleavage_opening_sdot_0(instance)* &
        (max(0.0_pReal, abs(traction_t) - traction_crit)/traction_crit)**kinematics_cleavage_opening_N(instance)
      if (abs(udott) > tol_math_check) then
        Ld = Ld + udott*lattice_Scleavage(1:3,1:3,2,index_myFamily+i,phase)
        dudott_dt = sign(1.0_pReal,traction_t)*udott*kinematics_cleavage_opening_N(instance)/ &
                    max(0.0_pReal, abs(traction_t) - traction_crit)  
        forall (k=1:3,l=1:3,m=1:3,n=1:3) &
          dLd_dTstar(k,l,m,n) = dLd_dTstar(k,l,m,n) + &
            dudott_dt*lattice_Scleavage(k,l,2,index_myFamily+i,phase)* &
                      lattice_Scleavage(m,n,2,index_myFamily+i,phase)
      endif                

      udotn = &
        sign(1.0_pReal,traction_n)* &
        kinematics_cleavage_opening_sdot_0(instance)* &
        (max(0.0_pReal, abs(traction_n) - traction_crit)/traction_crit)**kinematics_cleavage_opening_N(instance)
      if (abs(udotn) > tol_math_check) then
        Ld = Ld + udotn*lattice_Scleavage(1:3,1:3,3,index_myFamily+i,phase)
        dudotn_dt = sign(1.0_pReal,traction_n)*udotn*kinematics_cleavage_opening_N(instance)/ &
                    max(0.0_pReal, abs(traction_n) - traction_crit)
        forall (k=1:3,l=1:3,m=1:3,n=1:3) &
          dLd_dTstar(k,l,m,n) = dLd_dTstar(k,l,m,n) + &
            dudotn_dt*lattice_Scleavage(k,l,3,index_myFamily+i,phase)* &
                      lattice_Scleavage(m,n,3,index_myFamily+i,phase)
      endif                
    enddo
  enddo
 
end subroutine kinematics_cleavage_opening_LiAndItsTangent

end module kinematics_cleavage_opening
