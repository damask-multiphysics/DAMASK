
!************************************
!*      Module: CONSTITUTIVE        *
!************************************
!* contains:                        *
!* - constitutive equations         *
!* - parameters definition          *
!************************************

MODULE constitutive

!*** Include other modules ***
 use prec
 implicit none

 type(p_vec),   dimension(:,:,:), allocatable :: constitutive_state_old, &     ! pointer array to old state variables of each grain
                                                 constitutive_state_new        ! pointer array to new state variables of each grain
 integer(pInt), dimension(:,:,:), allocatable :: constitutive_sizeDotState, &  ! size of dotState array
                                                 constitutive_sizeState, &     ! size of state array per grain
                                                 constitutive_sizePostResults  ! size of postResults array per grain
 integer(pInt) constitutive_maxSizeDotState,constitutive_maxSizeState,constitutive_maxSizePostResults

CONTAINS
!****************************************
!* - constitutive_init
!* - constitutive_homogenizedC
!* - constitutive_microstructure
!* - constitutive_LpAndItsTangent
!* - constitutive_dotState
!* - constitutive_postResults
!****************************************


subroutine constitutive_init()
!**************************************
!*      Module initialization         *
!**************************************
 use prec, only: pReal,pInt
 use IO, only: IO_error, IO_open_file
 use mesh, only: mesh_maxNips,mesh_NcpElems,mesh_element,FE_Nips
 use material
 use constitutive_phenomenological
 use constitutive_j2
 use constitutive_dislobased

 integer(pInt), parameter :: fileunit = 200
 integer(pInt) e,i,g,myInstance

 if(.not. IO_open_file(fileunit,material_configFile)) call IO_error (100) ! corrupt config file

 call constitutive_phenomenological_init(fileunit)       ! parse all phases of this constitution
 call constitutive_j2_init(fileunit)
 call constitutive_dislobased_init(fileunit)

 close(fileunit)

 allocate(constitutive_state_old(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems))
 allocate(constitutive_state_new(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems))
 allocate(constitutive_sizeDotState(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems)) ;    constitutive_sizeDotState = 0_pInt
 allocate(constitutive_sizeState(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems)) ;       constitutive_sizeState = 0_pInt
 allocate(constitutive_sizePostResults(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems)) ; constitutive_sizePostResults = 0_pInt

 do e = 1,mesh_NcpElems                                  ! loop over elements
   do i = 1,FE_Nips(mesh_element(2,e))                   ! loop over IPs
     do g = 1,homogenization_Ngrains(mesh_element(3,e))  ! loop over grains
       myInstance = phase_constitutionInstance(material_phase(g,i,e))
       select case(phase_constitution(material_phase(g,i,e)))
         case (constitutive_phenomenological_label)
           allocate(constitutive_state_old(g,i,e)%p(constitutive_phenomenological_sizeState(myInstance)))
           allocate(constitutive_state_new(g,i,e)%p(constitutive_phenomenological_sizeState(myInstance)))
           constitutive_state_new(g,i,e)%p =        constitutive_phenomenological_stateInit(myInstance)
           constitutive_state_old(g,i,e)%p =        constitutive_phenomenological_stateInit(myInstance)
           constitutive_sizeDotState(g,i,e) =       constitutive_phenomenological_sizeDotState(myInstance)
           constitutive_sizeState(g,i,e) =          constitutive_phenomenological_sizeState(myInstance)
           constitutive_sizePostResults(g,i,e) =    constitutive_phenomenological_sizePostResults(myInstance)
         case (constitutive_j2_label)
           allocate(constitutive_state_old(g,i,e)%p(constitutive_j2_sizeState(myInstance)))
           allocate(constitutive_state_new(g,i,e)%p(constitutive_j2_sizeState(myInstance)))
           constitutive_state_new(g,i,e)%p =        constitutive_j2_stateInit(myInstance)
           constitutive_state_old(g,i,e)%p =        constitutive_j2_stateInit(myInstance)
           constitutive_sizeDotState(g,i,e) =       constitutive_j2_sizeDotState(myInstance)
           constitutive_sizeState(g,i,e) =          constitutive_j2_sizeState(myInstance)
           constitutive_sizePostResults(g,i,e) =    constitutive_j2_sizePostResults(myInstance)
         case (constitutive_dislobased_label)
           allocate(constitutive_state_old(g,i,e)%p(constitutive_dislobased_sizeState(myInstance)))
           allocate(constitutive_state_new(g,i,e)%p(constitutive_dislobased_sizeState(myInstance)))
           constitutive_state_new(g,i,e)%p =        constitutive_dislobased_stateInit(myInstance)
           constitutive_state_old(g,i,e)%p =        constitutive_dislobased_stateInit(myInstance)
           constitutive_sizeDotState(g,i,e) =       constitutive_dislobased_sizeDotState(myInstance)
           constitutive_sizeState(g,i,e) =          constitutive_dislobased_sizeState(myInstance)
           constitutive_sizePostResults(g,i,e) =    constitutive_dislobased_sizePostResults(myInstance)
         case default
           call IO_error(200,material_phase(g,i,e))      ! unknown constitution
       end select
     enddo
   enddo
 enddo
 constitutive_maxSizeDotState    = maxval(constitutive_sizeDotState)
 constitutive_maxSizeState       = maxval(constitutive_sizeState)
 constitutive_maxSizePostResults = maxval(constitutive_sizePostResults)

 return

end subroutine


function constitutive_homogenizedC(ipc,ip,el)
!*********************************************************************
!* This function returns the homogenized elacticity matrix           *
!* INPUT:                                                            *
!*  - state           : state variables                              *
!*  - ipc             : component-ID of current integration point    *
!*  - ip              : current integration point                    *
!*  - el              : current element                              *
!*********************************************************************
 use prec, only: pReal,pInt
 use material, only: phase_constitution,material_phase
 use constitutive_phenomenological
 use constitutive_j2
 use constitutive_dislobased
 implicit none

 !* Definition of variables
 integer(pInt) ipc,ip,el
 real(pReal), dimension(6,6) :: constitutive_homogenizedC

 select case (phase_constitution(material_phase(ipc,ip,el)))
   case (constitutive_phenomenological_label)
     constitutive_homogenizedC = constitutive_phenomenological_homogenizedC(constitutive_state_new,ipc,ip,el)
   case (constitutive_j2_label)
     constitutive_homogenizedC = constitutive_j2_homogenizedC(constitutive_state_new,ipc,ip,el)
   case (constitutive_dislobased_label)
     constitutive_homogenizedC = constitutive_dislobased_homogenizedC(constitutive_state_new,ipc,ip,el)

 end select

return
end function


subroutine constitutive_microstructure(Temperature,ipc,ip,el)
!*********************************************************************
!* This function calculates from state needed variables              *
!* INPUT:                                                            *
!*  - state           : state variables                              *
!*  - Tp              : temperature                                  *
!*  - ipc             : component-ID of current integration point    *
!*  - ip              : current integration point                    *
!*  - el              : current element                              *
!*********************************************************************
 use prec, only: pReal,pInt
 use material, only: phase_constitution,material_phase
 use constitutive_phenomenological
 use constitutive_j2
 use constitutive_dislobased
 implicit none

!* Definition of variables
integer(pInt) ipc,ip,el
real(pReal) Temperature

 select case (phase_constitution(material_phase(ipc,ip,el)))
   case (constitutive_phenomenological_label)
     call constitutive_phenomenological_microstructure(Temperature,constitutive_state_new,ipc,ip,el)
   case (constitutive_j2_label)
     call constitutive_j2_microstructure(Temperature,constitutive_state_new,ipc,ip,el)
   case (constitutive_dislobased_label)
     call constitutive_dislobased_microstructure(Temperature,constitutive_state_new,ipc,ip,el)

 end select

end subroutine


subroutine constitutive_LpAndItsTangent(Lp,dLp_dTstar, Tstar_v,Temperature,ipc,ip,el)
!*********************************************************************
!* This subroutine contains the constitutive equation for            *
!* calculating the velocity gradient                                 *
!* INPUT:                                                            *
!*  - Tstar_v         : 2nd Piola Kirchhoff stress tensor (Mandel)   *
!*  - ipc             : component-ID of current integration point    *
!*  - ip              : current integration point                    *
!*  - el              : current element                              *
!* OUTPUT:                                                           *
!*  - Lp              : plastic velocity gradient                    *
!*  - dLp_dTstar      : derivative of Lp (4th-order tensor)          *
!*********************************************************************
 use prec, only: pReal,pInt
 use material, only: phase_constitution,material_phase
 use constitutive_phenomenological
 use constitutive_j2
 use constitutive_dislobased
 implicit none

!* Definition of variables
 integer(pInt) ipc,ip,el
 real(pReal) Temperature
 real(pReal), dimension(6) :: Tstar_v
 real(pReal), dimension(3,3) :: Lp
 real(pReal), dimension(9,9) :: dLp_dTstar

 select case (phase_constitution(material_phase(ipc,ip,el)))
   case (constitutive_phenomenological_label)
     call constitutive_phenomenological_LpAndItsTangent(Lp,dLp_dTstar,Tstar_v,Temperature,constitutive_state_new,ipc,ip,el)
   case (constitutive_j2_label)
     call constitutive_j2_LpAndItsTangent(Lp,dLp_dTstar,Tstar_v,Temperature,constitutive_state_new,ipc,ip,el)
   case (constitutive_dislobased_label)
     call constitutive_dislobased_LpAndItsTangent(Lp,dLp_dTstar,Tstar_v,Temperature,constitutive_state_new,ipc,ip,el)

 end select

 return
end subroutine


function constitutive_dotState(Tstar_v,Temperature,ipc,ip,el)
!*********************************************************************
!* This subroutine contains the constitutive equation for            *
!* calculating the rate of change of microstructure                  *
!* INPUT:                                                            *
!*  - Tstar_v         : 2nd Piola Kirchhoff stress tensor (Mandel)   *
!*  - state           : current microstructure                       *
!*  - ipc             : component-ID of current integration point    *
!*  - ip              : current integration point                    *
!*  - el              : current element                              *
!* OUTPUT:                                                           *
!*  - constitutive_dotState : evolution of state variable            *
!*********************************************************************
 use prec, only: pReal,pInt
 use material, only: phase_constitution,material_phase
 use constitutive_phenomenological
 use constitutive_j2
 use constitutive_dislobased
 implicit none

!* Definition of variables
 integer(pInt) ipc,ip,el
 real(pReal) Temperature
 real(pReal), dimension(6) :: Tstar_v
 real(pReal), dimension(constitutive_sizeDotState(ipc,ip,el)) :: constitutive_dotState

 select case (phase_constitution(material_phase(ipc,ip,el)))
   case (constitutive_phenomenological_label)
     constitutive_dotState = constitutive_phenomenological_dotState(Tstar_v,Temperature,constitutive_state_new,ipc,ip,el)
   case (constitutive_j2_label)
     constitutive_dotState = constitutive_j2_dotState(Tstar_v,Temperature,constitutive_state_new,ipc,ip,el)
   case (constitutive_dislobased_label)
     constitutive_dotState = constitutive_dislobased_dotState(Tstar_v,Temperature,constitutive_state_new,ipc,ip,el)

 end select
 return
end function


pure function constitutive_postResults(Tstar_v,Temperature,dt,ipc,ip,el)
!*********************************************************************
!* return array of constitutive results                              *
!* INPUT:                                                            *
!*  - Tstar_v         : 2nd Piola Kirchhoff stress tensor (Mandel)   *
!*  - dt              : current time increment                       *
!*  - ipc             : component-ID of current integration point    *
!*  - ip              : current integration point                    *
!*  - el              : current element                              *
!*********************************************************************
 use prec, only: pReal,pInt
 use material, only: phase_constitution,material_phase
 use constitutive_phenomenological
 use constitutive_j2
 use constitutive_dislobased
 implicit none

!* Definition of variables
 integer(pInt), intent(in) :: ipc,ip,el
 real(pReal), intent(in) :: dt,Temperature
 real(pReal), dimension(6), intent(in) :: Tstar_v
 real(pReal), dimension(constitutive_sizePostResults(ipc,ip,el)) :: constitutive_postResults

 constitutive_postResults = 0.0_pReal
 select case (phase_constitution(material_phase(ipc,ip,el)))
   case (constitutive_phenomenological_label)
     constitutive_postResults = constitutive_phenomenological_postResults(Tstar_v,Temperature,dt,constitutive_state_new,ipc,ip,el)
   case (constitutive_j2_label)
     constitutive_postResults = constitutive_j2_postResults(Tstar_v,Temperature,dt,constitutive_state_new,ipc,ip,el)
   case (constitutive_dislobased_label)
     constitutive_postResults = constitutive_dislobased_postResults(Tstar_v,Temperature,dt,constitutive_state_new,ipc,ip,el)

 end select

return

end function

END MODULE