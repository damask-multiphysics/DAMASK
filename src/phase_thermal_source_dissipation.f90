!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine for thermal source due to plastic dissipation
!> @details to be done
!--------------------------------------------------------------------------------------------------
submodule(phase:thermal) source_dissipation

  type :: tParameters                                                                               !< container type for internal constitutive parameters
    real(pREAL) :: &
      kappa                                                                                         !< TAYLOR-QUINNEY factor
  end type tParameters

  type(tParameters), dimension(:),   allocatable :: param                                           !< containers of constitutive parameters (len Ninstances)


contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
module function source_dissipation_init(maxNsources) result(isMySource)

  integer, intent(in)                  :: maxNsources
  logical, dimension(:,:), allocatable :: isMySource

  type(tDict), pointer :: &
    phases, &
    phase, &
    thermal, &
    src
  class(tList), pointer :: &
    sources
  character(len=:), allocatable :: refs
  integer :: ph,Nmembers,so,Nsources


  isMySource = thermal_active('dissipation',maxNsources)
  if (count(isMySource) == 0) return

  print'(/,1x,a)', '<<<+-  phase:thermal:source_dissipation init  -+>>>'
  print'(/,1x,a,1x,i0)', '# phases:',count(isMySource); flush(IO_STDOUT)


  phases => config_material%get_dict('phase')
  allocate(param(phases%length))

  do ph = 1, phases%length
    Nsources = count(isMySource(:,ph))
    if (Nsources == 0) cycle
    if (Nsources > 1) call IO_error(600,ext_msg='dissipation')
    Nmembers = count(material_ID_phase == ph)
    phase => phases%get_dict(ph)
    thermal => phase%get_dict('thermal')
    sources => thermal%get_list('source')
    do so = 1, sources%length
      if (isMySource(so,ph)) then
        associate(prm  => param(ph))
          src => sources%get_dict(so)
          print'(/,1x,a,1x,i0,1x,a,1x,a,1x,i0)', 'phase',ph,'('//phases%key(ph)//')','source',so
          refs = config_listReferences(src,indent=3)
          if (len(refs) > 0) print'(/,1x,a)', refs

          prm%kappa = src%get_asReal('kappa')
          call phase_allocateState(thermalState(ph)%p(so),Nmembers,0,0,0)
        end associate
        exit
      end if
    end do
  end do


end function source_dissipation_init


!--------------------------------------------------------------------------------------------------
!> @brief Ninstancess dissipation rate
!--------------------------------------------------------------------------------------------------
module function source_dissipation_f_T(ph,en) result(f_T)

  integer, intent(in) :: ph, en
  real(pREAL) :: &
    f_T
  real(pREAL), dimension(3,3) :: &
    Mp                                                                                              !< Mandel stress work conjugate with Lp

  Mp = matmul(matmul(transpose(mechanical_F_i(ph,en)),mechanical_F_i(ph,en)),mechanical_S(ph,en))

  associate(prm => param(ph))
    f_T = prm%kappa*sum(abs(Mp*mechanical_L_p(ph,en)))
  end associate

end function source_dissipation_f_T

end submodule source_dissipation
