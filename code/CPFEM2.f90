!--------------------------------------------------------------------------------------------------
! $Id: CPFEM.f90 4761 2016-01-17 13:29:42Z MPIE\m.diehl $
!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief needs a good name and description
!--------------------------------------------------------------------------------------------------
module CPFEM2

 implicit none
 private

 public :: &
   CPFEM_general, &
   CPFEM_initAll

contains


!--------------------------------------------------------------------------------------------------
!> @brief call (thread safe) all module initializations
!--------------------------------------------------------------------------------------------------
subroutine CPFEM_initAll(el,ip)
 use prec, only: &
   pInt
 use prec, only: &
   prec_init
 use numerics, only: &
   numerics_init
 use debug, only: &
   debug_init
 use FEsolving, only: &
   FE_init
 use math, only: &
   math_init
 use mesh, only: &
   mesh_init
 use lattice, only: &
   lattice_init
 use material, only: &
   material_init
 use constitutive, only: &
   constitutive_init
 use crystallite, only: &
   crystallite_init
 use homogenization, only: &
   homogenization_init, &
materialpoint_postResults
 use IO, only: &
   IO_init
 use DAMASK_interface
#ifdef FEM
 use FEZoo, only: &
   FEZoo_init
#endif

 implicit none
 integer(pInt), intent(in) ::                        el, &                                          !< FE el number
                                                     ip                                             !< FE integration point number

 call DAMASK_interface_init                                                                    ! Spectral and FEM interface to commandline
 call prec_init
 call IO_init
#ifdef FEM
 call FEZoo_init
#endif
 call numerics_init
 call debug_init
 call math_init
 call FE_init
 call mesh_init(ip, el)                                                                        ! pass on coordinates to alter calcMode of first ip
 call lattice_init
 call material_init
 call constitutive_init
 call crystallite_init
 call homogenization_init
 call materialpoint_postResults

end subroutine CPFEM_initAll


!--------------------------------------------------------------------------------------------------
!> @brief allocate the arrays defined in module CPFEM and initialize them
!--------------------------------------------------------------------------------------------------
subroutine CPFEM_init
 use, intrinsic :: iso_fortran_env                                                                 ! to get compiler_version and compiler_options (at least for gfortran 4.6 at the moment)
 use prec, only: &
   pInt
 use IO, only: &
   IO_read_realFile,&
   IO_read_intFile, &
   IO_timeStamp, &
   IO_error
 use numerics, only: &
   worldrank
 use debug, only: &
   debug_level, &
   debug_CPFEM, &
   debug_levelBasic, &
   debug_levelExtensive
 use FEsolving, only: &
   restartRead, &
   modelName
 use material, only: &
   material_phase, &
   homogState, &
   phase_plasticity, &
   plasticState, &
   material_Nhomogenization
 use crystallite, only: &
   crystallite_F0, &
   crystallite_Fp0, &
   crystallite_Lp0, &
   crystallite_Fi0, &
   crystallite_Li0, &
   crystallite_dPdF0, &
   crystallite_Tstar0_v

 implicit none
 integer(pInt) :: k,l,m,ph,homog
 character(len=1024) :: rankStr

 mainProcess: if (worldrank == 0) then 
   write(6,'(/,a)')   ' <<<+-  CPFEM init  -+>>>'
   write(6,'(a)')     ' $Id: CPFEM.f90 4761 2016-01-17 13:29:42Z MPIE\m.diehl $'
   write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"
 endif mainProcess

 ! *** restore the last converged values of each essential variable from the binary file
 if (restartRead) then
   if (iand(debug_level(debug_CPFEM), debug_levelExtensive) /= 0_pInt) then
     write(6,'(a)') '<< CPFEM >> restored state variables of last converged step from binary files'
     flush(6)
   endif
   
   write(rankStr,'(a1,i0)')'_',worldrank

   call IO_read_intFile(777,'recordedPhase'//trim(rankStr),modelName,size(material_phase))
   read (777,rec=1) material_phase
   close (777)

   call IO_read_realFile(777,'convergedF'//trim(rankStr),modelName,size(crystallite_F0))
   read (777,rec=1) crystallite_F0
   close (777)

   call IO_read_realFile(777,'convergedFp'//trim(rankStr),modelName,size(crystallite_Fp0))
   read (777,rec=1) crystallite_Fp0
   close (777)

   call IO_read_realFile(777,'convergedFi'//trim(rankStr),modelName,size(crystallite_Fi0))
   read (777,rec=1) crystallite_Fi0
   close (777)

   call IO_read_realFile(777,'convergedLp'//trim(rankStr),modelName,size(crystallite_Lp0))
   read (777,rec=1) crystallite_Lp0
   close (777)

   call IO_read_realFile(777,'convergedLi'//trim(rankStr),modelName,size(crystallite_Li0))
   read (777,rec=1) crystallite_Li0
   close (777)

   call IO_read_realFile(777,'convergeddPdF'//trim(rankStr),modelName,size(crystallite_dPdF0))
   read (777,rec=1) crystallite_dPdF0
   close (777)

   call IO_read_realFile(777,'convergedTstar'//trim(rankStr),modelName,size(crystallite_Tstar0_v))
   read (777,rec=1) crystallite_Tstar0_v
   close (777)

   call IO_read_realFile(777,'convergedStateConst'//trim(rankStr),modelName)
   m = 0_pInt
   readPlasticityInstances: do ph = 1_pInt, size(phase_plasticity)
     do k = 1_pInt, plasticState(ph)%sizeState
       do l = 1, size(plasticState(ph)%state0(1,:))
         m = m+1_pInt
         read(777,rec=m) plasticState(ph)%state0(k,l)
     enddo; enddo
   enddo readPlasticityInstances
   close (777)

   call IO_read_realFile(777,'convergedStateHomog'//trim(rankStr),modelName)
   m = 0_pInt
   readHomogInstances: do homog = 1_pInt, material_Nhomogenization
     do k = 1_pInt, homogState(homog)%sizeState
       do l = 1, size(homogState(homog)%state0(1,:))
         m = m+1_pInt
         read(777,rec=m) homogState(homog)%state0(k,l)
     enddo; enddo
   enddo readHomogInstances
   close (777)

   restartRead = .false.
 endif
 flush(6)

end subroutine CPFEM_init


!--------------------------------------------------------------------------------------------------
!> @brief perform initialization at first call, update variables and call the actual material model
!--------------------------------------------------------------------------------------------------
subroutine CPFEM_general(age, dt)
 use prec, only: &
   pReal, &
   pInt
 use numerics, only: &
   worldrank
 use debug, only: &
   debug_level, &
   debug_CPFEM, &
   debug_levelBasic, &
   debug_levelExtensive, &
   debug_levelSelective
 use FEsolving, only: &
   terminallyIll, &
   restartWrite
 use math, only: &
   math_identity2nd, &
   math_mul33x33, &
   math_det33, &
   math_transpose33, &
   math_I3, &
   math_Mandel3333to66, &
   math_Mandel66to3333, &
   math_Mandel33to6, &
   math_Mandel6to33
 use material, only: &
   plasticState, &
   sourceState, &
   homogState, &
   thermalState, &
   damageState, &
   vacancyfluxState, &
   hydrogenfluxState, &
   material_phase, &
   phase_plasticity, &
   phase_Nsources, &
   material_Nhomogenization
 use crystallite, only: &
   crystallite_partionedF,&
   crystallite_F0, &
   crystallite_Fp0, &
   crystallite_Fp, &
   crystallite_Fi0, &
   crystallite_Fi, &
   crystallite_Lp0, &
   crystallite_Lp, &
   crystallite_Li0, &
   crystallite_Li, &
   crystallite_dPdF0, &
   crystallite_dPdF, &
   crystallite_Tstar0_v, &
   crystallite_Tstar_v
 use homogenization, only: &  
   materialpoint_F, &
   materialpoint_F0, &
   materialpoint_stressAndItsTangent, &
   materialpoint_postResults
 use IO, only: &
   IO_write_jobRealFile, &
   IO_warning
 use DAMASK_interface

 implicit none
 real(pReal), intent(in) ::                          dt                                         !< time increment
 logical, intent(in) ::                              age                                        !< age results

 integer(pInt) ::                                    i, k, l, m, ph, homog, mySource
 character(len=1024) :: rankStr

 !*** age results and write restart data if requested
 if (age) then
   crystallite_F0  = crystallite_partionedF                                                    ! crystallite deformation (_subF is perturbed...)
   crystallite_Fp0 = crystallite_Fp                                                            ! crystallite plastic deformation
   crystallite_Lp0 = crystallite_Lp                                                            ! crystallite plastic velocity
   crystallite_Fi0 = crystallite_Fi                                                            ! crystallite intermediate deformation
   crystallite_Li0 = crystallite_Li                                                            ! crystallite intermediate velocity
   crystallite_dPdF0 = crystallite_dPdF                                                        ! crystallite stiffness
   crystallite_Tstar0_v = crystallite_Tstar_v                                                  ! crystallite 2nd Piola Kirchhoff stress

   forall ( i = 1:size(plasticState    )) plasticState(i)%state0     = plasticState(i)%state   ! copy state in this lenghty way because: A component cannot be an array if the encompassing structure is an array
   do i = 1, size(sourceState)
     do mySource = 1,phase_Nsources(i)
       sourceState(i)%p(mySource)%state0 = sourceState(i)%p(mySource)%state                    ! copy state in this lenghty way because: A component cannot be an array if the encompassing structure is an array
   enddo; enddo
   if (iand(debug_level(debug_CPFEM), debug_levelBasic) /= 0_pInt) &
     write(6,'(a)') '<< CPFEM >> aging states'

   do homog = 1_pInt, material_Nhomogenization
     homogState       (homog)%state0 =  homogState       (homog)%state
     thermalState     (homog)%state0 =  thermalState     (homog)%state
     damageState      (homog)%state0 =  damageState      (homog)%state
     vacancyfluxState (homog)%state0 =  vacancyfluxState (homog)%state
     hydrogenfluxState(homog)%state0 =  hydrogenfluxState(homog)%state
   enddo


   if (restartWrite) then
     if (iand(debug_level(debug_CPFEM), debug_levelBasic) /= 0_pInt) &
        write(6,'(a)') '<< CPFEM >> writing state variables of last converged step to binary files'
   
     write(rankStr,'(a1,i0)')'_',worldrank

     call IO_write_jobRealFile(777,'recordedPhase'//trim(rankStr),size(material_phase))
     write (777,rec=1) material_phase
     close (777)

     call IO_write_jobRealFile(777,'convergedF'//trim(rankStr),size(crystallite_F0))
     write (777,rec=1) crystallite_F0
     close (777)

     call IO_write_jobRealFile(777,'convergedFp'//trim(rankStr),size(crystallite_Fp0))
     write (777,rec=1) crystallite_Fp0
     close (777)

     call IO_write_jobRealFile(777,'convergedFi'//trim(rankStr),size(crystallite_Fi0))
     write (777,rec=1) crystallite_Fi0
     close (777)

     call IO_write_jobRealFile(777,'convergedLp'//trim(rankStr),size(crystallite_Lp0))
     write (777,rec=1) crystallite_Lp0
     close (777)

     call IO_write_jobRealFile(777,'convergedLi'//trim(rankStr),size(crystallite_Li0))
     write (777,rec=1) crystallite_Li0
     close (777)

     call IO_write_jobRealFile(777,'convergeddPdF'//trim(rankStr),size(crystallite_dPdF0))
     write (777,rec=1) crystallite_dPdF0
     close (777)

     call IO_write_jobRealFile(777,'convergedTstar'//trim(rankStr),size(crystallite_Tstar0_v))
     write (777,rec=1) crystallite_Tstar0_v
     close (777)

     call IO_write_jobRealFile(777,'convergedStateConst'//trim(rankStr))
     m = 0_pInt
     writePlasticityInstances: do ph = 1_pInt, size(phase_plasticity)
       do k = 1_pInt, plasticState(ph)%sizeState
         do l = 1, size(plasticState(ph)%state0(1,:))
           m = m+1_pInt
           write(777,rec=m) plasticState(ph)%state0(k,l)
       enddo; enddo
     enddo writePlasticityInstances
     close (777)

     call IO_write_jobRealFile(777,'convergedStateHomog'//trim(rankStr))
     m = 0_pInt
     writeHomogInstances: do homog = 1_pInt, material_Nhomogenization
       do k = 1_pInt, homogState(homog)%sizeState
         do l = 1, size(homogState(homog)%state0(1,:))
           m = m+1_pInt
           write(777,rec=m) homogState(homog)%state0(k,l)
       enddo; enddo
     enddo writeHomogInstances
     close (777)

   endif
 endif

 if (.not. terminallyIll) &
   call materialpoint_stressAndItsTangent(.True., dt)

end subroutine CPFEM_general

end module CPFEM2
