!********************************************************************
! Material subroutine for MSC.Marc Version 0.1
!
! written by F. Roters, P. Eisenlohr, L. Hantcherli, W.A. Counts
! MPI fuer Eisenforschung, Duesseldorf
!
! last modified: 28.03.2007
!********************************************************************
!     Usage:
!             - choose material as hypela2
!             - set statevariable 2 to index of material
!             - set statevariable 3 to index of texture
!             - choose output of user variables if desired
!             - make sure the file "mattex.mpie" exists in the working
!               directory
!             - use nonsymmetric option for solver (e.g. direct 
!               profile or multifrontal sparse, the latter seems
!               to be faster!)
!********************************************************************
!     Marc subroutines used:
!             - hypela2
!             - plotv
!             - quit
!********************************************************************
!     Marc common blocks included:
!             - concom: lovl, ncycle, inc, incsub
!             - creeps: timinc
!********************************************************************
!
 include "prec.f90"
 include "math.f90"
 include "IO.f90"
 include "mesh.f90"
 include "constitutive.f90"
 include "CPFEM.f90"
!
 SUBROUTINE hypela2(d,g,e,de,s,t,dt,ngens,n,nn,kcus,matus,ndi,&
                    nshear,disp,dispt,coord,ffn,frotn,strechn,eigvn,ffn1,&
                    frotn1,strechn1,eigvn1,ncrd,itel,ndeg,ndm,&
                    nnode,jtype,lclass,ifr,ifu)
!********************************************************************
! This is the Marc material routine
!********************************************************************
!
! *************   user subroutine for defining material behavior  **************
!
!
! CAUTION : Due to calculation of the Deformation gradients, Stretch Tensors and
!         Rotation tensors at previous and current states, the analysis can be
!         computationally expensive. Please use the user subroutine ->  hypela
!         if these kinematic quantities are not needed in the constitutive model
!
!
! IMPORTANT NOTES :
!
! (1) F,R,U are only available for continuum and membrane elements (not for
!     shells and beams).
!
! (2) For total Lagrangian formulation use the -> 'Elasticity,1' card(=
!     total Lagrange with large disp) in the parameter section of input deck.
!     For updated Lagrangian formulation use the -> 'Plasticity,3' card(=
!     update+finite+large disp+constant d) in the parameter section of
!     input deck.
!
!
!     d            stress strain law to be formed
!     g            change in stress due to temperature effects
!     e            total elastic strain
!     de           increment of strain
!     s            stress - should be updated by user
!     t            state variables (comes in at t=n, must be updated
!                                   to have state variables at t=n+1)
!     dt           increment of state variables
!     ngens        size of stress - strain law
!     n            element number
!     nn           integration point number
!     kcus(1)      layer number
!     kcus(2)      internal layer number
!     matus(1)     user material identification number
!     matus(2)     internal material identification number
!     ndi          number of direct components
!     nshear       number of shear components
!     disp         incremental displacements
!     dispt        displacements at t=n   (at assembly,        lovl=4) and
!                  displacements at t=n+1 (at stress recovery, lovl=6)
!     coord        coordinates
!     ncrd         number of coordinates
!     ndeg         number of degrees of freedom
!     itel         dimension of F and R, either 2 or 3
!     nnode        number of nodes per element
!     jtype        element type
!     lclass       element class
!     ifr          set to 1 if R has been calculated
!     ifu          set to 1 if strech has been calculated
!
!     at t=n   :
!
!     ffn          deformation gradient
!     frotn        rotation tensor
!     strechn      square of principal stretch ratios, lambda(i)
!     eigvn(i,j)   i principal direction components for j eigenvalues
!
!     at t=n+1 :
!
!     ffn1         deformation gradient
!     frotn1       rotation tensor
!     strechn1     square of principal stretch ratios, lambda(i)
!     eigvn1(i,j)  i principal direction components for j eigenvalues
!
!     The following operation obtains U (stretch tensor) at t=n+1 :
!
!     call scla(un1,0.d0,itel,itel,1)
!     do 3 k=1,3
!      do 2 i=1,3
!       do 1 j=1,3
!        un1(i,j)=un1(i,j)+dsqrt(strechn1(k))*eigvn1(i,k)*eigvn1(j,k)
!1      continue
!2     continue
!3    continue
!
!
 use prec,  only: pReal,pInt
 use CPFEM, only: CPFEM_stress_all, CPFEM_jaco_old
 use math,  only: invnrmMandel, nrmMandel
 implicit real(pReal) (a-h,o-z)
!
! Marc common blocks are in fixed format so they have to be pasted in here
! Beware of changes in newer Marc versions -- these are from 2005r3
! concom is needed for inc, subinc, ncycle, lovl
! include 'concom'
 common/concom/ &
     iacous, iasmbl, iautth,    ibear,  icompl,     iconj,  icreep, ideva(50), idyn,   idynt,&
     ielas,  ielcma, ielect,    iform,  ifour,      iharm,  ihcps,  iheat,     iheatt, ihresp,&
     ijoule, ilem,   ilnmom,    iloren, inc,        incext, incsub, ipass,     iplres, ipois,&
     ipoist, irpflo, ismall,    ismalt, isoil,      ispect, ispnow, istore,    iswep,  ithcrp,&
     itherm, iupblg, iupdat,    jacflg, jel,        jparks, largst, lfond,     loadup, loaduq,&
     lodcor, lovl,   lsub,      magnet, ncycle,     newtnt, newton, noshr,     linear, ivscpl,&
     icrpim, iradrt, ipshft,    itshr,  iangin,     iupmdr, iconjf, jincfl,    jpermg, jhour,&
     isolvr, jritz,  jtable,    jshell, jdoubl,     jform,  jcentr, imini,     kautth, iautof,&
     ibukty, iassum, icnstd,    icnstt, kmakmas,    imethvp,iradrte,iradrtp,   iupdate,iupdatp,&
     ncycnt, marmen ,idynme,    ihavca, ispf,       kmini,  imixed, largtt,    kdoela, iautofg,&
     ipshftp,idntrc, ipore,     jtablm, jtablc,     isnecma,itrnspo,imsdif,    jtrnspo,mcnear,&
     imech,  imecht, ielcmat,   ielectt,magnett,    imsdift,noplas, jtabls,    jactch, jtablth,&
     kgmsto ,jpzo,   ifricsh,   iremkin,iremfor,    ishearp,jspf,   machining, jlshell,icompsol,&
     iupblgfo,jcondir,nstcrp,   nactive,ipassref,   nstspnt,ibeart,icheckmpc,  noline, icuring,&
     ishrink,ioffsflg,isetoff,  iharmt, inc_incdat, iautspc,ibrake
! creeps is needed for timinc (time increment)
! include 'creeps'
 common/creeps/ &
     cptim,timinc,timinc_p,timinc_s,timincm,timinc_a,timinc_b,creept(33),icptim,icfte,icfst,&
     icfeq,icftm,icetem,mcreep,jcreep,icpa,icftmp,icfstr,icfqcp,icfcpm,icrppr,icrcha,icpb,iicpmt,iicpa
!
 integer(pInt) cp_en, i
!
 dimension e(*),de(*),t(*),dt(*),g(*),d(ngens,*),s(*), n(2),coord(ncrd,*),disp(ndeg,*),matus(2),dispt(ndeg,*),ffn(itel,*),&
           frotn(itel,*),strechn(itel),eigvn(itel,*),ffn1(itel,*),frotn1(itel,*),strechn1(itel),eigvn1(itel,*),kcus(2)
!
! call general material routine only in increment 0 and for lovl==6 (stress recovery)
          
!     subroutine cpfem_general(mpie_ffn, mpie_ffn1, mpie_cn, mpie_tinc, mpie_enp, mpie_in)
!********************************************************************
!     This routine calculates the material behaviour
!********************************************************************
!     mpie_ffn         deformation gradient for t=t0
!     mpie_ffn1        deformation gradient for t=t1
!     mpie_cn          number of cycle
!     mpie_tinc        time increment
!     mpie_en          element number
!     mpie_in          intergration point number
!********************************************************************
 cp_en = mesh_FEasCP('elem', n(1))
 if ((lovl==6).or.(inc==0)) then
    call CPFEM_general(ffn, ffn1, inc, incsub, ncycle, timinc, cp_en, nn)
 endif
! return stress and jacobi
!     Mandel: 11, 22, 33, 12, 23, 13 
!     Marc:   11, 22, 33, 12, 23, 13
 s(1:ngens)=invnrmMandel(1:ngens)*CPFEM_stress_all(1:ngens, nn, cp_en)
 d(1:ngens,1:ngens)=CPFEM_jaco_old(1:ngens,1:ngens, nn, cp_en)
 forall(i=1:ngens) d(i,1:ngens)=d(i,1:ngens)*invnrmMandel(1:ngens)
 return
 
 END SUBROUTINE
!
!
 SUBROUTINE plotv(v,s,sp,etot,eplas,ecreep,t,m,nn,layer,ndi,nshear,jpltcd)
!********************************************************************
!     This routine sets user defined output variables for Marc
!********************************************************************
!
!     select a variable contour plotting (user subroutine).
!
!     v            variable
!     s (idss)         stress array
!     sp           stresses in preferred direction
!     etot          total strain (generalized)
!     eplas         total plastic strain
!     ecreep        total creep strain
!     t             current temperature
!     m            element number
!     nn           integration point number
!     layer        layer number
!     ndi (3)       number of direct stress components
!     nshear (3)    number of shear stress components
!
!********************************************************************
 use prec,  only: pReal,pInt
 use CPFEM, only: CPFEM_results, CPFEM_Nresults
 use constitutive, only: constitutive_maxNresults
 use mesh,  only: mesh_FEasCP
 implicit none
!
 real(pReal) s(*),etot(*),eplas(*),ecreep(*),sp(*)
 real(pReal) v, t(*)
 integer(pInt) m, nn, layer, ndi, nshear, jpltcd
!
! assign result variable
 v=CPFEM_results(mod(jpltcd, CPFEM_Nresults+constitutive_maxNresults),&
                 int(jpltcd/(CPFEM_Nresults+constitutive_maxNresults)),&
                 nn, mesh_FEasCP('elem', m))
 return
 END SUBROUTINE
!
!
! subroutine utimestep(timestep,timestepold,icall,time,timeloadcase)
!********************************************************************
!     This routine modifies the addaptive time step of Marc
!********************************************************************
! use prec, only: pReal,pInt
! use CPFEM, only : CPFEM_timefactor_max
! implicit none
!
! real(pReal) timestep, timestepold, time,timeloadcase 
! integer(pInt) icall
!
! user subroutine for modifying the time step in auto step
!
!   timestep    :  the current time step as suggested by marc
!                  to be modified in this routine
!   timestepold :  the current time step before it was modified by marc
!   icall       :  =1 for setting the initial time step
!                  =2 if this routine is called during an increment
!                  =3 if this routine is called at the beginning
!                     of the increment
!   time        :  time at the start of the current increment
!   timeloadcase:  time period of the current load case
!
!   it is in general not recommended to increase the time step
!   during the increment.
!   this routine is called right after the time step has (possibly)
!   been updated by marc.
!
!  user coding
!     reduce timestep during increment in case mpie_timefactor is too large
! if(icall==2_pInt) then
!    if(mpie_timefactor_max>1.25_pReal) then
!        timestep=min(timestep,timestepold*0.8_pReal)
!    end if
! return
!     modify timestep at beginning of new increment
! else if(icall==3_pInt) then
!    if(mpie_timefactor_max<=0.8_pReal) then
!        timestep=min(timestep,timestepold*1.25_pReal)
!    else if (mpie_timefactor_max<=1.0_pReal) then
!        timestep=min(timestep,timestepold/mpie_timefactor_max)
!    else if (mpie_timefactor_max<=1.25_pReal) then
!        timestep=min(timestep,timestepold*1.01_pReal)
!    else
!        timestep=min(timestep,timestepold*0.8_pReal)
!    end if
! end if
! return
! end
