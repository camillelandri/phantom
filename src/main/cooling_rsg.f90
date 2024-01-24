!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2023 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module cooling_rsg
!
! Simple beta-cooling prescription used for experiments on gravitational
!  instability in discs
!
! :References:
!   Gammie (2001), ApJ 553, 174-183
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - beta_cool : *beta factor in Gammie (2001) cooling*
!
! :Dependencies: infile_utils, io
!
 implicit none
 real, public :: Tdust  = 1100.
 real, public :: pdust  = -0.9
 real, public :: Tstar  = 3500.
 real, public :: Rstar  = 1500.
 real, public :: Mstar  = 20.
 real, public :: radacc  = 1.

contains
!-----------------------------------------------------------------------
!+
!   Explicit cooling based on equilibrium temperature
!+
!-----------------------------------------------------------------------
subroutine cooling_rsg_explicit(xi,yi,zi,ui,dudti,dt,T_on_u)
 use io,                only:fatal, warning
 real, intent(in)    :: ui,xi,yi,zi
 real, intent(in)    :: T_on_u,dt
 real, intent(inout) :: dudti
 

 real :: r2,ueq,Teq

 r2 = xi*xi + yi*yi + zi*zi
 ! Get temperature due to stellar radiation with geometrical dilution factor (7.36 in Lamers, Cassenilli 1999)
 ! for now p=-0.9, Tc = 1100K for silicate grains (Bladh Hofner 2012)
 Teq = Tstar * (0.5 * (1 - sqrt(1 - Rstar**2/r2)))**(1/(4+pdust)) ! eq tempature when taking stellar radation into account with a dilution factor
 ueq = Teq / T_on_u
 if (ui < 0.99*ueq) then
   dudti  = dudti  - (ui-ueq) / dt
   !call warning('cooling','ueq larger than ui',var='dudti',val=dudti)
 else
     dudti  = dudti - (ui-ueq) / dt
 endif
end subroutine cooling_rsg_explicit

!-----------------------------------------------------------------------
!+
!  writes input options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_cooling_rsg(iunit)
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit

 call write_inopt(Tdust,'Tdust','condensation temperature of dust, in K, (Bladh Hofner 2012 for values)',iunit)
 call write_inopt(pdust,'pdust','opacity exponent for dust [Kappa(lambda^-pdust)], (Bladh Hofner 2012 for values)',iunit)
 call write_inopt(Tstar,'Tstar','Photosphere temperature of primary, in K',iunit)
 call write_inopt(Rstar,'Rstar','Photosphere radius of primary, in Rsun',iunit)
 call write_inopt(Mstar,'Mstar','Mass of primary, in Msun',iunit)
 call write_inopt(radacc,'radacc','factor for the wind acceleration dur to rad. pres. (free wind = 1)',iunit)                                    

end subroutine write_options_cooling_rsg

!-----------------------------------------------------------------------
!+
!  reads input options from the input file
!+
!-----------------------------------------------------------------------
subroutine read_options_cooling_rsg(name,valstring,imatch,igotall,ierr)
 use io, only:fatal, warning
 character(len=*), intent(in)  :: name,valstring
 logical,          intent(out) :: imatch,igotall
 integer,          intent(out) :: ierr
 integer, save :: ngot = 0

 imatch  = .true.
 igotall = .true. ! none of the cooling options are compulsory
 select case(trim(name))
 case('Tdust')
    read(valstring,*,iostat=ierr) Tdust
    ngot = ngot + 1
    if (Tdust < 0.) call warning('read_options','dust temperature must be >= 0')
 case('pdust')
      read(valstring,*,iostat=ierr) pdust
      ngot = ngot + 1
 case('Tstar')
      read(valstring,*,iostat=ierr) Tstar
      ngot = ngot + 1
      if (Tstar < 0.) call warning('read_options','star temperature must be >= 0')
 case('Rstar')
      read(valstring,*,iostat=ierr) Rstar
      ngot = ngot + 1
      if (Rstar < 0.) call warning('read_options','star radius must be >= 0')
 case('Mstar')
      read(valstring,*,iostat=ierr) Mstar
      ngot = ngot + 1
      if (Mstar < 0.) call warning('read_options','star mass must be >= 0')
 case('radacc')
          read(valstring,*,iostat=ierr) radacc
          ngot = ngot + 1
          if (radacc < 0.) call warning('read_options','radiative acc must be >= 0')
 case default
    imatch = .false.
 end select
 if (ngot >= 1) igotall = .true.

end subroutine read_options_cooling_rsg

end module cooling_rsg
