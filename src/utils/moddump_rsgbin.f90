!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2023 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module moddump
!
! RSG binary with secondary grazing the envelope; adapted from moddump_binary
! and moddump_sink
!
! :References: None
!
! :Owner: Camille Landri
!
! :Runtime parameters: None
!
! :Dependencies: centreofmass, dim, extern_corotate, externalforces,
!   infile_utils, io, options, part, physcon, prompting, readwrite_dumps,
!   rho_profile, setbinary, table_utils, timestep, units, vectorutils
!
 implicit none

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use part,              only:nptmass,xyzmh_ptmass,vxyz_ptmass,ihacc,ihsoft,igas,&
                             delete_dead_or_accreted_particles,mhd,rhoh,shuffle_part,&
                             kill_particle,copy_particle
 use setbinary,         only:set_binary
 use units,             only:umass,udist,utime
 use physcon,           only:au,solarm,solarr,gg,pi
 use centreofmass,      only:reset_centreofmass,get_centreofmass
 use prompting,         only:prompt
 use options,           only:iexternalforce
 use externalforces,    only:omega_corotate,iext_corotate
 use extern_corotate,   only:icompanion_grav,companion_xpos,companion_mass,primarycore_xpos,&
                             primarycore_mass,primarycore_hsoft,hsoft
 use infile_utils,      only:open_db_from_file,inopts,read_inopt,close_db
 use table_utils,       only:yinterp
 use rho_profile,       only:read_mesa
 use dim,               only:maxptmass,maxp
 use io,                only:fatal,idisk1,iprint
 use timestep,          only:tmax,dtmax
 use readwrite_dumps,   only:read_dump

 integer, intent(inout)    :: npart
 integer, intent(inout)    :: npartoftype(:)
 real,    intent(inout)    :: massoftype(:)
 real,    intent(inout)    :: xyzh(:,:),vxyzu(:,:)
 integer                   :: i,ierr,setup_case,ioption=1,irhomax,n
 integer                   :: iremove = 2
 integer                   :: nstar1,nstar2,nptmass1
 real                      :: primary_mass,companion_mass_1,companion_mass_2,mass_ratio,m1,a,hsoft2,pmass1,pmass2
 real                      :: mass_donor,separation,newCoM,period,m2,primarycore_xpos_old
 real                      :: a1,a2,e,vr,hsoft_default = 3.
 real                      :: hacc1,hacc2,hacc3,hsoft_primary,mcore,comp_shift=100,sink_dist,vel_shift
 real                      :: mcut,rcut,Mstar,radi,rhopart,rhomax = 0.0
 real                      :: time2,hfact2,Rstar
 real, allocatable         :: r(:),den(:),pres(:),temp(:),enitab(:),Xfrac(:),Yfrac(:),m(:)
 logical                   :: corotate_answer,iprimary_grav_ans
 character(len=20)         :: filename = 'binary.in'
 character(len=100)        :: densityfile,dumpname
 type(inopts), allocatable :: db(:)
 integer                   :: isinkpart
 real                      :: racc,mass,mass_old,newx
 logical                   :: iresetCM


 if (nptmass > 3) then
    call fatal('moddump_binary','Number of sink particles > 3')
 elseif (nptmass == 3) then
    print*, 'Three sink particles are present. Choose option below:'
    print "(1(/,a))",'1) Remove a sink from the simulation'
    call prompt('Select option above : ',ioption)
    select case(ioption)

    case(1)
       do i=1,nptmass
          write(*,'(A,I2,A,ES10.3,A,ES10.3)') 'Point mass ',i,': M = ',xyzmh_ptmass(4,i),&
                                              ' and radial position = ',sqrt(dot_product(xyzmh_ptmass(1:3,i),xyzmh_ptmass(1:3,i)))
       enddo
       call prompt('Which sink would you like to remove : ',iremove)
       if (iremove == 3) then
          xyzmh_ptmass(:,iremove) = 0.
          vxyz_ptmass(:,iremove) = 0.
          nptmass = 2
       elseif (iremove == 2) then
          xyzmh_ptmass(:,2) = xyzmh_ptmass(:,3)
          vxyz_ptmass(:,2) = vxyz_ptmass(:,3)
          nptmass = 2
       endif
    end select

 elseif (nptmass == 2) then
    print*, 'Two sinks particles are present. Choose option below:'
    print "(4(/,a))",'1) Transform from corotating frame to inertial frame', &
                     '2) Shift companion position in the co-rotating frame', &
                     '3) Add velocity to companion', &
                     '4) (Re)set sink properties'
    call prompt('Select option above : ',ioption)
    select case(ioption)
    case(1)
       call prompt('Please write the name of the input file : ',filename)
       call open_db_from_file(db,filename,20,ierr)
       call read_inopt(omega_corotate,'omega_corotate',db)
       call close_db(db)
       call transform_from_corotating_to_inertial_frame(xyzh,vxyzu,npart,nptmass,omega_corotate,xyzmh_ptmass,vxyz_ptmass)

    case(2)
       call prompt('How many code units to shift companion (+ve is towards primary)?',comp_shift)
       sink_dist = sqrt((xyzmh_ptmass(1,1)-xyzmh_ptmass(1,2))**2 &
                      + (xyzmh_ptmass(2,1)-xyzmh_ptmass(2,2))**2 &
                      + (xyzmh_ptmass(3,1)-xyzmh_ptmass(3,2))**2)

       xyzmh_ptmass(1,2) = -(comp_shift/sink_dist * (xyzmh_ptmass(1,2)-xyzmh_ptmass(1,1)) - xyzmh_ptmass(1,2))
       xyzmh_ptmass(2,2) = -(comp_shift/sink_dist * (xyzmh_ptmass(2,2)-xyzmh_ptmass(2,1)) - xyzmh_ptmass(2,2))
       xyzmh_ptmass(3,2) = -(comp_shift/sink_dist * (xyzmh_ptmass(3,2)-xyzmh_ptmass(3,1)) - xyzmh_ptmass(3,2))

       call reset_centreofmass(npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)
       iexternalforce = iext_corotate
       omega_corotate = sqrt((sink_dist-comp_shift)*(xyzmh_ptmass(4,1)+xyzmh_ptmass(4,2)))/(sink_dist-comp_shift)**2

       do i=1,npart
          vxyzu(1,i) = 0.0
          vxyzu(2,i) = 0.0
          vxyzu(3,i) = 0.0
       enddo

       do i=1,nptmass
          vxyz_ptmass(1,i) = 0.0
          vxyz_ptmass(1,i) = 0.0
          vxyz_ptmass(1,i) = 0.0
       enddo

    case(3)
       sink_dist = sqrt((xyzmh_ptmass(1,1)-xyzmh_ptmass(1,2))**2 &
                 + (xyzmh_ptmass(2,1)-xyzmh_ptmass(2,2))**2 &
                 + (xyzmh_ptmass(3,1)-xyzmh_ptmass(3,2))**2)

       print*, utime, umass, udist
       vel_shift = sink_dist/(5.0*60*60*24*365/utime)

       call prompt('Give velocity to add in direction of the primary : ',vel_shift, 0.0)

       vxyz_ptmass(1,2) = -(vel_shift/sink_dist * (xyzmh_ptmass(1,2)-xyzmh_ptmass(1,1)) - vxyz_ptmass(1,2))
       vxyz_ptmass(2,2) = -(vel_shift/sink_dist * (xyzmh_ptmass(2,2)-xyzmh_ptmass(2,1)) - vxyz_ptmass(2,2))
       vxyz_ptmass(3,2) = -(vel_shift/sink_dist * (xyzmh_ptmass(3,2)-xyzmh_ptmass(3,1)) - vxyz_ptmass(3,2))

       call reset_centreofmass(npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)

    case(4)
       call set_sinkproperties(xyzmh_ptmass)
    end select

 else  ! One or fewer point masses
    !choose what to do with the star: set a binary or setup a magnetic field
   print*, 'Set up a binary system by adding a sink companion'
    setup_case = 1

    select case(setup_case)
    case(1,3)
       ! set binary defaults
       companion_mass_1 = 1.0
       a1 = 7000.
       e = 0.8
       mcore = 0.
       hacc1 = 0.
       hacc2 = 20.
       vr = 0.

       ! find current stellar radius
       Rstar = 0.
       do i = 1,npart
          Rstar  = max(Rstar,sqrt(dot_product(xyzh(1:3,i),xyzh(1:3,i))))
       enddo
       print*, 'ptmass number:', nptmass
       print*, 'Current mass unit is ', umass,'g):'
       pmass1 = massoftype(igas)
       print*, 'Current particle mass in code units are ', pmass1
       call prompt('Enter companion mass in code units',companion_mass_1,0.) ! For case 8, eventually want to read mass of star 2 from header instead of prompting it
       print*, 'Current length unit is ', udist ,'cm):'
       print*, 'Current stellar radius in code units is ', Rstar
       call prompt('Enter orbit semi-major axis in code units', a1, 0.)
       call prompt('Enter orbit eccentricity', e, 0., 1.)
       call prompt('Enter companion radial velocity', vr)

       if (nptmass == 1) then ! there is a sink stellar core
          mcore = xyzmh_ptmass(4,1)
          hacc1 = xyzmh_ptmass(ihacc,1)
          print*, 'Current accretion radius of primary core is ', hacc1,' code units'
          call prompt('Enter accretion radius for the primary core in code units', hacc1, 0.)
          hacc2 = 0.
          call prompt('Enter accretion radius for the companion in code units', hacc2, 0.)
       endif

       corotate_answer = .false.
       call prompt('Do you want to transform to a corotating frame and simulate corotating binary?', corotate_answer)
       call reset_centreofmass(npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)

       !removes the dead or accreted particles for a correct total mass computation
       call delete_dead_or_accreted_particles(npart,npartoftype)
       print*,' Got ',npart,npartoftype(igas),' after deleting accreted particles'

       !sets up the binary system orbital parameters
       hsoft_primary = 0.
       if (nptmass == 1) then
          mcore = xyzmh_ptmass(4,1)
          hsoft_primary = xyzmh_ptmass(ihsoft,1)  ! stash primary core hsoft before calling set_binary, which resets the softening lengths
       elseif (nptmass == 0) then
          mcore = 0.
       else
          call fatal('moddump_binary', 'sink particle not specified (nptmass > 1)')
       endif

       primary_mass = npartoftype(igas) * massoftype(igas) + mcore
       print*, 'Current primary mass in code units is ',primary_mass

       ! set the binary
       if (corotate_answer) then ! corotating frame
          iexternalforce = iext_corotate  !turns on corotation
          call set_binary(primary_mass,companion_mass_1,a1,e,hacc1,hacc2,xyzmh_ptmass,vxyz_ptmass,nptmass,ierr,omega_corotate)
          print "(/,a,es18.10,/)", ' The angular velocity in the corotating frame is: ', omega_corotate

          ! set all the gas velocities in corotating frame to 0, implying that the binary is corotating
          ! at the moment, only a corotating binary can be set up in the corotating frame
          do i=1,npart
             vxyzu(1:3,i) = 0.
          enddo
       else ! non corotating frame
          call set_binary(primary_mass,companion_mass_1,a1,e,hacc1,hacc2,xyzmh_ptmass,vxyz_ptmass,nptmass,ierr)
          ! sink no. 2 & 3 are created by "set_binary" in the ptmass arrays
       endif

       if (nptmass == 3) then ! if original star has a point mass core
          !move primary core from pos 2 to 1
          xyzmh_ptmass(1:3,1) = xyzmh_ptmass(1:3,2)
          vxyz_ptmass(1:3,1) = vxyz_ptmass(1:3,2)

          !move companion point mass from pos 3 to 2
          xyzmh_ptmass(:,2) = xyzmh_ptmass(:,3)
          vxyz_ptmass(1:3,2) = vxyz_ptmass(1:3,3)
          vxyz_ptmass(1,2) = vxyz_ptmass(1,2) + vr

          if (setup_case == 1) then  ! Assume companion should be a sink particle
             nptmass = nptmass - 1  ! Delete point mass 3 (duplicate of companion)
          elseif (setup_case == 8) then  ! Companion does not contain a sink particle and is read from second dumpfile
             nptmass = nptmass - 2  ! Delete point masses 2 and 3, leaving just the primary core
          endif

          !takes necessary inputs from user 2 (the softening lengths for the sinks have to be taken in input after using the "set_binary" function since it resets them)
          xyzmh_ptmass(ihsoft,1) = hsoft_primary
          print*, 'Current softening length of the primary core is ', xyzmh_ptmass(ihsoft,1),' code units'
          if (setup_case == 1) call prompt('Enter softening length for companion',xyzmh_ptmass(ihsoft,2),0.)

       elseif (nptmass == 2) then ! if original star is coreless
          ! Just need to delete both point masses
          nptmass = 0
       endif

       if (setup_case == 1) then
          !shifts gas to the primary point mass created in 'set_binary'
          do i=1,npart
             xyzh(1:3,i) = xyzh(1:3,i) + xyzmh_ptmass(1:3,1)
             vxyzu(1:3,i) = vxyzu(1:3,i) + vxyz_ptmass(1:3,1)
          enddo
       elseif (setup_case == 8) then
          nstar1 = npart ! save npart in star 1
          nptmass1 = nptmass  ! stash nptmass for dump 1, as read_dumps overwrites it
          dumpname = ''
          call prompt('Enter name of second dumpfile',dumpname)
          nstar2 = nstar1
          call prompt('Enter no. of particles in second dumpfile',nstar2)

          ! Move star 1 particles to avoid getting overwritten when reading second dump file.
          if (2*nstar1 > maxp) then  ! Check if particle array is large enough to provide particle-copying buffer
             call fatal('moddump_binary','Two times number of particles in star 1 exceeds MAXP. Need to compile with larger MAXP')
          endif
          if (nstar1 > nstar2) then ! Move ith particle of star 1 to nstar1+i
             do i=1,nstar1
                call copy_particle(i,nstar1+i,.false.)
             enddo
          else ! Move ith particle of star 1 to nstar2+i
             do i=1,nstar1
                call copy_particle(i,nstar2+i,.false.)
             enddo
          endif

          ! read dump file containing star 2
          call read_dump(trim(dumpname),time2,hfact2,idisk1+1,iprint,0,1,ierr)
          nptmass = nptmass1 + nptmass  ! set nptmass to be sum of nptmass in dump 1 and dump 2
          pmass2 = massoftype(igas)
          if (ierr /= 0) call fatal('read_dump','error reading second dump file')
          if ( abs(1.-pmass2/pmass1) > 1.e-3) then
             call fatal('moddump_binary','unequal mass particles between dumps 1 and 2, pmass2 /= pmass1')
          endif
          print*,'Setting gas mass to be that from first dump,',pmass1
          massoftype(igas) = pmass1

          if (nstar1 > nstar2) then ! Move ith particle of star 1 to nstar2+i
             do i=1,nstar1
                call copy_particle(nstar1+i,nstar2+i,.false.)
             enddo
          endif

          npart = nstar1 + nstar2
          npartoftype(igas) = npart

          ! shift star2 to secondary point mass (deleted)
          do i=1,nstar2
             xyzh(1:3,i) = xyzh(1:3,i) + xyzmh_ptmass(1:3,2)
             vxyzu(1:3,i) = vxyzu(1:3,i) + vxyz_ptmass(1:3,2)
          enddo
          ! shift star1 to primary point mass (deleted)
          do i=nstar2+1,npart
             xyzh(1:3,i) = xyzh(1:3,i) + xyzmh_ptmass(1:3,1)
             vxyzu(1:3,i) = vxyzu(1:3,i) + vxyz_ptmass(1:3,1)
          enddo

       endif
       call reset_centreofmass(npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)
       
    end select
 endif

 return
end subroutine modify_dump

subroutine set_sinkproperties(xyzmh_ptmass)

 use part,       only:nptmass,ihacc,ihsoft,igas,imacc,ilum,ireff,imacc,ihacc,ihsoft,xyzmh_ptmass_label
 use units,      only:umass,udist,utime,unit_energ
 use physcon,    only:solarm,solarr,solarl
 use prompting,  only:prompt
 use dim,        only:nsinkproperties
 use io,         only:iprint
 integer :: i,j,iselect,ioption
 real    :: fac,var
 real,    intent(inout) :: xyzmh_ptmass(:,:)
 character(len=100)        :: dumpname

 do i = 1,nptmass
    print '("sink properties for #",i2," (in code units)")',i
    do j = 1,nsinkproperties
       write(iprint,"(3x,i2,1x,a,es10.3)")  j,xyzmh_ptmass_label(j),xyzmh_ptmass(j,i)
    enddo
 enddo
 if (nptmass == 1) then
    iselect = 1
 else
    iselect = 1
    call prompt('Select sink particle : ',iselect,1,nptmass)
    if (iselect < 1 .or. iselect > nptmass) stop 'wrong sink particle number'
 endif

 ioption =1
 do while (ioption > 0 .and. ioption < 17)
    call prompt('Select sink property (0 to exit): ',ioption,0,nsinkproperties)
    if (ioption == 0) exit
    var = xyzmh_ptmass(ioption,iselect)
    dumpname = '  o what value for ' // trim(xyzmh_ptmass_label(ioption)) // ' (in solar unit)'
    call prompt(dumpname,var)
    select case (ioption)
    case (ihacc,ihsoft,iReff)
       fac =  solarr / udist
    case (ilum)
       fac =  solarl * utime / unit_energ
    case (imacc,4)
       fac = solarm / umass
    case default
       fac = 1.
    end select
    xyzmh_ptmass(ioption,iselect) = var*fac
 enddo
 print *,'summary'
 do j = 1,nsinkproperties
    write(iprint,"(3x,i2,1x,a,es10.3)")  j,xyzmh_ptmass_label(j),xyzmh_ptmass(j,iselect)
 enddo

end subroutine set_sinkproperties

subroutine transform_from_corotating_to_inertial_frame(xyzh,vxyzu,npart,nptmass,omega_corotate,xyzmh_ptmass,vxyz_ptmass)
 use options,     only:iexternalforce
 use vectorutils, only:cross_product3D
 integer, intent(in) :: npart,nptmass
 real, intent(in) :: omega_corotate,xyzh(:,:),xyzmh_ptmass(:,:)
 real, intent(inout) :: vxyzu(:,:),vxyz_ptmass(:,:)
 real, dimension(3) :: omega_vec,omegacrossr
 integer :: i

 iexternalforce = 0
 omega_vec = (/ 0.,0.,omega_corotate /)
 do i=1,npart
    call cross_product3D(omega_vec,xyzh(1:3,i),omegacrossr)
    vxyzu(1,i) = vxyzu(1,i) + omegacrossr(1)
    vxyzu(2,i) = vxyzu(2,i) + omegacrossr(2)
    vxyzu(3,i) = vxyzu(3,i) + omegacrossr(3)
 enddo
 do i=1,nptmass
    call cross_product3D(omega_vec,xyzmh_ptmass(1:3,i),omegacrossr)
    vxyz_ptmass(1,i) = vxyz_ptmass(1,i) + omegacrossr(1)
    vxyz_ptmass(2,i) = vxyz_ptmass(2,i) + omegacrossr(2)
    vxyz_ptmass(3,i) = vxyz_ptmass(3,i) + omegacrossr(3)
 enddo

end subroutine transform_from_corotating_to_inertial_frame

subroutine set_trinary(mprimary,msecondary,mtertiary,semimajoraxis12,semimajoraxis13,&
                      accretion_radius1,accretion_radius2,accretion_radius3,&
                      xyzmh_ptmass,vxyz_ptmass,nptmass)
 real,    intent(in)    :: mprimary,msecondary,mtertiary
 real,    intent(in)    :: semimajoraxis12,semimajoraxis13
 real,    intent(in)    :: accretion_radius1,accretion_radius2,accretion_radius3
 real,    intent(inout) :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 integer, intent(inout) :: nptmass

 integer :: i1,i2,i3
 real    :: m1,m2,m3,mtot,dx12(3),dx13(3),dv12(3),dv13(3)
 real    :: x1(3),x2(3),x3(3),v1(3),v2(3),v3(3)

 i1 = nptmass + 1
 i2 = nptmass + 2
 i3 = nptmass + 3
 nptmass = nptmass + 3

 ! masses
 m1 = mprimary
 m2 = msecondary
 m3 = mtertiary
 mtot = m1 + m2 + m3

!
!--check for stupid parameter choices
!
! if (mprimary <= 0.)      stop 'ERROR: primary mass <= 0'
! if (massratio < 0.)      stop 'ERROR: binary mass ratio < 0'
! if (semimajoraxis <= 0.) stop 'ERROR: semi-major axis <= 0'
! if (eccentricity > 1. .or. eccentricity < 0.) &
!    stop 'ERROR: eccentricity must be between 0 and 1'

 dx12 = (/semimajoraxis12,0.,0./)
 dv12 = (/0.,sqrt((m1+m2)/dx12(1)),0./)

 dx13 = (/semimajoraxis13,0.,0./)
 dv13 = (/0.,sqrt(mtot/dx13(1)),0./)

 ! positions of each star so centre of mass is at zero
 x1 = -(dx12*m2 + dx13*m3)/mtot
 x2 = (dx12*m1 + dx12*m3 - dx13*m3)/mtot
 x3 = (dx13*m1 + dx13*m2 - dx12*m2)/mtot

 ! velocities
 v1 = -(dv12*m2 + dv13*m3)/mtot
 v2 = (dv12*m1 + dv12*m3 - dv13*m3)/mtot
 v3 = (dv13*m1 + dv13*m2 - dv12*m2)/mtot

!
!--positions and accretion radii
!
 xyzmh_ptmass(:,i1:i3) = 0.
 xyzmh_ptmass(1:3,i1) = x1
 xyzmh_ptmass(1:3,i2) = x2
 xyzmh_ptmass(1:3,i3) = x3
 xyzmh_ptmass(4,i1) = m1
 xyzmh_ptmass(4,i2) = m2
 xyzmh_ptmass(4,i3) = m3
 xyzmh_ptmass(5,i1) = accretion_radius1
 xyzmh_ptmass(5,i2) = accretion_radius2
 xyzmh_ptmass(5,i3) = accretion_radius3
 xyzmh_ptmass(6,i1) = 0.0
 xyzmh_ptmass(6,i2) = 0.0
 xyzmh_ptmass(6,i3) = 0.0
!
!--velocities
!
 vxyz_ptmass(:,i1) = v1
 vxyz_ptmass(:,i2) = v2
 vxyz_ptmass(:,i3) = v3

end subroutine set_trinary


end module moddump
