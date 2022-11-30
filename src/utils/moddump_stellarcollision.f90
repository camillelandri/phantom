!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2022 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module moddump
!
! test common envelope - put point source star next to gas sphere
!
! :References: None
!
! :Owner: Mike Lau
!
! :Runtime parameters: None
!
! :Dependencies: centreofmass, dim, extern_corotate, externalforces,
!   infile_utils, io, options, part, physcon, prompting, readwrite_dumps,
!   rho_profile, setbinary, table_utils, timestep, units
!
 implicit none
 integer, parameter :: prec=kind(1d0), i64=selected_int_kind(15)
contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use part,              only:nptmass,xyzmh_ptmass,vxyz_ptmass,ihacc,ihsoft,igas,&
                             delete_dead_or_accreted_particles,mhd,rhoh,shuffle_part,&
                             kill_particle,copy_particle,igas
 use setbinary,         only:set_binary
 use units,             only:umass,udist,utime
 use physcon,           only:au,solarm,solarr,gg,pi
 use centreofmass,      only:reset_centreofmass,get_centreofmass
 use prompting,         only:prompt
 use externalforces,    only:iext_corotate
 use infile_utils,      only:open_db_from_file,inopts,read_inopt,close_db
 use table_utils,       only:yinterp
 use rho_profile,       only:read_mesa
 use dim,               only:maxptmass
 use io,                only:fatal,idisk1,iprint
 use readwrite_dumps,   only:read_dump
 use timestep,          only:tmax

 integer, intent(inout)    :: npart
 integer, intent(inout)    :: npartoftype(:)
 real,    intent(inout)    :: massoftype(:)
 real,    intent(inout)    :: xyzh(:,:),vxyzu(:,:)
 integer                   :: i,ierr
 integer                   :: nstar1,nstar2
 real                      :: time2,hfact2
 real                      :: x1com(3), v1com(3), x2com(3), v2com(3), xcom(3), vcom(3)
 real                      :: m, m1,m2,rad1,rad2, mpart1, mpart2, halfrad1,halfrad2
 real                      :: a, b, p, lambda, Vinf, Vcontact, Vratio
 real                      :: r, vr, theta, theta0, vtheta, e
 real                      :: dx, dy, vx, vy
 character(len=100)        :: dumpname


!takes necessary inputs from user 1
print*, 'Current length unit is ', udist ,'cm):'
print*, 'Current mass unit is ', umass,'g):'
print*, 'Current time unit is ', utime,'s):'
mpart1 = massoftype(igas)
call reset_centreofmass(npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)

!removes the dead or accreted particles for a correct total mass computation
call delete_dead_or_accreted_particles(npart,npartoftype)
print*,' Got ',npart,npartoftype(igas),' after deleting accreted particles'

! remove point masses just in case
nptmass = 0

   nstar1 = npart ! save npart in star 1
   dumpname = ''
   call prompt('Enter name of second dumpfile',dumpname)
   nstar2 = nstar1
   call prompt('Enter no. of particles in second dumpfile',nstar2)
   
   ! Get mass of star 1 before adding star2 otherwise it gets screwed up
   call get_centreofmass(x1com, v1com, nstar1, xyzh(:,1:nstar1), vxyzu(:,1:nstar1), mass=m1)

   ! Move star 1 particles to avoid getting overwritten when reading second dump file.
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
   if (ierr /= 0) stop 'error reading second dump file'

   if (nstar1 > nstar2) then ! Move ith particle of star 1 to nstar2+i
      do i=1,nstar1
         call copy_particle(nstar1+i,nstar2+i,.false.)
      enddo
   endif

   npart = nstar1 + nstar2
   npartoftype(igas) = npart
   mpart2 = massoftype(igas)  

! calculate radii of each star
   call get_centreofmass(x2com, v2com, nstar2, xyzh(:,1:nstar2), vxyzu(:,1:nstar2), mass=m2)
   call get_radii(npart,xyzh,nstar1,nstar2,x1com,x2com,rad1,rad2,1.,mpart1,mpart2,m1,m2)
   call get_radii(npart,xyzh,nstar1,nstar2,x1com,x2com,halfrad1,halfrad2,.5,mpart1,mpart2,m1,m2)
   
   print*, 'Star 1 has mass ', m1 ,', radius ', rad1, ' and half-mass radius ', halfrad1
   print*, 'Star 2 has mass ', m2 ,',radius ', rad2, ' and half-mass radius ', halfrad2

   !call prompt('Enter radius of star 1:', rad1)
   !call prompt('Enter radius of star 2:', rad2)

   !call prompt('Enter half mass radius of star 1:', halfrad1)
   !call prompt('Enter half mass radius of star 2:', halfrad2)


! calculate contact velocity and ask for relative velocity (G=1 here)
   Vcontact = sqrt((2*(m1+m2))/(rad1+rad2))
   Vratio = 0.07
   print*,'Contact velocity is ', Vcontact
   call prompt('Enter ratio of relative infinity velocity over contact velocity:', Vratio)
   Vinf = Vratio*Vcontact
   call prompt('Enter relative infinity velocity:', Vinf)

   
! ask for periastron distance and deduce impact parameter
   lambda = 0.34
   m = m1 +m2
   call prompt('Enter ratio of periastron distance over sum of the radii:', lambda)
   b =(rad1 + rad2) * sqrt(lambda**2 + lambda/(Vratio**2))
   call prompt('Enter impact parameter:', b)

! Get start postition of star 1 in rest frame of star 2   
! Hyperbolic trajectory parameters - De Vittori et al. 2012
   a = - m / Vinf**2
   e = sqrt(1 + (b/a)**2)
   theta0 = acos(-1/e)
   r = 3*(rad1 + rad2)
   p = abs(a*(e**2-1))
   call prompt('Enter distance between stars at t0:', r)
   theta = pi - theta0 + acos(p/(r*e)-1/e)
   print*, 'theta0=', theta0, 'theta=', theta
   dx = r*cos(theta)
   dy = r*sin(theta)
   print*, 'star 1 starts at x=', dx, ', y=',dy,' in the rest frame of star 2.'

! Get start velocity of star 1
   vr = -sqrt(Vinf**2 + 2*m/r - b**2*Vinf**2/r**2)
   vtheta = -b*Vinf/r
   vx = (vr*cos(theta) - vtheta*sin(theta))
   vy = vr*sin(theta) + vtheta*cos(theta)

   print*, 'star 1 starts with velocity vx=', vx, ', vy=', vy,' in the rest frame of star 2.'
   print*, 'Orbit parameters : a=', a,', b=', b,', e=', e, ', vr=', vr,', vtheta=', vtheta, ', theta=', theta
   ! shift star1 to position in star 2 CoM
   do i=nstar2+1,npart
      xyzh(1,i) = xyzh(1,i) + dx
      xyzh(2,i) = xyzh(2,i) + dy
      vxyzu(1,i) = vxyzu(1,i) + vx
      vxyzu(2,i) = vxyzu(2,i) + vy
   enddo
   ! Shift point mass 1 position in CoM of star 2
   xyzmh_ptmass(1,1) = xyzmh_ptmass(1,1) + dx
   xyzmh_ptmass(2,1) = xyzmh_ptmass(2,1) + dy
   vxyz_ptmass(1,1) = vxyz_ptmass(1,1) + vx
   vxyz_ptmass(2,1) = vxyz_ptmass(2,1) + vy


! Go back to CoM frame
   ! Get CoM
   call get_centreofmass(xcom, vcom, npart, xyzh(:,:), vxyzu(:,:))
   print*, 'CoM at ',  xcom(1), xcom(2), xcom(3)
   print*, 'vCoM at ',  vcom(1), vcom(2), vcom(3)
   ! Shift pt masses to CoM
   xyzmh_ptmass(1:3,1) = xyzmh_ptmass(1:3,1) - xcom(1:3)
   xyzmh_ptmass(1:3,2) = xyzmh_ptmass(1:3,2) - xcom(1:3)
   vxyz_ptmass(1:3,1) = vxyz_ptmass(1:3,1) - vcom(1:3)
   vxyz_ptmass(1:3,2) = vxyz_ptmass(1:3,2) - vcom(1:3)
   ! Shift stars to CoM
   ! Star 1
   do i=nstar2+1,npart
      xyzh(1:3,i) = xyzh(1:3,i) - xcom(1:3)
      vxyzu(1:3,i) = vxyzu(1:3,i) - vcom(1:3)
   enddo
   ! Star 2
   do i=1,nstar2
      xyzh(1:3,i) = xyzh(1:3,i) - xcom(1:3)
      vxyzu(1:3,i) = vxyzu(1:3,i) - vcom(1:3)
   enddo

print*, 'star 1 starts at x=',  xyzmh_ptmass(1,1), ', y=', xyzmh_ptmass(2,1),' in the CoM.'
print*, 'star 1 starts with velocity vx=', vxyz_ptmass(1,1), ', vy=', vxyz_ptmass(2,1),' in the CoM.'
print*, 'star 2 starts at x=',  xyzmh_ptmass(1,2), ', y=', xyzmh_ptmass(2,2),' in the CoM.'
print*, 'star 2 starts with velocity vx=', vxyz_ptmass(1,2), ', vy=', vxyz_ptmass(2,2),' in the CoM.'
 return
end subroutine modify_dump

!
! Determine radius of each star based on particles - modified from moddump_binarystar
!
subroutine get_radii(npart,xyzh,nstar1,nstar2,x1com,x2com,rad1,rad2,percent,mpart1,mpart2,m1,m2)
   integer, intent(in)    :: npart
   real,    intent(inout) :: xyzh(:,:)
   integer, intent(in)    :: nstar1, nstar2
   real,    intent(in)    :: x1com(:),x2com(:)
   real,    intent(out)   :: rad1, rad2
   real      :: rads1(nstar1), rads2(nstar2)
   integer :: i
   real :: dx, dy, dz, dr, percent
   real :: mpart1, mpart2, m1, m2

   rad2 = 0.0
   do i = 1, nstar2
      dx   = xyzh(1,i) - x2com(1)
      dy   = xyzh(2,i) - x2com(2)
      dz   = xyzh(3,i) - x2com(3)
      dr   = sqrt(dx*dx + dy*dy + dz*dz)
      rads2(i) = dr
   enddo
   call dpquicksort(rads2)
   i = 0
   do while (i*mpart2 < percent*m2)
      i = i+1
   enddo
   rad2 = max(rads2(i),rads2(i-1))

   rad1 = 0.0
   do i = nstar2+1, npart
      dx   = xyzh(1,i) - x1com(1)
      dy   = xyzh(2,i) - x1com(2)
      dz   = xyzh(3,i) - x1com(3)
      dr   = sqrt(dx*dx + dy*dy + dz*dz)
      rads1(i-nstar2) = dr
   enddo
   call dpquicksort(rads1)
   i = 0
   do while (i*mpart1 < percent*m1)
     i = i+1
   enddo
   rad1 = max(rads1(i),rads1(i-1))
  
  end subroutine get_radii
  

 ! dual pivot quicksort - from https://www.mjr19.org.uk/
  recursive subroutine dpquicksort(array)
    real(prec), intent(inout)::array(:)
    real(prec) :: temp,p1,p2
    integer :: i,j,last,l,k,g

    last=size(array)

    if (last.lt.40) then ! use insertion sort on small arrays
       do i=2,last
          temp=array(i)
          do j=i-1,1,-1
             if (array(j).le.temp) exit
             array(j+1)=array(j)
          enddo
          array(j+1)=temp
       enddo
       return
    endif
    p1=array(last/3)
    p2=array(2*last/3)
    if (p2.lt.p1) then
       temp=p1
       p1=p2
       p2=temp
    endif
    array(last/3)=array(1)
    array(1)=p1
    array(2*last/3)=array(last)
    array(last)=p2

    g=last
    l=2
    do while (array(l).lt.p1)
       l=l+1
    enddo
    k=l

    do while(k.lt.g)
       temp=array(k)
       if (temp.lt.p1) then
          array(k)=array(l)
          array(l)=temp
          l=l+1
       else if (temp.gt.p2) then
          do while(array(g-1).gt.p2)
             g=g-1
          enddo
          if (k.ge.g) exit
          g=g-1
          if (array(g).lt.p1) then
             array(k)=array(l)
             array(l)=array(g)
             array(g)=temp
             l=l+1
          else
             array(k)=array(g)
             array(g)=temp
          endif
       endif
       k=k+1
    enddo
    if (l.gt.2) then
       array(1)=array(l-1)
       array(l-1)=p1
       call dpquicksort(array(1:l-2))
    endif
    call dpquicksort(array(l:g-1))
    if (g.lt.last) then
       array(last)=array(g)
       array(g)=p2
       call dpquicksort(array(g+1:last))
    endif

  end subroutine dpquicksort



end module moddump