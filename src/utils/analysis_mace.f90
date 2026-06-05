!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module analysis
!
! Analysis routine to run MACE on phantom dumps
! 
! :References: None
!
! :Owner: Camille Landri
!
! :Runtime parameters: None
!
! :Dependencies: None
!
 use dvode_module
 use part,       only: maxp
 use raytracer,  only: get_all_tau_single
 use ftorch, only : torch_model, torch_tensor, torch_kCPU, torch_delete, &
                    torch_tensor_from_array, torch_model_load, torch_model_forward
 use omp_lib, only : omp_get_max_threads, omp_get_thread_num, omp_set_num_threads
 implicit none
 character(len=20), parameter, public :: analysistype = 'mace'
 public :: do_analysis

 integer, parameter :: abs_size = 468
 integer, parameter :: real_size = 472 ! 4 physical parameters + 468 abundances, input tensor size
 integer :: latent_size = 16
 character(len=8) :: epoch = '14'
 character(len=256) :: encoder_file, decoder_file

 real, allocatable    :: abundances(:,:), abundance_prev(:,:), one(:)
 character(len=16)    :: abundance_label(abs_size)
 character(len=8), dimension(30) :: saved_labels = ['H2','E','H','CO','HCN','N2','SIC2','CS','SiS','SiO','CH4','H2O', 'C2H4','HCL','NH3','H2S','C2H2','C2H','OH','CN','SIC','SIN','HC3N','C2','SINC','CH3','NH2','HS','CL','CO2']
 integer, allocatable :: saved_labels_i(:)
 integer, allocatable :: abundance_label_i(:)

 real, allocatable, target :: in_data(:,:), latent_abs(:, :), latent_abs_evolved(:, :), out_data(:, :)
 integer(8), allocatable :: iorig_old(:)
 integer, allocatable :: iprev(:)
 logical :: done_init = .false.
 real :: AuvAv = 4.65, albedo = 0.5

type(torch_model) :: encoder, decoder

 character(len=256) :: dir

 ! Minmaxes for scaling
 real, dimension(12) :: minmax
 real :: rho_min, rho_max
 real :: T_min, T_max
 real :: delta_min, delta_max
 real :: Av_min, Av_max
 real :: dt_max, dt_fract
 real :: n_min, n_max

 ! ODE parameters
 real           :: atol, rtol
 real, allocatable :: atol_arr(:), rtol_arr(:)
 real, allocatable :: C(:)
 real, allocatable :: A(:, :)
 real, allocatable :: B(:, :, :)
 character(len=16) :: label

 ! ODE solver variables
 type(dvode_t) :: solver
 integer  :: mf, istate,iopt, itol, itask, neq
 integer, parameter  :: lrw = 4096, liw = 46
 real :: t, tout, dt
 integer, dimension(liw)  :: iwork
 real, dimension(lrw) :: rwork
 real, allocatable :: y(:)

 ! Flag for testing
 logical :: test_pass
 logical :: verbose

 ! Timing variables
 real :: startTime, stopTime

 ! OpenMP variables
 integer :: thread, thread_count

 ! Model name
 character(len=256) :: model_file = '20260602_071233'


contains

subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
 use part,       only: isdead_or_accreted, iorig, rhoh, nptmass, xyzmh_ptmass, iReff, iboundary, igas, iphase, iamtype, maxp
 use linklist,   only: set_linklist
 use units,      only: utime,unit_density,udist
 use physcon,    only: atomic_mass_unit
 use eos,        only: get_temperature, ieos, gamma,gmw, init_eos
 use io,         only: fatal, iverbose

 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: num,npart,iunit
 real,             intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in) :: particlemass,time
 real, save    :: tprev = 0.
 integer, save :: nprev = 0
 type(torch_tensor), dimension(:), allocatable :: in_tensors, latent_tensors, latent_evolved_tensors, out_tensors
 real          :: dt_cgs, rho_cgs, numberdensity, T_gas, gammai, mui, AUV, xi
 real          :: column_density(npart), xyzh_copy(4,npart)
 real          :: max_radius, radius, time_count
 integer       :: i, j, k, i_radius, ierr, npart_copy = 0
 integer       :: iu=10,ios
 character(len=9) :: filename
 integer :: isize

 if (.not.done_init) then
    done_init = .true.
    print*, "Initialising MACE"
    ! read in model metadata and species labels
    call read_metadata(trim(model_file)//'/'//trim('meta.txt'), latent_size, dt_fract, epoch)
    allocate(saved_labels_i(size(saved_labels)))
    call read_abundance_labels(trim(model_file)//'/'//trim('abundance_label.txt'), abundance_label, saved_labels, saved_labels_i)
    ! Path to torchscript models
    encoder_file = trim(model_file)//'/'//trim(trim(adjustl(epoch))//'_encoder.pt')
    decoder_file = trim(model_file)//'/'//trim(trim(adjustl(epoch))//'_decoder.pt')
    ! load autoencoder
    print*, " - Loading Torch models"
    call torch_model_load(encoder, encoder_file, torch_kCPU)
    call torch_model_load(decoder, decoder_file, torch_kCPU)
    ! load minmax values
    call read_minmax(trim(model_file)//'/'//trim('minmax.txt'), 12, minmax)
    rho_min = log10(minmax(1))
    rho_max = log10(minmax(2))
    T_min = log10(minmax(3))
    T_max = log10(minmax(4))
    delta_min = log10(minmax(5))
    delta_max = log10(minmax(6))
    Av_min = log10(minmax(7))
    Av_max = log10(minmax(8))
    n_min = log10(minmax(9))
    n_max = log10(minmax(10))
    dt_max = minmax(11)

    ! Initialise abundances arrays
    print*, " - Allocating arrays"
    maxp = maxp
    allocate(abundances(abs_size,maxp))
    abundances = 0.
    allocate(abundance_prev(abs_size,maxp))
    abundance_prev = 0.

    ! Initialise other arrays
    allocate(one(maxp))
    one = 1.
    allocate(iorig_old(maxp))
    iorig_old = 0
    allocate(iprev(maxp))
    iprev = 0

    ! Allocate ODE solver arrays
    allocate(atol_arr(latent_size))
    allocate(rtol_arr(latent_size))
    allocate(C(latent_size))
    allocate(A(latent_size, latent_size))
    allocate(B(latent_size, latent_size, latent_size))
    allocate(y(latent_size))

    ! Initialise ODE parameters
    print*, " - Initialising ODE parameters"
    call read_ode_params(model_file, epoch, atol, rtol)
    neq    = latent_size ! number of first order odes.
    itol   = 2 !or 2 according as atol (below) is a scalar or array.
    itask  = 1 !for normal computation of output values of y at t = tout.
    istate = 1 !integer flag (input and output).  set istate = 1.
    iopt   = 0 !to indicate no optional input used.
    mf = 22
    atol_arr = atol
    rtol_arr = rtol

    print*, " - Setting abundances"
    !$omp parallel do default(none) &
    !$omp shared(npart,xyzh,abundances, abundance_label) &
    !$omp private(i)
    do i=1, npart
       if (.not.isdead_or_accreted(xyzh(4,i))) then
          call chem_init(abundances(:,i),abundance_label)
       endif
    enddo
    call init_eos(ieos, ierr)
    if (ierr /= 0) call fatal(analysistype, "Failed to initialise EOS")
    print*, "MACE initialisation complete"
 else
    dt_cgs = (time - tprev)*utime
    dt = dt_cgs / dt_max * dt_fract           ! scale dt for latent space evolution
    print*, " - Not first step data, timestep = ",dt_cgs, "npart = ",npart, "nprev = ",nprev
    thread_count = omp_get_max_threads()
    print*, " - Running on ", thread_count, " threads"
    call cpu_time(startTime)
    xyzmh_ptmass(iReff,1) = 2.
    npart_copy = npart
    xyzh_copy = xyzh(:,:npart)
    call set_linklist(npart_copy,npart_copy,xyzh_copy,vxyzu)
    ! temporary fix to get column density without companion (buggy when particle is aligned with the two stars, will need proper fix in get_all_tau_companion)
    call get_all_tau_single(npart, xyzmh_ptmass(1:3,1), xyzmh_ptmass(iReff,1), xyzh, one, xyzmh_ptmass(iReff,1), 5, .false., column_density)
    max_radius = 0.0
    do i = 1, npart
       if (.not.isdead_or_accreted(xyzh(4, i))) then
          radius = sqrt(xyzh(1, i)**2 + xyzh(2, i)**2 + xyzh(3, i)**2)
          if (radius > max_radius) then
             max_radius = radius
             i_radius = i
          endif
       endif
    enddo
    column_density = column_density + rhoh(xyzh(4,i_radius),particlemass)*unit_density * max_radius * udist
    call cpu_time(stopTime)
    print*, " - Column density calculation complete, time taken = ", (stopTime - startTime)/thread_count, " seconds"
    time_count = (stopTime - startTime)/thread_count
    call cpu_time(startTime)
    print*, "building data arrays for MACE"
    allocate(in_data(npart,real_size))
    allocate(latent_abs(npart,latent_size))
    allocate(latent_abs_evolved(npart,latent_size))
    allocate(out_data(npart,abs_size))

    !$omp parallel do default(none) &
    !$omp shared(npart,xyzh,vxyzu,dt_cgs,nprev,iorig,iorig_old,iprev,abundances,abundance_label,abundance_prev) &
    !$omp shared(particlemass,unit_density,ieos,gamma,gmw,column_density,albedo,AuvAv,in_data) &
    !$omp shared(rho_min, rho_max, T_min, T_max, delta_min, delta_max, Av_min, Av_max, n_min, n_max) &
    !$omp private(i,j,k,rho_cgs,numberdensity,T_gas,gammai,mui,AUV,xi, thread,radius)
    outer: do i=1,npart
    ! Loop over particles to build input arrays for MACE.
       if (.not.isdead_or_accreted(xyzh(4,i))) then
          inner: do j=1,nprev
             ! get previous index of particle
             if (iorig(i) == iorig_old(j)) then
                iprev(i) = j
                exit inner
             endif
          enddo inner
          ! Thermodynamic quantities
          rho_cgs = rhoh(xyzh(4,i),particlemass)*unit_density
          gammai = gamma
          mui    = gmw
          numberdensity = rho_cgs / (mui * atomic_mass_unit)
          T_gas = get_temperature(ieos,xyzh(1:3, i),rhoh(xyzh(4,i),particlemass),vxyzu(:,i),gammai,mui)
          T_gas = max(T_gas,20.0d0)
          radius = sqrt(xyzh(1, i)**2 + xyzh(2, i)**2 + xyzh(3, i)**2)
          ! Radiation quantities
          AUV = AuvAv * column_density(i) / (mui * atomic_mass_unit) / 1.87e21
          call get_xi(AUV, xi)
          in_data(i,1:4) = [numberdensity, T_gas, xi, AUV]
         if (iorig(i) == 5001) then
             print*, " - Particle ", iorig(i), " has input data ", in_data(i,1:4)
          end if
          if (j == iprev(i)) then
             thread = omp_get_thread_num()
             ! if particle existed in previous dump, get previous abundances
             in_data(i,5:real_size) = abundance_prev(:,iprev(i))
          else
             !print*, "Particle ", iorig(i), " not found in previous step, initializing abundances"
             call chem_init(abundances(:,i),abundance_label)
             in_data(i,5:real_size) = abundances(:,i)
          endif  
          ! Scale input data
          call log_and_scale(in_data(i,1), rho_min, rho_max) ! density
          call log_and_scale(in_data(i,2), T_min, T_max)     ! temperature
          call log_and_scale(in_data(i,3), delta_min, delta_max) ! delta
          call log_and_scale(in_data(i,4), Av_min, Av_max)   ! Av
          do k = 5, real_size
             call log_and_scale(in_data(i,k), n_min, n_max) ! abundances
          enddo
         if (iorig(i) == 5001) then
             print*, " - Particle ", iorig(i), " has input data scaled ", in_data(i,1:4)
          end if
       endif
    enddo outer
    !$omp end parallel do
    call cpu_time(stopTime)
    print*, " - Array building done, time taken = ", (stopTime - startTime)/thread_count, " seconds"
    time_count = time_count + (stopTime - startTime)/thread_count
   
    call cpu_time(startTime)
    write(*,*), "Creating Torch tensors"
    allocate(in_tensors(1))
    allocate(latent_tensors(1))
    allocate(latent_evolved_tensors(1))
    allocate(out_tensors(1))
    call torch_tensor_from_array(in_tensors(1), in_data, torch_kCPU)
    call torch_tensor_from_array(latent_tensors(1), latent_abs, torch_kCPU)
    call torch_tensor_from_array(latent_evolved_tensors(1), latent_abs_evolved, torch_kCPU)
    call torch_tensor_from_array(out_tensors(1), out_data, torch_kCPU)

    write(*,*), "Running encoder"
    call torch_model_forward(encoder, in_tensors, latent_tensors)
    call cpu_time(stopTime)
    print*, " - Encoder run complete, time taken = ", (stopTime - startTime)/thread_count, " seconds"
    time_count = time_count + (stopTime - startTime)/thread_count
    call cpu_time(startTime)
    !$omp parallel do default(none) &
    !$omp shared(npart,mf,iopt,itol,itask,neq,rtol_arr,atol_arr,latent_size,iorig,latent_abs, latent_abs_evolved, dt) &
    !$omp private(i,y,t,tout,solver,istate,iwork,rwork)
    do i=1,npart
       ! Evolve latent space
       ! Initialise ODE solver variables
       istate = 1
       y      = latent_abs(i,:) ! initial value of the dependent variable.
       t      = 0.0d0 ! initial value of the independent variable.
       tout   = dt ! first point where output is desired (/= t).
       call solver%initialize(f=ode)
       ! Call the ODE solver to evolve latent_abs
       call solver%solve(neq,y,t,tout,itol,rtol_arr,atol_arr,itask,istate,&
                         iopt,rwork,lrw,iwork,liw,mf)
       if (.not. istate == 2) then
          write(*,*) "Error in ODE solver for particle ", iorig(i), " istate = ", istate
       endif
       latent_abs_evolved(i,:) = y
       !print*, "Particle ", iorig(i), " evolved in latent space, y=", y
    end do
    !$omp end parallel do
    call cpu_time(stopTime)
    print*, " - Latent space evolution complete, time taken = ", (stopTime - startTime)/thread_count, " seconds"
    time_count = time_count + (stopTime - startTime)/thread_count
    call cpu_time(startTime)
   
    ! Decode
    write(*,*), "Running decoder"
    call torch_model_forward(decoder, latent_evolved_tensors, out_tensors)
    call cpu_time(stopTime)
    print*, " - Decoder run complete, time taken = ", (stopTime - startTime)/thread_count, " seconds"
    time_count = time_count + (stopTime - startTime)/thread_count
    call cpu_time(startTime)
    write(*,*), "Unscaling output and updating abundances array"

    !$omp parallel do default(none) &
    !$omp shared(npart,iorig,out_data,abundances,n_min, n_max,iphase) &
    !$omp private(i,j)
    do i=1,npart
       ! Unscale abundances
       if (iorig(i) == 5001) then
          print*, " - Scaled abundances for particle ", iorig(i)
          print*, "   > for H2: ", out_data(i,72)
          print*, "   > for HE: ", out_data(i,173)
          print*, "   > for CO: ", out_data(i,36)
       end if
       do j = 1, abs_size
          call unscale_and_unlog(out_data(i,j), n_min, n_max)
       end do
      if (iorig(i) == 5001) then
          print*, " - Unscaled abundances for particle ", iorig(i)
          print*, "   > for H2: ", out_data(i,72)
          print*, "   > for HE: ", out_data(i,173)
          print*, "   > for CO: ", out_data(i,37)
       end if
       abundances(:,i) = out_data(i,:)

    end do
    call cpu_time(stopTime)
    print *, 'Unscaling time, s : ',  (stopTime - startTime)/thread_count
    time_count = time_count + (stopTime - startTime)/thread_count
    ! clean up tensors
    deallocate(in_tensors)
    deallocate(latent_tensors)
    deallocate(latent_evolved_tensors)
    deallocate(out_tensors)
    deallocate(in_data)
    deallocate(latent_abs)
    deallocate(latent_abs_evolved)
    deallocate(out_data)
    print*, "done"
    print*,"Total MACE time for this dump, s: ", time_count
   endif

 ! store current step data before moving on to next step
 call write_chem(npart, dumpfile, saved_labels, saved_labels_i)
 nprev = npart
 tprev = time
 iorig_old = iorig
 abundance_prev = abundances
 end subroutine do_analysis

 subroutine get_xi(AUV, xi)
   use physcon, only: pi
   real, intent(in) :: AUV
   real, intent(out) :: xi
   real :: W(6), GA(6), ceta
   integer :: i

   W(1) = 0.17132449
   W(2) = 0.36076157
   W(3) = 0.46791393
   W(4) = W(1)
   W(5) = W(2)
   W(6) = W(3)
   GA(1) = 0.93246951
   GA(2) = 0.66120939
   GA(3) = 0.23861919
   GA(4) = -GA(1)
   GA(5) = -GA(2)
   GA(6) = -GA(3)

   xi = 0.0
   do i=1,6
      ceta = (pi*GA(i)+pi)/2.0
      xi=xi+(W(i)*(sin(ceta)*exp((-AUV*ceta)/sin(ceta))))
   enddo
   xi = (pi/4.0)*xi
 
   end subroutine get_xi

 subroutine write_chem(npart, dumpfile, saved_labels, saved_labels_i)
    integer, intent(in)          :: npart
    character(len=*), intent(in) :: dumpfile
    character(len=*), dimension(:), intent(in) :: saved_labels
    integer, dimension(:), intent(in) :: saved_labels_i
    character(len=256) :: header
    integer :: iu, isize, i, k
       write(*,*) " - Writing chemical abundances to file"
       open(newunit=iu,file=dumpfile//'.comp',status='replace',action='write')
       header = '# '
       do k = 1, size(saved_labels)
          header = trim(header) // trim(saved_labels(k)) // ', '
       end do
       write(iu, *) header
       do i=1, npart
          do k = 1, size(saved_labels)
             write(iu, '(ES14.7,1x)', advance="no") abundances(saved_labels_i(k),i)
          end do
          write(iu,*) ! new line
       enddo
       close(iu)
   end subroutine write_chem

 subroutine ode(me, neq, t, y, ydot)
      ! Defines the ODE system dyi/dt = Ci + Aij*yj + Bijk*yi*yk
      class(dvode_t),intent(inout) :: me
      integer :: neq
      real :: t
      real :: y(neq)
      real :: ydot(neq)
      integer :: i, j, k

      do j = 1, neq
         ydot(j) = C(j)
         do i = 1, neq
               ydot(j) = ydot(j) + A(j, i) * y(i)
         end do
         do i = 1, neq
               do k = 1, neq
                  ydot(j) = ydot(j) + B(j, i, k) * y(i) * y(k)
               end do
         end do
      end do
   end subroutine ode

 subroutine read_ode_params(model_file, epoch, atol, rtol)
      ! Subroutine to read ODE parameters from files
      implicit none
      character(len=*), intent(in), optional :: model_file, epoch
      real, intent(out) :: atol, rtol
      integer :: i, j, k

      ! Load ODE parameters from file
      print*, " - Reading ODE parameters from " // trim(trim(adjustl(epoch)) //'_ODE_params.txt')
      open(unit=10, file=trim(model_file)//'/'//trim(trim(adjustl(epoch))//'_ODE_params.txt'), status='old')
      read(10, *) label, atol
      read(10, *) label, rtol
      close(10)
      write(*,*) "      > atol =", atol
      write(*,*) "      > rtol =", rtol

      ! Read C (vector)
      print*, " - Reading Coefficient C from " // trim(trim(adjustl(epoch)) //'_ODE_C.txt')
      open(unit=11, file=trim(model_file)//'/'//trim(trim(adjustl(epoch))//'_ODE_C.txt'), status='old')
      do i = 1, latent_size
         read(11, *) C(i)
      end do
      close(11)

      ! Read A (matrix)
      print*, " - Reading Coefficient A from " // trim(trim(adjustl(epoch)) //'_ODE_A.txt')
      open(unit=12, file=trim(model_file)//'/'//trim(trim(adjustl(epoch))//'_ODE_A.txt'), status='old')
      do i = 1, latent_size
         do j = 1, latent_size
               read(12, *) A(i, j)
         end do
      end do
      close(12)

      ! Read B (tensor)
      print*, " - Reading Coefficient B from " // trim(trim(adjustl(epoch)) //'_ODE_B.txt')
      open(unit=13, file=trim(model_file)//'/'//trim(trim(adjustl(epoch))//'_ODE_B.txt'), status='old')
      do i = 1, latent_size
         do j = 1, latent_size
               do k = 1, latent_size
                  read(13, *) B(i, j, k)
               end do
         end do
      end do
      close(13)
   end subroutine read_ode_params

 subroutine log_and_scale (x, xmin, xmax)
      ! scale parameters
      ! First clip to min
      ! log transform then normalise
      ! xmin and xmax are log10 of min and max values for parameter

      real, intent(inout) :: x
      real, intent(in) :: xmin, xmax
      ! clip to min
      if (x < 10.0d0**xmin) then
         x = 10.0d0**xmin
      end if
      x = log10(x)
      x = (x - xmin) / ABS(xmax - xmin)
   end subroutine log_and_scale


 subroutine unscale_and_unlog (x, xmin, xmax)
      ! unscale parameters 
      ! denormalise then exp10 transform
      ! xmin and xmax are log10 of min and max values for parameter
      real, intent(inout) :: x
      real, intent(in) :: xmin, xmax
      x = x * ABS(xmax - xmin) + xmin
      x = 10.0d0**x
   end subroutine unscale_and_unlog

 subroutine read_minmax(filename, filesize, minmax)
      ! Subroutine to read minmax values from a file
      implicit none
      character(len=*), intent(in) :: filename
      character(len=128) :: dummy
      integer, intent(in) :: filesize
      real, dimension(filesize), intent(out) :: minmax
      integer :: i, unit_num,ios

      write(*,*) " - Reading minmax values from ", filename
      unit_num = 20
      open(unit=unit_num, file=filename, status='old')
         do i = 1, filesize
            read(unit_num, *, iostat=ios) dummy, minmax(i)
            if (ios /= 0) then
               write(*,*) "Error reading minmax file."
               stop 1
            end if
         end do
      close(unit_num)
   end subroutine read_minmax

 subroutine read_abundance_labels(filename, labels, saved_labels, saved_labels_i)
      ! Subroutine to read abundances labels from a file
      implicit none
      character(len=*), intent(in) :: filename
      character(len=*), dimension(:), intent(out) :: labels
      character(len=8), dimension(:), intent(in) :: saved_labels
      integer, dimension(:), intent(out) :: saved_labels_i
      integer :: i, j, unit_num, ios
      write(*,*) " - Reading ", size(labels), "abundance_label labels from ", filename
      unit_num = 30
      open(unit=unit_num, file=filename, status='old')
      do i = 1, size(labels)
         read(unit_num, *, iostat=ios) labels(i)
         if (ios /= 0) then
            write(*,*) "Error reading abundance_label file."
            stop 1
         end if
         do j = 1, size(saved_labels)
            if (trim(labels(i)) == trim(saved_labels(j))) then
               saved_labels_i(j) = i
               exit
            end if
         end do
      end do
      close(unit_num)
   end subroutine read_abundance_labels

 subroutine read_metadata(filename, latent_size, dt_fract, epoch)
      ! Subroutine to read metadata from a file
      implicit none
      character(len=*), intent(in) :: filename
      integer, intent(inout) :: latent_size
      character(len=8), intent(inout) :: epoch
      real, intent(inout) :: dt_fract
      character(len=128) :: label, metadata
      integer :: i, unit_num, ios, metadata_size
      write(*,*) " - Reading metadata from ", filename
      unit_num = 40
      open(unit=unit_num, file=filename, status='old')
      ! get length of metadata
      metadata_size = 0
      do
         read(unit_num, *, iostat=ios)
         if (ios /= 0) exit
         metadata_size = metadata_size + 1
      end do
      rewind(unit_num)
      ! read metadata
      do i = 1, metadata_size
         read(unit_num, *, iostat=ios) label, metadata
         if (trim(label) == 'z_dim:') then
            read(metadata, *) latent_size
            write(*,*) "      > latent_size =", latent_size
         else if (trim(label) == 'dt_fract:') then
            read(metadata, *) dt_fract
            write(*,*) "      > dt_fract =", dt_fract
         else if (trim(label) == 'epoch:') then
            read (metadata, *) epoch
            write(*,*) "      > epoch =", epoch
         end if

         if (ios /= 0) then
            write(*,*) "Error reading metadata file."
            stop 1
         end if
      end do
      close(unit_num)
   end subroutine read_metadata

 subroutine chem_init(abundance, abundance_label)
    ! Subroutine to initialise chemical abundance
    implicit none
    character(len=16), dimension(:), intent(in) :: abundance_label
    real, dimension(size(abundance_label)), intent(out) :: abundance
    integer :: i
    ! Initial abundance for the krome model taken from Agúndez et al. (2020)
    ! H2, He, CO, C2H2, HCN, N2, SiC2, CS, SiS, SiO, CH4, H2O, HCl, C2H4, NH3, HCP, HF, H2S, E
    do i = 1, size(abundance_label)
       if (abundance_label(i) == 'H2') then
          abundance(i) = 0.5d0
       else if (abundance_label(i) == 'HE') then
          abundance(i) = 8.5d-2
       else if (abundance_label(i) == 'CO') then
          abundance(i) = 4.0d-4
       else if (abundance_label(i) == 'C2H2') then
          abundance(i) = 2.19d-5
       else if (abundance_label(i) == 'HCN') then
          abundance(i) = 2.045d-5
       else if (abundance_label(i) == 'N2') then
          abundance(i) = 2.0d-5
       else if (abundance_label(i) == 'SIC2') then
          abundance(i) = 9.35d-6
       else if (abundance_label(i) == 'CS') then
          abundance(i) = 5.3d-6
       else if (abundance_label(i) == 'SIS') then
          abundance(i) = 2.99d-6
       else if (abundance_label(i) == 'SIO') then
          abundance(i) = 2.51d-6
       else if (abundance_label(i) == 'CH4') then
          abundance(i) = 1.75d-6
       else if (abundance_label(i) == 'H2O') then
          abundance(i) = 1.275d-6
       else if (abundance_label(i) == 'HCL') then
          abundance(i) = 1.625d-7
       else if (abundance_label(i) == 'C2H4') then
          abundance(i) = 3.425d-8
       else if (abundance_label(i) == 'NH3') then
          abundance(i) = 3.0d-8
       else if (abundance_label(i) == 'HCP') then
          abundance(i) = 1.25d-8
       else if (abundance_label(i) == 'HF') then
          abundance(i) = 8.5d-9
       else if (abundance_label(i) == 'H2S') then
          abundance(i) = 2.0d-9
       else
          abundance(i) = 0.0
       end if
    end do
   end subroutine chem_init

end module analysis