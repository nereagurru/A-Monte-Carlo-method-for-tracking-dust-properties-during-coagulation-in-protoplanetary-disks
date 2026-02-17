
! this module performes collisions between representative bodies
! in this module all the quantities should be computed at the center of current cell
! in the matrices: first index -> representative particle, 2nd -> physical
module collisions

   use constants,  only: pi, third, mH2, AH2, AU, year
   use discstruct, only: cs, omegaK, densg, Pg, vgas, alpha
   use grid,       only: g
   use initproblem,only: m0
   use parameters, only: dmmax, vfrag 
#ifdef MULTI_COMPONENT
   use parameters, only: matdenssi, matdensw, tolerance
#else
   use parameters, only: matdens, con1, con2, tolerance
#endif
#ifdef EROSION
   use parameters, only: erosion_mass_ratio
#endif
   use types,      only: swarm, list_of_swarms
   use hdf5
   use hdf5output, only: hdf5_file_write, hdf5_file_t
#ifdef KERNEL_TEST
use parameters, only: repeat, path, datadir
#endif
   implicit none

   private
#ifdef KERNEL_TEST
   public :: mc_collisions_test
#else
   public :: mc_collisions
#endif

#ifdef TESTCOLLISIONS
   public :: col_rates, rel_vels, vel_vs_centr, vel_rd_centr, stokes_nr_centr
#endif

   contains

#ifdef KERNEL_TEST
#else
   ! the routine performes collisional evolution on swarms located in the cell nr,ni with the aid of the MC algorithm
   subroutine mc_collisions(nr, ni, bin, swrm, dtime, realtime, kk, first_idx)
      implicit none
      type(swarm), dimension(:), allocatable, target            :: swrm
      type(list_of_swarms), dimension(:,:), allocatable, target :: bin
      type(swarm), dimension(:), pointer              :: swarms      ! local rps array
      integer, intent(in)                             :: nr, ni      ! indices of the cell
      real, intent(in)                                :: dtime       ! time step
      real, intent(in)                                :: realtime    ! physical time
      integer, intent(in)                             :: first_idx   ! first index of swrm in rbin
      integer                                         :: nsws        ! number of representative particles in given cell
      integer                                         :: nri, nrk    ! indices of physical and representative particles choosen to the next collision
      real, dimension(:), allocatable                 :: colrates    ! collision rates matrix
      real, dimension(:), allocatable                 :: colrates_rp
      real, dimension(:), allocatable                 :: relvels     ! relative velocities matrix
      real, dimension(:), allocatable                 :: accelncol   ! coagulation acceleration matrix (in the case of high mass ratio,
                                                                     ! instead of performing every collision separately, we group the collisions)
      real, dimension(:), allocatable                 :: stokesnr    ! stokes numbers of particles
      real, dimension(:), allocatable                 :: vs, vr      ! vertical settilng and radial drift velocities
      real                                            :: local_time, dt ! physical time spent in the cell, time step between two collisions
      real                                            :: rand        ! random number
      real                                            :: totr        ! total collision rate
      real                                            :: vn          ! maximum radial velocity from the pressure gradient
      real                                            :: Reynolds, v0, Vg2, veta, tL, teta ! relative velocities stuff
      real                                            :: lmfp, gasdens ! mean free path, gas density
      integer                                         :: i, k, j, l, ij, Ncol, j_loop, idnr, imax
      integer, dimension(:,:), allocatable            :: ij_back
      integer, intent(out)                            :: kk ! collisions counter
      real                                            :: Kepler_freq, cs_speed
      real :: deltar, deltaz   ! radial and vertical size of grid (for adaptative dmmax)
      real :: mass_variation = 0. !0.01
      real :: weightw, weightk
      real :: r, mmax
      real :: mean_mswarm
      integer, allocatable :: idx(:)
      integer :: w
!#ifdef MULTI_COMPONENT
      real:: rpos, zpos ! this is needed to avoid particles ending up in the same position during merging
!#endif

      deltar = g%rup(nr) - g%rlo(nr)
      deltaz = g%zup(nr,ni) - g%zlo(nr,ni)


      ! calculation of some values needed for calculations
      gasdens = densg(g%rce(nr),g%zce(nr,ni),realtime)         ! gas density in the center of cell
      lmfp = mH2 / ( gasdens * AH2 )                           ! mean free path in gas in the center of cell
      Kepler_freq = omegaK(g%rce(nr))                          ! keplerian frequency at the radial centre of cell
      cs_speed = cs(g%rce(nr))
      Reynolds = sqrt(0.5 * pi) * alpha(g%rce(nr)) * cs_speed / (Kepler_freq * lmfp) ! Reynolds number
      v0 = sqrt(alpha(g%rce(nr))) * cs_speed                                 ! velocity of the largest eddy
      veta = v0 * Reynolds**(-0.25)                                               ! velocity of the smallest eddy
      tL = 1. / Kepler_freq
      teta = Reynolds**(-0.5) * tL                                                ! overturn time of the smallest eddy
      vn = 0.25 * (Pg(g%rce(nr)+1., g%zce(nr,ni),realtime) - Pg(g%rce(nr)-1., g%zce(nr,ni), realtime)) / &  ! maximum radial velocity from the pressure gradient
               gasdens / Kepler_freq
      Vg2 = 1.5 *v0**2.0                                                          ! turbulent velocity of gas squared

      !-------------------------------------------------------------------------------------

      swarms => bin(nr,ni)%p       ! points to the swarms array that contains only swarms located in the current cell nr,ni
      swarms(:)%coll_f = 0         ! setting collision flag to zero, will be updated to 1 if collisions happen
      nsws = size(swarms)          ! number of swarms in the current cell
#ifdef MULTI_COMPONENT
      mean_mswarm = sum(swarms(:)%mswarm*(1.+swarms(:)%w))/nsws
#else
      mean_mswarm = sum(swarms(:)%mswarm)/nsws
#endif
      Ncol = nsws*(nsws+1)/2       ! number of all potential collisions
      allocate (colrates(Ncol), accelncol(Ncol), relvels(Ncol), ij_back(Ncol, 2))
      allocate (stokesnr(nsws), vs(nsws), vr(nsws), colrates_rp(nsws))

      do i=1, nsws
         do j=i, nsws
            ij = (i-1)*nsws - (i-1)*(i-2)/2 + (j - i + 1)!(i-1)*nsws -(i-1)*i/2 + (j-i)
            ij_back(ij,1) = i
            ij_back(ij,2) = j 
         enddo
      enddo
      ! calculates initial Stokes number for particles and their settling and radial drift velocities
      do i = 1, nsws
         call stokes_nr_centr(i, swarms, stokesnr, lmfp, gasdens, Kepler_freq, cs_speed)
         call vel_vs_centr(nr, ni, i, stokesnr, Kepler_freq, vs)
         call vel_rd_centr(nr, i, stokesnr, vr, vn, realtime)
      enddo

      colrates_rp(:) = 0.

      ! calculates relative velocities and collision rates matrix
      do ij = 1, Ncol
         i = ij_back(ij,1)
         j = ij_back(ij,2)
         call rel_vels(i, j, swarms, stokesnr, vr, vs, relvels(ij), vn, Reynolds, veta, &
                     & Vg2, tL, teta, Kepler_freq, cs_speed)
         call col_rates(i, j, swarms, relvels(ij), colrates(ij), accelncol(ij), g%vol(nr,ni), deltar, deltaz)
         colrates_rp(i) = colrates_rp(i) + colrates(ij)
      enddo

      !------------ MAIN COLLISIONS LOOP ----------------------------------------------------
      local_time = 0.0
      kk = 0
      do while (local_time < dtime)
         totr = sum(colrates_rp)                 ! total collision rate
         call random_number(rand)

         dt = -1. * log(rand) / totr       ! time step between this and the next collision
         !write(*,*) 'dt vs. dtime ', dt/year, dtime/year, maxval(colrates), minval(colrates)

         if (dt > dtime) then ! 0 or 1 collisions, decided by a random number
            call random_number(rand)
            if (rand> dtime/dt) then
               local_time = dtime
               cycle
            endif
         endif
         local_time = local_time + dt      ! update of the local time
         call choose_swarms(nri, nrk, swarms, colrates_rp, colrates, totr)
         ! nri saves the one with lowest npar, so the one that gains mass during sticking
         call collision(nri, nrk, swarms, relvels, accelncol)

         ! if one of them got emptied, refill it; nri cannot get emptied
#ifdef MULTI_COMPONENT
         if (swarms(nrk)%mswarm*(1.+swarms(nrk)%w) < mean_mswarm*tolerance) then
#else
         if (swarms(nrk)%mswarm < mean_mswarm*tolerance) then
#endif
            if (swarms(nrk)%mswarm > 0.) then ! if lower, simply empty

               ! find most similar particle to nrk
               w = minloc(abs(swarms(:)%mass-swarms(nrk)%mass), dim=1, mask=[(i/= nrk, i = 1, size(swarms))])
               ! merge (not collide, compute average)
#ifdef MULTI_COMPONENT
               weightw = 1./(1.+swarms(nrk)%mswarm/swarms(w)%mswarm*(1.+swarms(nrk)%w)/(1.+swarms(w)%w))
               weightk = 1./(1.+swarms(w)%mswarm/swarms(nrk)%mswarm*(1.+swarms(w)%w)/(1.+swarms(nrk)%w))
#else
               weightw = 1./(1.+swarms(nrk)%mswarm/swarms(w)%mswarm)
               weightk = 1./(1.+swarms(w)%mswarm/swarms(nrk)%mswarm)
#endif
               swarms(w)%npar = swarms(w)%npar + swarms(nrk)%npar
#ifdef MULTI_COMPONENT
               ! multicomponent is 0D so we do not update position
               swarms(w)%fw = swarms(w)%fw*weightw + swarms(nrk)%fw*weightk
               swarms(w)%w = swarms(w)%fw/(1.-swarms(w)%fw)
               swarms(w)%mswarm = swarms(w)%mswarm+swarms(nrk)%mswarm 
               swarms(w)%mass = swarms(w)%mswarm*(1.+swarms(w)%w)/swarms(w)%npar
#else
               !swarms(w)%rdis = swarms(w)%rdis*weightw + swarms(nrk)%rdis*weightk
               !swarms(w)%zdis = swarms(w)%zdis*weightw + swarms(nrk)%zdis*weightk
               swarms(w)%mswarm = swarms(w)%mswarm+swarms(nrk)%mswarm ! this is not best practise, because they are large numbers
               swarms(w)%mass = swarms(w)%mswarm/swarms(w)%npar

               
#endif
               
               swarms(w)%coll_f =  1 ! otherwise it won't get updated
               swarms(nrk)%mswarm = 0.
               swarms(nrk)%npar = 0.
               ! update w group
               call stokes_nr_centr(w, swarms, stokesnr, lmfp, gasdens, Kepler_freq, cs_speed)
               call vel_vs_centr(nr, ni, w, stokesnr, Kepler_freq, vs)
               call vel_rd_centr(nr, w, stokesnr, vr, vn, realtime)

               do j_loop=1, nsws
                  !if (imax == j_loop) cycle
                  i = min(w, j_loop)
                  j = max(w, j_loop)
                  ij = (i-1)*nsws - (i-1)*(i-2)/2 + (j - i + 1)
                  call rel_vels(i, j, swarms, stokesnr, vr, vs, relvels(ij), vn, Reynolds, veta, &
                              & Vg2, tL, teta, Kepler_freq, cs_speed)

                  colrates_rp(i) = colrates_rp(i) - colrates(ij)
                  call col_rates(i, j, swarms, relvels(ij), colrates(ij), accelncol(ij), g%vol(nr,ni), deltar, deltaz)
                  colrates_rp(i) = colrates_rp(i) + colrates(ij)
               enddo
            endif


            ! find the right group to split: within the 90% from the total mass
#ifdef MULTI_COMPONENT
            mmax = maxval(swarms(:)%mswarm*(1.+swarms(:)%w), mask = swarms(:)%npar >= 2.0)
            idx  = pack([(i, i=1,size(swarms))],(swarms(:)%mswarm*(1.+swarms(:)%w) >= 0.9*mmax) .and. (swarms(:)%npar >= 2.))

#else
            mmax = maxval(swarms(:)%mswarm, mask = swarms(:)%npar >= 2.0)
            idx  = pack([(i, i=1,size(swarms))],(swarms(:)%mswarm >= 0.9*mmax) .and. (swarms(:)%npar >= 2.))
#endif
            ! choose randomly one 
            call random_number(r)
            imax = idx(1 + int(r*size(idx)))

            swarms(imax)%coll_f =  1
            idnr = swarms(nrk)%idnr
            swarms(imax)%npar = swarms(imax)%npar/2.
            swarms(imax)%mswarm = swarms(imax)%mswarm/2.
            ! copy i in k
!#ifdef MULTI_COMPONENT
            rpos = swarms(nrk)%rdis
            zpos = swarms(nrk)%zdis
!#endif
            swarms(nrk) = swarms(imax) 
            swarms(nrk)%idnr = idnr
!#ifdef MULTI_COMPONENT
            swarms(nrk)%rdis = rpos
            swarms(nrk)%zdis = zpos
!#endif
            ! then we need to calculate all properties but if this was the case we could do it once. To do
            call stokes_nr_centr(imax, swarms, stokesnr, lmfp, gasdens, Kepler_freq, cs_speed)
            call vel_vs_centr(nr, ni, imax, stokesnr, Kepler_freq, vs)
            call vel_rd_centr(nr, imax, stokesnr, vr, vn, realtime)

            do j_loop=1, nsws
               if (imax == j_loop) cycle
               i = min(imax, j_loop)
               j = max(imax, j_loop)
               ij = (i-1)*nsws - (i-1)*(i-2)/2 + (j - i + 1)
               call rel_vels(i, j, swarms, stokesnr, vr, vs, relvels(ij), vn, Reynolds, veta, &
                           & Vg2, tL, teta, Kepler_freq, cs_speed)

               colrates_rp(i) = colrates_rp(i) - colrates(ij)
               call col_rates(i, j, swarms, relvels(ij), colrates(ij), accelncol(ij), g%vol(nr,ni), deltar, deltaz)
               colrates_rp(i) = colrates_rp(i) + colrates(ij)
            enddo
         endif

         call stokes_nr_centr(nri, swarms, stokesnr, lmfp, gasdens, Kepler_freq, cs_speed)
         call stokes_nr_centr(nrk, swarms, stokesnr, lmfp, gasdens, Kepler_freq, cs_speed)
         call vel_vs_centr(nr, ni, nri, stokesnr, Kepler_freq, vs)
         call vel_vs_centr(nr, ni, nrk, stokesnr, Kepler_freq, vs)
         call vel_rd_centr(nr, nri, stokesnr, vr, vn, realtime)
         call vel_rd_centr(nr, nrk, stokesnr, vr, vn, realtime)
         do j_loop=1, nsws
            if (nri == j_loop) cycle
            i = min(nri, j_loop)
            j = max(nri, j_loop)
            ij = (i-1)*nsws - (i-1)*(i-2)/2 + (j - i + 1)
            call rel_vels(i, j, swarms, stokesnr, vr, vs, relvels(ij), vn, Reynolds, veta, &
                        & Vg2, tL, teta, Kepler_freq, cs_speed)

            colrates_rp(i) = colrates_rp(i) - colrates(ij)
            call col_rates(i, j, swarms, relvels(ij), colrates(ij), accelncol(ij), g%vol(nr,ni), deltar, deltaz)
            colrates_rp(i) = colrates_rp(i) + colrates(ij)
         enddo

         do j_loop=1, nsws
            if (nrk == j_loop) cycle
            i = min(nrk, j_loop)
            j = max(nrk, j_loop)
            ij =  (i-1)*nsws - (i-1)*(i-2)/2 + (j - i + 1)
            call rel_vels(i, j, swarms, stokesnr, vr, vs, relvels(ij), vn, Reynolds, veta, &
                        & Vg2, tL, teta, Kepler_freq, cs_speed)
            colrates_rp(i) = colrates_rp(i) - colrates(ij)
            call col_rates(i, j, swarms, relvels(ij), colrates(ij), accelncol(ij), g%vol(nr,ni), deltar, deltaz)
            colrates_rp(i) = colrates_rp(i) + colrates(ij)
         enddo

         kk = kk + 1
      enddo
      !-------------------------------------------------------------------------------------
      !write(*,*) '       collisions in zone',nr,ni,'done: total',kk, 'collisions'
      deallocate (colrates, accelncol, relvels, ij_back)
      deallocate (stokesnr, vs, vr, colrates_rp)
      
      do k = 1, nsws
         if(swarms(k)%coll_f /= 0) then
            l = first_idx
            do while (swrm(l)%idnr /= swarms(k)%idnr)
               l = l + 1
            enddo
            swrm(l) = swarms(k)
         else
            cycle
         endif   
      enddo
      !write(*,*) '       swrm updated!'
      nullify(swarms)

      return
   end subroutine mc_collisions

#endif





#ifdef KERNEL_TEST
   ! the routine performes collisional evolution on swarms located in the cell nr,ni with the aid of the MC algorithm
   subroutine mc_collisions_test(nr, ni, bin, swrm, dtime, realtime, kk, first_idx)
      implicit none
      type(swarm), dimension(:), allocatable, target            :: swrm
      type(list_of_swarms), dimension(:,:), allocatable, target :: bin
      type(swarm), dimension(:), allocatable              :: swarms      ! local rps array
      integer, intent(in)                             :: nr, ni      ! indices of the cell
      real, intent(inout)                                :: dtime       ! time step
      real, intent(in)                                :: realtime    ! physical time
      integer, intent(in)                             :: first_idx   ! first index of swrm in rbin
      integer                                         :: nsws        ! number of representative particles in given cell
      integer                                         :: nri, nrk    ! indices of physical and representative particles choosen to the next collision
      real, dimension(:), allocatable                 :: colrates, colrates_rp, accelncol    ! collision rates matrix
      real                                            :: local_time, dt ! physical time spent in the cell, time step between two collisions
      real                                            :: rand        ! random number
      real                                            :: totr        ! total collision rate
      integer                                         :: i, k, j, l, ij, Ncol, j_loop, idnr, imax
      integer, dimension(:,:), allocatable            :: ij_back
      integer, intent(out)                            :: kk ! collisions counter
      real :: mass_variation = 0. !0.01
      real, dimension(:), allocatable :: t_arr
      integer :: cont_time
      real :: fin
      type(hdf5_file_t)                                               :: file
      real :: r, mmax
      integer, allocatable :: idx(:)
      integer :: w, loops
      character(len=100)               :: loops_str
      character(len=100)               :: command
      real :: mswarm_init 
      mswarm_init = 10.e20
#ifdef TEST1
      dtime = 100000.
      allocate(t_arr(6))
      t_arr = (/1., 10., 100., 1000., 10000., 100000./)
#endif
#ifdef TEST2
      dtime = 20.
      allocate(t_arr(5))
      t_arr = (/4., 8., 12., 16., 20./)
#endif
#ifdef TEST3
      dtime = 0.9
      allocate(t_arr(3))
      t_arr = (/0.4, 0.7, 0.9/)
#endif
      allocate(swarms(size(swrm)))
      do loops=1, repeat
         swrm(:)%mass = 1.
         swrm(:)%mswarm = mswarm_init
         swrm(:)%npar = swrm(:)%mswarm/swrm(:)%mass


         !-------------------------------------------------------------------------------------
         
         swarms = swrm       ! points to the swarms array that contains only swarms located in the current cell nr,ni
         swarms(:)%coll_f = 0         ! setting collision flag to zero, will be updated to 1 if collisions happen
         nsws = size(swarms)          ! number of swarms in the current cell


         Ncol = nsws*(nsws+1)/2 !nsws*(nsws-1)/2       ! number of all potential collisions
         allocate (colrates(Ncol), accelncol(Ncol), ij_back(Ncol, 2), colrates_rp(nsws))


         do i=1, nsws
            do j=i, nsws ! i+1, nsws
               ij = (i-1)*nsws - (i-1)*(i-2)/2 + (j - i + 1)!(i-1)*nsws -(i-1)*i/2 + (j-i)
               ij_back(ij,1) = i
               ij_back(ij,2) = j 
            enddo
         enddo

         colrates_rp(:) = 0.
         ! calculates relative velocities and collision rates matrix
         do ij = 1, Ncol
            i = ij_back(ij,1)
            j = ij_back(ij,2)
            call col_rates(i, j, swarms, colrates(ij), accelncol(ij))
            colrates_rp(i) = colrates_rp(i) + colrates(ij)
         enddo

         write(*,*) 'ntot is ', nsws
         !------------ MAIN COLLISIONS LOOP ----------------------------------------------------
         local_time = 0.0
         kk = 0
         cont_time = 1

         write(*,*) 'initially total sum is ', sum(swarms(:)%mswarm), sum(swarms(:)%mass*swarms(:)%npar)
         do 
            if (t_arr(cont_time)<local_time) then
               write(loops_str, '(I0)') loops
               write(command, '(A,A,A)') './directory.sh ', trim(datadir), trim(loops_str)
               CALL SYSTEM(command)
               open(unit=2, file='outputs/path.txt', action='read')
               read(2,'(A)') path
               close(2)
               CALL SYSTEM('rm -rf outputs/path.txt')
               call hdf5_file_write(file, swarms, local_time, cont_time-1, local_time)
               write(*,*) local_time
               cont_time = cont_time + 1
               if (cont_time>size(t_arr)) exit
            endif

            totr = sum(colrates_rp)            ! total collision rate

            call random_number(rand)
            dt = -1. * log(rand) / totr       ! time step between this and the next collision

            local_time = local_time + dt      ! update of the local time

            ! ----- calculate which particles collide!
            ! select representative particle nri to undergo the collision
            call random_number(rand)
            rand = rand * totr
            j = 1
            fin = colrates_rp(1)
            do while (rand > fin)
               fin = fin + colrates_rp(j+1)
               j = j + 1
            enddo
            nri = j

            ! select physical particle nrk to undergo the collision
            call random_number(rand)
            rand = rand * colrates_rp(nri)

            nrk = nri ! we know that nri>=1
            i = min(nri, nrk)
            j = max(nri, nrk)
            ij = (i-1)*nsws - (i-1)*(i-2)/2 + (j - i + 1)
            fin = colrates(ij)
            do while ((rand > fin) .and. (nrk < nsws))
               nrk = nrk + 1
               i = min(nri, nrk)
               j = max(nri, nrk)
               ij = (i-1)*nsws - (i-1)*(i-2)/2 + (j - i + 1)
               fin = fin + colrates(ij)
            enddo

            ! the one that has less npar is nri
            if (swarms(i)%npar < swarms(j)%npar) then
               nri = i
               nrk = j
            else 
               nri = j
               nrk = i
            endif

            call collision(nri, nrk, swarms, accelncol)

            ! if one of them got emptied, refill it; nri cannot get emptied
            !
            if (swarms(nrk)%mswarm < mswarm_init*tolerance) then
               if (swarms(nrk)%npar < 1.) then
                  swarms(nrk)%mass = 0. ! to avoid it to split
                  if (swarms(nrk)%npar<0.)  then
                     write(*,*) swarms(nrk)%npar, swarms(nrk)%mass, swarms(nrk)%mswarm
                     stop
                  endif
               else
                  ! find most similar particle to nrk
                  w = minloc(abs(swarms(:)%mass-swarms(nrk)%mass), dim=1, mask=[(i/= nrk, i = 1, size(swarms))])
                  ! merge (not collide, compute average)
                  swarms(w)%mass = (swarms(w)%mass*swarms(w)%mswarm + swarms(nrk)%mass*swarms(nrk)%mswarm)/&
                                 & (swarms(w)%mswarm+swarms(nrk)%mswarm)
                  swarms(w)%npar = (swarms(w)%mswarm + swarms(nrk)%mswarm)/swarms(w)%mass
                  swarms(w)%mswarm = swarms(w)%mass*swarms(w)%npar
                  swarms(w)%coll_f =  1 ! we want this particle to update as well
                  swarms(nrk)%mass = 0.
                  swarms(nrk)%npar = 0.
                  swarms(nrk)%mswarm = 0.
                  do j_loop=1, nsws
                     !if (imax == j_loop) cycle
                     i = min(w, j_loop)
                     j = max(w, j_loop)
                     ij = (i-1)*nsws - (i-1)*(i-2)/2 + (j - i + 1)
                     colrates_rp(i) = colrates_rp(i) - colrates(ij)
                     call col_rates(i, j, swarms, colrates(ij), accelncol(ij))
                     colrates_rp(i) = colrates_rp(i) + colrates(ij)
                  enddo

               endif

               ! find the right group to split: within the 90% from the total mass
               mmax = maxval(swarms(:)%mswarm, mask = swarms(:)%npar >= 2.0)
               idx  = pack([(i, i=1,size(swarms))],(swarms(:)%mswarm >= 0.9*mmax) .and. (swarms(:)%npar >= 2.))
               ! choose randomly one 
               call random_number(r)
               imax = idx(1 + int(r*size(idx)))

               swarms(imax)%coll_f =  1 ! we want this particle to update as well
               idnr = swarms(nrk)%idnr

               swarms(imax)%npar = swarms(imax)%npar/2.
               swarms(imax)%mswarm = swarms(imax)%mswarm/2.
               ! copy i in k
               swarms(nrk) = swarms(imax)
               swarms(nrk)%idnr = idnr


               do j_loop=1, nsws
                  !if (imax == j_loop) cycle
                  i = min(imax, j_loop)
                  j = max(imax, j_loop)
                  ij = (i-1)*nsws - (i-1)*(i-2)/2 + (j - i + 1)
                  colrates_rp(i) = colrates_rp(i) - colrates(ij)
                  call col_rates(i, j, swarms, colrates(ij), accelncol(ij))
                  colrates_rp(i) = colrates_rp(i) + colrates(ij)
               enddo
            endif
      
            do j_loop=1, nsws
               !if (nri == j_loop) cycle
               i = min(nri, j_loop)
               j = max(nri, j_loop)
               ij = (i-1)*nsws - (i-1)*(i-2)/2 + (j - i + 1)
               colrates_rp(i) = colrates_rp(i) - colrates(ij)
               call col_rates(i, j, swarms, colrates(ij), accelncol(ij))
               colrates_rp(i) = colrates_rp(i) + colrates(ij)
            enddo

            do j_loop=1, nsws
               !if (nrk == j_loop) cycle
               i = min(nrk, j_loop)
               j = max(nrk, j_loop)
               ij = (i-1)*nsws - (i-1)*(i-2)/2 + (j - i + 1)
               colrates_rp(i) = colrates_rp(i) - colrates(ij)
               call col_rates(i, j, swarms, colrates(ij), accelncol(ij))
               colrates_rp(i) = colrates_rp(i) + colrates(ij)
            enddo

            kk = kk + 1
         enddo

         deallocate (colrates, ij_back, accelncol, colrates_rp)
      enddo
      
      write(*,*) 'Finally total sum is ', sum(swarms(:)%mswarm), sum(swarms(:)%mass*swarms(:)%npar)
      !-------------------------------------------------------------------------------------
      !write(*,*) '       collisions in zone',nr,ni,'done: total',kk,'collisions'

   

      !write(*,*) '      Updating swrm:'            ! TODO: update only the modified swarms
      do k = 1, nsws
         if(swarms(k)%coll_f /= 0) then
            l = FINDLOC(swrm(:)%idnr,swarms(k)%idnr,dim=1)
            swrm(l) = swarms(k)
         else
            cycle
         endif   
      enddo
      !write(*,*) '       swrm updated!'
      deallocate(swarms)
      deallocate(t_arr)

      return
   end subroutine mc_collisions_test

#endif


   ! calculating Stokes numbers of particle "i" IN THE CENTER OF CELL nr,ni
   subroutine stokes_nr_centr(i, swarms, stokesnr, lmfp, gasdens, Kepler_freq, cs_speed)
      implicit none
      type(swarm), dimension(:)                       :: swarms
      integer, intent(in)                             :: i
      real, dimension(:), allocatable                 :: stokesnr
      real, intent(in)                                :: lmfp, gasdens, Kepler_freq, cs_speed
      real                                            :: rad
#ifdef MULTI_COMPONENT
      rad = (0.75 / pi / swarms(i)%rhoi * swarms(i)%mass)**third                       ! particle radius
#else
      rad = con2 * swarms(i)%mass**third                       ! particle radius
#endif
      if (rad > 2.25 * lmfp) then ! Stokes regime
#ifdef MULTI_COMPONENT
         stokesnr(i) = sqrt(2.*pi) * swarms(i)%rhoi * AH2 * rad**2. * Kepler_freq / (9. * mH2 * cs_speed )
#else
         stokesnr(i) = sqrt(2.*pi) * matdens * AH2 * rad**2. * Kepler_freq / (9. * mH2 * cs_speed )
#endif
      else                        ! Epstein regime
#ifdef MULTI_COMPONENT
         stokesnr(i) = rad * swarms(i)%rhoi / (sqrt(8./pi) * cs_speed  * gasdens) * Kepler_freq
#else
         stokesnr(i) = rad * matdens / (sqrt(8./pi) * cs_speed  * gasdens) * Kepler_freq
#endif
      endif

      return
   end subroutine stokes_nr_centr

   ! calculation of relative velocities of bodies:
   ! we take 5 sources:
   ! Brownian motion vB
   ! turbulence vT (Ormel & Cuzzi 2007), implementation stolen from Til Birnstiel
   ! radial drift vr
   ! vertical settling vs
   ! azimuthal drift vtan
   subroutine rel_vels(ni, nj, swarms, stokesnr, vr, vs, relvels, vn, &
                     & Reynolds, veta, Vg2, tL, teta, Kepler_freq, cs_speed)
      implicit none
      type(swarm), dimension(:), pointer              :: swarms
      integer, intent(in)                             :: ni, nj
      real                                            :: relvels
      real, dimension(:), allocatable                 :: stokesnr, vr, vs
      real, intent(in)                                :: vn
      real                                            :: vB2, vT2, vtan
      real                                            :: gts, lts
      real, intent(in)                                :: tL, teta
      real, intent(in)                                :: Reynolds
      real, intent(in)                                :: veta, Vg2
      real                                            :: St1, St2
      real                                            :: y_star, c1, c2, c3, c0, eps, hulp1, hulp2
      real, parameter                                 :: ya = 1.6
      real, intent(in)                                :: Kepler_freq, cs_speed

      ! Brownian motions
      vB2 = 8.* mH2 * cs_speed**2 * (swarms(ni)%mass+swarms(nj)%mass) / (pi*swarms(ni)%mass*swarms(nj)%mass)

      ! turbulence
      if (stokesnr(nj) > stokesnr(ni)) then
         gts = stokesnr(nj) / Kepler_freq
         lts = stokesnr(ni) / Kepler_freq
         St1 = stokesnr(nj)
         St2 = stokesnr(ni)
      else
         gts = stokesnr(ni) / Kepler_freq
         lts = stokesnr(nj) / Kepler_freq
         St1 = stokesnr(ni)
         St2 = stokesnr(nj)
      endif

      if (gts < 0.2*teta) then
            vT2 = 1.5 *(veta/teta *(gts - lts))**2.0

      elseif (gts < teta/ya) then
            vT2 = Vg2 *(St1-St2)/(St1+St2)*(St1**2.0/(St1+Reynolds**(-0.5)) - St2**2.0/(St2+Reynolds**(-0.5)))
      elseif (gts < 5.0*teta) then
         !Eq. 17 of OC07. The second term with St_i**2.0 is negligible (assuming !Re>>1)
         !hulp1 = Eq. 17; hulp2 = Eq. 18

         hulp1 = ( (St1-St2)/(St1+St2) * (St1**2.0/(St1+ya*St1) - St2**2.0/(St2+ya*St1)) )!note the -sign
         hulp2 = 2.0*(ya*St1-Reynolds**(-0.5)) + St1**2.0/(ya*St1+St1) - St1**2.0/(St1+Reynolds**(-0.5)) +&
                                          St2**2.0/(ya*St1+St2) - St2**2.0/(St2+Reynolds**(-0.5))
         vT2 = Vg2 *(hulp1 + hulp2)

      elseif (gts < tL*0.2)  then
         eps=St2/St1!stopping time ratio
         vT2 = Vg2 *( St1*(2.0*ya - (1.0+eps) + 2.0/(1.0+eps) *(1.0/(1.0+ya) + eps**3.0/(ya+eps) )) )

      elseif (gts < tL) then
         !now y* lies between 1.6 (St1 << 1) and 1.0 (St1>=1). The fit below fits ystar to less than 1%
         c3 =-0.29847604
         c2 = 0.32938936
         c1 =-0.63119577
         c0 = 1.6015125
         y_star = c0 + c1*St1 + c2*St1**2.0 + c3*St1**3.0
         !we can then employ the same formula as before
         eps=St2/St1
         vT2 = Vg2 *( St1*(2.0*y_star - (1.0+eps) + 2.0/(1.0+eps) *(1.0/(1.0+y_star) + eps**3.0/(y_star+eps) )) )

      else
         vT2 = Vg2 *( 1.0/(1.0+St1) + 1.0/(1.0+St2) )
      endif


      ! tangential
      vtan = vn * ( 1. /(1.+stokesnr(nj)**2.) - 1. / (1.+stokesnr(ni)**2.) )

      ! total
      relvels = sqrt(vB2  + vT2 + (vs(nj) - vs(ni))**2 + (vr(nj) - vr(ni))**2 + vtan**2)


      return
   end subroutine rel_vels


#ifdef KERNEL_TEST
   real function kernel(m1, m2)
      implicit none
      real :: m1, m2

#ifdef TEST1
      kernel = 1.
#endif
#ifdef TEST2
      kernel = 0.5*(m1 + m2)

#endif

#ifdef TEST3
      kernel = m1*m2
#endif

      return
   end function

   subroutine col_rates(ni, nk, swarms, colrates, accelncol)
      implicit none
      integer, intent(in)                             :: ni, nk       ! number of line of colrates to update
      type(swarm), dimension(:), allocatable              :: swarms
      real                               :: vol      ! volume of the cell
      real :: N_frac, colrates, accelncol, dmmax_frac
      integer :: i, j 

      if (swarms(ni)%npar<swarms(nk)%npar) then
         i = ni
         j = nk
      else
         i = nk
         j = ni
      endif
      
      N_frac = swarms(j)%npar

      if (ni == nk) then
         if (N_frac>=2.) then ! there needs to be at least two bodies to collide with themselves
            N_frac = 0.5*N_frac
         else
            N_frac = 0.
         endif
      endif


      vol = real(size(swarms)) * 10.e20 
      colrates = N_frac*kernel(swarms(ni)%mass, swarms(nk)%mass)/vol

      ! group collisions because otherwise col rate too high when Ni<< Nj
 
      dmmax_frac = max(swarms(i)%npar/swarms(j)%npar, swarms(j)%mass/swarms(i)%mass) ! do not let to group if there are not enough particles

      ! adaptive dmmax
      if (dmmax_frac < dmmax) then
         accelncol = dmmax/dmmax_frac
      else
         accelncol = 1.
      endif

      colrates = colrates/accelncol
      ! so the collision rate is supressed
      if (swarms(nk)%npar < 1.0) colrates = 0.0
      if (swarms(ni)%npar < 1.0) colrates = 0.0
      return
   end subroutine col_rates


#else
   ! calculation of the collision rates between representative particle nl and all physical particles
   subroutine col_rates(ni, nk, swarms, relvels, colrates, accelncol, vol, deltar, deltaz)
      implicit none
      integer, intent(in)                             :: ni, nk       ! number of line of colrates to update
      type(swarm), dimension(:), pointer              :: swarms
      real                                            :: colrates
      real                                            :: accelncol
      real                                            :: relvels
      real, intent(in)                                :: vol      ! volume of the cell
      real, intent(in)                                :: deltar, deltaz ! size of radial and vertical grid
      real :: N_frac, dmmax_frac
      integer :: i,j

      if (swarms(ni)%npar<swarms(nk)%npar) then
         i = ni
         j = nk
      else
         i = nk
         j = ni
      endif
      
      N_frac = swarms(j)%npar

      if (ni == nk) then
         if (N_frac>=2.) then ! there needs to be at least two bodies to collide with themselves
            N_frac = 0.5*N_frac
         else
            N_frac = 0.
         endif
      endif

      
#ifdef MULTI_COMPONENT
      colrates = N_frac * relvels * pi**third * (0.75)**(2. * third) * &
                  ((swarms(ni)%mass/swarms(ni)%rhoi)**third + (swarms(nk)%mass/swarms(nk)%rhoi)**third)**2./vol
#else
      colrates = N_frac * relvels * con1 * (swarms(ni)%mass**third + swarms(nk)%mass**third)**2./vol
#endif
      dmmax_frac = max(swarms(i)%npar/swarms(j)%npar, swarms(j)%mass/swarms(i)%mass) ! do not let to group if there are not enough particles

      ! adaptive dmmax
      if (dmmax_frac< dmmax) then
         accelncol = dmmax/dmmax_frac
      else
         accelncol = 1.
      endif
      
      colrates = colrates / accelncol


      if (swarms(nk)%npar < 1.0) colrates = 0.0
      if (swarms(ni)%npar < 1.0) colrates = 0.0 

      return
   end subroutine col_rates
#endif


   ! choosing particles to the next collision
   ! nri -> representative
   ! nrk -> physical
   subroutine choose_swarms(nri, nrk, swarms, colrates_rp, colrates, totrate)
      implicit none
      integer, intent(out)                            :: nri, nrk
      type(swarm), dimension(:), pointer              :: swarms
      real, dimension(:), allocatable, intent(in)     :: colrates_rp
      real, dimension(:), allocatable, intent(in)     :: colrates
      integer                                         :: ij, j, i
      real, intent(in)                                :: totrate
      real, dimension(2)                              :: rand
      real                                            :: fin
      real :: nsws

      nsws = size(swarms)

      ! ----- calculate which particles collide!
      ! select representative particle nri to undergo the collision
      call random_number(rand)
      rand(1) = rand(1) * totrate
      j = 1

      fin = colrates_rp(1)
      do while (rand(1) > fin)
         fin = fin + colrates_rp(j+1)
         j = j + 1
      enddo
      nri = j

      ! select physical particle nrk to undergo the collision
      rand(2) = rand(2) * colrates_rp(nri)

      nrk = nri ! we know that nri>=1
      i = min(nri, nrk)
      j = max(nri, nrk)
      ij = (i-1)*nsws - (i-1)*(i-2)/2 + (j - i + 1)
      fin = colrates(ij)
      do while ((rand(2) > fin) .and. (nrk < nsws))
         nrk = nrk + 1
         i = min(nri, nrk)
         j = max(nri, nrk)
         ij = (i-1)*nsws - (i-1)*(i-2)/2 + (j - i + 1)
         fin = fin + colrates(ij)
      enddo

      ! the one that has less npar is nri
      if (swarms(i)%npar < swarms(j)%npar) then
         nri = i
         nrk = j
      else 
         nri = j
         nrk = i
      endif

      return
   end subroutine choose_swarms

#ifdef KERNEL_TEST
   ! performing the collision: deciding the collision outcome - put your collision model here
   ! only the representative particle is updated
subroutine collision(nri,nrk,swarms, accelncol)
   implicit none
   integer, intent(in)                             :: nri, nrk
   type(swarm), dimension(:), allocatable              :: swarms
   real, dimension(:), allocatable                 :: accelncol
   real                                            :: rvel
   integer :: ij, i, j
   i = min(nri, nrk)
   j = max(nri, nrk)
   ij = (i-1)*size(swarms) - (i-1)*(i-2)/2 + (j - i + 1)!(i-1)*size(swarms) -(i-1)*i/2 + (j-i)


   call hit_and_stick(nri,nrk, ij, swarms, accelncol)

   swarms(nri)%coll_f =  1
   swarms(nrk)%coll_f =  1
   return
end subroutine collision


#else
   ! performing the collision: deciding the collision outcome - put your collision model here
   ! only the representative particle is updated
   subroutine collision(nri,nrk,swarms,relvels, accelncol)
      implicit none
      integer, intent(in)                             :: nri, nrk
      type(swarm), dimension(:), pointer              :: swarms
      real, dimension(:), allocatable                 :: accelncol
      real, dimension(:), allocatable                 :: relvels
      real                                            :: rvel
      integer :: ij, i, j
      real :: ran

      i = min(nri, nrk)
      j = max(nri, nrk)
      ij = (i-1)*size(swarms) - (i-1)*(i-2)/2 + (j - i + 1)!(i-1)*size(swarms) -(i-1)*i/2 + (j-i)
      rvel = relvels(ij)

      if (rvel < vfrag) then
         call hit_and_stick(nri,nrk, ij, swarms, accelncol)
#ifdef EROSION
      ! this need to be fixed
      else if (swarms(nri)%mass/swarms(nrk)%mass .ge. erosion_mass_ratio) then
         call erosion(nri,nrk,swarms,accelncol)
#endif
      else ! for comparison with mcdust, I should use the same frag model
      
         ! assume they both fragment
         call fragmentation(nri, swarms)
         call random_number(ran)
         if (swarms(nri)%npar/swarms(nrk)%npar>ran) then
            call fragmentation(nrk, swarms)
#ifdef MULTI_COMPONENT
            swarms(nrk)%npar = swarms(nrk)%mswarm/swarms(nrk)%mass*(1.+swarms(nrk)%w)
         endif
         swarms(nri)%npar = swarms(nri)%mswarm/swarms(nri)%mass*(1.+swarms(nri)%w)
#else
            swarms(nrk)%npar = swarms(nrk)%mswarm/swarms(nrk)%mass
         endif
         swarms(nri)%npar = swarms(nri)%mswarm/swarms(nri)%mass
#endif
      endif
      swarms(nri)%coll_f =  1
      swarms(nrk)%coll_f =  1
      return
      if((swarms(nri)%mass < 0.) .or.(swarms(nrk)%mass < 0.)) then
         write(*,*) 'negative mass', swarms(nri)%mass, 'in id', &
            swarms(nri)%idnr, swarms(nri)%zdis/AU, swarms(nri)%rdis/AU
         write(*,*) 'negative mass', swarms(nrk)%mass, 'in id', &
            swarms(nrk)%idnr, swarms(nrk)%zdis/AU, swarms(nrk)%rdis/AU
         stop
      endif
   end subroutine collision
#endif
   ! sticking collision
   subroutine hit_and_stick(nri,nrk, ij, swarms, accelncol)
      implicit none
      integer, intent(in)                             :: nri, nrk, ij
#ifdef KERNEL_TEST
      type(swarm), dimension(:), allocatable              :: swarms
#else
      type(swarm), dimension(:), pointer              :: swarms
#endif
      real, dimension(:), allocatable               :: accelncol
      ! we are going to save in particle where N1 is smaller the sticked part
      ! in the other the remanent

      ! nri is the one that is updated as npari < = npark
      ! the one with smallest mass is multiplied by accelnco

      if (nri == nrk) then
         swarms(nri)%mass = 2.*swarms(nri)%mass
         swarms(nri)%npar = 0.5*swarms(nri)%npar
      else


#ifdef MULTI_COMPONENT
         swarms(nri)%fw = (swarms(nri)%mass*swarms(nri)%fw + accelncol(ij) * swarms(nrk)%mass * swarms(nrk)%fw) / &
                              & (swarms(nri)%mass + accelncol(ij) * swarms(nrk)%mass)
         swarms(nri)%w = swarms(nri)%fw/(1.-swarms(nri)%fw)
         swarms(nri)%rhoi = 1. / ( (swarms(nri)%fw/matdensw) + ((1.-swarms(nri)%fw)/matdenssi))
         swarms(nri)%mass = swarms(nri)%mass + accelncol(ij) * swarms(nrk)%mass
         !npar for i1 stays the same
         swarms(nri)%mswarm = swarms(nri)%mass*swarms(nri)%npar/(1.+swarms(nri)%w) ! this is msi
   
         ! ----- update particle 2 -----
         ! mass stays the same
         swarms(nrk)%npar = swarms(nrk)%npar - swarms(nri)%npar*accelncol(ij)
         swarms(nrk)%mswarm = swarms(nrk)%mass*swarms(nrk)%npar/(1.+swarms(nrk)%w)
#else
         swarms(nri)%mass = swarms(nri)%mass + accelncol(ij) * swarms(nrk)%mass
         !npar for i1 stays the same
         swarms(nri)%mswarm = swarms(nri)%mass*swarms(nri)%npar
   
         ! ----- update particle 2 -----
         ! mass stays the same
         swarms(nrk)%npar = swarms(nrk)%npar - swarms(nri)%npar*accelncol(ij)
         swarms(nrk)%mswarm = swarms(nrk)%mass*swarms(nrk)%npar
#endif
      endif

      return
   end subroutine hit_and_stick

   ! fragmentation collision
   ! put your fragment size distribution here:
   ! n(m) ~ m^(kappa - 2)
   subroutine fragmentation(nri,swarms)
      implicit none
      integer, intent(in)                             :: nri
      real                                            :: ran
      type(swarm), dimension(:), pointer              :: swarms
      real, parameter                                 :: kappa = 1./6.  ! n(m) ~ m^(kappa - 2)

      call random_number(ran)
      swarms(nri)%mass = (ran * (swarms(nri)%mass**kappa - m0**kappa) +  m0**kappa )**(1./kappa)
      swarms(nri)%mass = max(swarms(nri)%mass, m0)

      return
   end subroutine fragmentation

#ifdef EROSION
   !erosion collision as in dustpy
   !when small aggregate hits large one, it fragments
   !and big aggregate loses a chunk mass same as target particle

   subroutine erosion(nri, nrk, swarms, accelncol)
      implicit none
      integer, intent(in)                             :: nri, nrk
      type(swarm), dimension(:), pointer              :: swarms
      real, dimension(:), allocatable                 :: accelncol
      real                                            :: ran, p
      

      p = swarms(nrk)%mass/swarms(nri)%mass
      call random_number(ran)
      if (ran .ge. p) then
            swarms(nri)%mass = swarms(nri)%mass - accelncol(nri, nrk) * swarms(nrk)%mass
      else
         swarms(nri)%mass = swarms(nrk)%mass
      endif
      return
   end subroutine erosion
#endif

   ! vertical settling velocity
   subroutine vel_vs_centr(nr, ni, i, stokesnr, Kepler_freq, vs)
      implicit none
      integer, intent(in)                             :: nr, ni
      real, dimension(:), allocatable, intent(in)     :: stokesnr
      real, dimension(:), allocatable                 :: vs
      integer, intent(in)                             :: i       ! index of particle
      real, intent (in)                               :: Kepler_freq
      vs(i) = g%zce(nr,ni) * Kepler_freq * min(stokesnr(i), 0.5) !stokesnr(i) / (1. + stokesnr(i)**2.) 
      !vs(i) = min(g%zce(nr,ni) * Kepler_freq * stokesnr(i), g%zce(nr,ni) * Kepler_freq) !stokesnr(i) / (1. + stokesnr(i)**2.) 

      return
   end subroutine vel_vs_centr


   ! velocity of radial drift
   subroutine vel_rd_centr(nr, i, stokesnr, vr, vn, realtime)
      implicit none
      integer, intent(in)                             :: nr
      real, dimension(:), allocatable                 :: stokesnr
      real, dimension(:), allocatable                 :: vr
      real, intent(in)                                :: vn
      integer, intent(in)                             :: i      ! index of particle
      real, intent(in)                                :: realtime

      vr(i) = (2. * vn * stokesnr(i)  + &
               vgas(g%rce(nr), realtime)) / (1. + stokesnr(i) * stokesnr(i))

      return
   end subroutine vel_rd_centr


end

