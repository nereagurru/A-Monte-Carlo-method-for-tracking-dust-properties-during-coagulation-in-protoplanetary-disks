! this module should contain all the initial conditions for dust
module initproblem

   use constants,  only: pi, third, AU
   use discstruct, only: gasmass, cs, omegaK
   use parameters, only: dtg, minrad0, maxrad0, r0, alpha_t, vfrag
#ifdef MULTI_COMPONENT
   use parameters, only: matdenssi, matdensw
   use discstruct, only: sigmag
#else
   use parameters, only: matdens
#endif
   use types,      only: swarm
   use advection,  only: stokesnr
   
   implicit none
   
   private
   public   :: init_swarms, init_random_seed, m0, nmonom0
#ifdef MULTI_COMPONENT
   public :: f_dens
#endif
   real     :: m0, nmonom0
   
   contains
   
   ! initializing the swrm array and some variables
   subroutine init_swarms(Ntot, swrm &
#ifdef MULTI_COMPONENT
                        &, vol_0D &
#endif
                        &)
      implicit none
      type(swarm), dimension(:), allocatable, target   :: swrm        ! list of all swarms in the simulation
      integer, intent(in)                              :: Ntot        ! number of all swarms
#ifdef MULTI_COMPONENT
      real :: vol_0D
#endif
      real                                             :: mdust       ! mass of the dust in the simulatated domain
      real, dimension(2)                               :: rand        ! random numbers
      real                                             :: Hg          ! pressure height scale of the gas
      integer                                          :: i
      real, parameter                                  :: s = 0.0     ! to initialize the particles density r slope (for MMSN s = -0.25)
      real                                             :: mswarm      ! initially all particles have the same mswarm + size
#ifdef MULTI_COMPONENT
      real                                             :: dz_0D, rmean, St0, a0
      dz_0D = 0.00001*AU
      rmean = 0.5*(minrad0+maxrad0)*AU
      vol_0D = 2.*pi*rmean*(maxrad0-minrad0)*AU*dz_0D
      a0 = 1. ! cm-sized initially
      St0 = pi/2.*a0*matdenssi/sigmag(rmean, 0.)  ! Epstein regime
      write(*,*) 'St0 and a0 are ', St0, a0
      write(*,*) 'midplane dust-to-gas ratio is ', dtg/(sqrt(alpha_t/(alpha_t+St0)))
      mdust = dtg/(sqrt(alpha_t/(alpha_t+St0))) * gasmass(minrad0*AU,maxrad0*AU,dz_0D, 0.0)
#else
      ! total mass of dust =  dust to gas ratio x mass of the gas               
      mdust = dtg * gasmass(minrad0*AU,maxrad0*AU,0.0)
#endif

      ! mass of one swarm
      mswarm = mdust / real(Ntot) ! this is for total
      ! monomer mass
#ifdef MULTI_COMPONENT
      m0 = 4. * third * pi * r0**3 * matdenssi
#else
      m0 = 4. * third * pi * r0**3 * matdens
#endif
      ! initial number of monomers
      nmonom0 = mdust/m0
      
      if (.not.allocated(swrm)) allocate( swrm(Ntot) )

      ! initializing the particles
      do i = 1, Ntot
         swrm(i)%idnr = i
#ifdef MULTI_COMPONENT
         swrm(i)%rdis = (minrad0 + (maxrad0-minrad0)*real(i-1)/real(Ntot-1))*AU
         swrm(i)%fw = 0.
         swrm(i)%w = swrm(i)%fw/(1.-swrm(i)%fw)
         swrm(i)%rhoi = f_dens(swrm(i)%w)
         swrm(i)%mass = 4. * third * pi * a0**3 * swrm(i)%rhoi
         swrm(i)%mswarm = mswarm ! we track mswarm_silicates; this is valid because initially all water vapor
         swrm(i)%npar = swrm(i)%mswarm * (1.+swrm(i)%w)/ swrm(i)%mass 
         swrm(i)%zdis = (-0.5*dz_0D + dz_0D*real(i-1)/real(Ntot-1))
#else
         call random_number(rand)
         swrm(i)%rdis = (((maxrad0*AU)**(s+1.) - (minrad0*AU)**(s+1.))*rand(1) + (minrad0*AU)**(s+1.))**(1./(s+1.))
#ifdef KERNEL_TEST
         swrm(i)%mass = 1.
         swrm(i)%mswarm = 10.e20
         swrm(i)%npar = mswarm / swrm(i)%mass
         call random_number(rand)
         swrm(i)%zdis = 0.001*AU* sqrt(-2.*log(rand(1))) * cos(2.*pi*rand(2)) 
#else
         ! global disk
         swrm(i)%mass = m0
         swrm(i)%mswarm = mswarm
         swrm(i)%npar = swrm(i)%mswarm/ swrm(i)%mass 

#ifdef VERTICALMODEL
         swrm(i)%rdis = minrad0 * AU + (maxrad0 - minrad0) * AU * (real(i)-0.5) / real(Ntot)   ! BE CAREFUL WITH 1D vertical column tests!!!!!!; you don't want the s then
#else
         call random_number(rand)
         swrm(i)%rdis = (((maxrad0*AU)**(s+1.) - (minrad0*AU)**(s+1.))*rand(1) + (minrad0*AU)**(s+1.))**(1./(s+1.))
#endif
         Hg = cs(swrm(i)%rdis) / omegaK(swrm(i)%rdis)
         call random_number(rand)
         swrm(i)%zdis = Hg * sqrt(-2.*log(rand(1))) * cos(2.*pi*rand(2))
#endif
#endif
         swrm(i)%stnr = stokesnr(swrm(i), 0.0)
         swrm(i)%velr = 0.0
         swrm(i)%velz = 0.0
         swrm(i)%coll_f = 0
      enddo
          
      return
   end subroutine init_swarms
   
   ! initialize the random number generator
   subroutine init_random_seed
      implicit none
      integer                            :: i, n, clock
      integer, dimension(:), allocatable :: seed

      call random_seed(size = n)
      allocate(seed(n))

      call system_clock(count=clock)

      seed = 37 * [(i - 1, i = 1, n)]
      seed = seed + clock

      call random_seed(put = seed)
      deallocate(seed)

   end subroutine init_random_seed

#ifdef MULTI_COMPONENT
   real function f_dens(w)
      implicit none
      real, intent(in) :: w

      f_dens = (1. + w)/(w/matdensw + 1./matdenssi)

      return
   end function 
#endif
end
