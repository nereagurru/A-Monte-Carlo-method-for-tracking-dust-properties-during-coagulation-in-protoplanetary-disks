! this module should contain all the initial conditions for dust
module initproblem

   use constants,  only: pi, third, AU
   use discstruct, only: gasmass, cs, omegaK
   use parameters, only: dtg, minrad0, maxrad0, r0
#ifdef MULTI_COMPONENT
   use parameters, only: matdenssi, matdensw
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
   subroutine init_swarms(Ntot, swrm)
      implicit none
      type(swarm), dimension(:), allocatable, target   :: swrm        ! list of all swarms in the simulation
      integer, intent(in)                              :: Ntot        ! number of all swarms
      real                                             :: mdust       ! mass of the dust in the simulatated domain
      real, dimension(2)                               :: rand        ! random numbers
      real                                             :: Hg          ! pressure height scale of the gas
      integer                                          :: i
      real, parameter                                  :: s = 0.0     ! to initialize the particles density r slope (for MMSN s = -0.25)
      real                                             :: mswarm      ! initially all particles have the same mswarm + size
      ! total mass of dust =  dust to gas ratio x mass of the gas               
      mdust = dtg * gasmass(minrad0*AU,maxrad0*AU,0.0)
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
         call random_number(rand)
         swrm(i)%rdis = (((maxrad0*AU)**(s+1.) - (minrad0*AU)**(s+1.))*rand(1) + (minrad0*AU)**(s+1.))**(1./(s+1.))
#ifdef MULTI_COMPONENT
         swrm(i)%fw = 0.
         swrm(i)%w = 0.
         swrm(i)%rhoi = f_dens(swrm(i)%w)
         swrm(i)%mass = 4. * third * pi * r0**3 * swrm(i)%rhoi
         swrm(i)%mswarm = mswarm/(1.+swrm(i)%w) ! we track mswarm_silicates
         swrm(i)%npar = swrm(i)%mswarm * (1.+swrm(i)%w)/ swrm(i)%mass 

         ! 0D simulation; but do not set all particles at same r or z. Otherwise grid volume 0.
         call random_number(rand)
         swrm(i)%zdis = 0.001 *AU* sqrt(-2.*log(rand(1))) * cos(2.*pi*rand(2))
#else

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
   real function f_dens(fw)
      implicit none
      real, intent(in) :: fw

      f_dens = (1. + fw)/(fw/matdensw + 1./matdenssi)

      return
   end function 
#endif
end
