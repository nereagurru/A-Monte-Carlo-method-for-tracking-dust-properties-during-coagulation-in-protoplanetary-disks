
#ifdef MULTI_COMPONENT
module phase_change

   use constants,    only: mH2O, kB, pi, year, third, AU
   use types,        only: swarm
   use initproblem,  only: m0, f_dens
   use discstruct,  only: Temp
   
   implicit none

   private
   public   :: sublimation, condensation

   contains

      ! sublimation event; remove all water
   subroutine sublimation(swrm, dtime, vapour_mass, rhovap, temper, vol)  
      implicit none
      type(swarm), dimension(:), allocatable, target   :: swrm        ! list of all swarms in the simulation
      real, intent(in)                                 :: dtime, rhovap, temper, vol
      real, intent(inout) :: vapour_mass
      integer :: i
      real :: vth, mwater_remove, mwater_remove_fw, total_mass_loss

      vth = sqrt(8*kB*temper/pi/mH2O)
      total_mass_loss = 0.

      do i=1, size(swrm)
         mwater_remove = dtime*pi*(swrm(i)%mass/swrm(i)%rhoi*3./4./pi)**(2.*third)*vth*(-vapour_mass/vol+rhovap)
         mwater_remove_fw = min(swrm(i)%mass*swrm(i)%w/(swrm(i)%w+1.) , mwater_remove) ! we cannot remove more water than there is


         swrm(i)%w = max(0., swrm(i)%w - (1.+swrm(i)%w)*mwater_remove_fw/swrm(i)%mass)
         swrm(i)%fw = swrm(i)%w/(1.+swrm(i)%w)
         total_mass_loss = total_mass_loss + mwater_remove_fw*swrm(i)%npar
         swrm(i)%mass = swrm(i)%mass - mwater_remove_fw 

         swrm(i)%npar = swrm(i)%mswarm / swrm(i)%mass *(1. + swrm(i)%w)
         swrm(i)%rhoi = f_dens(swrm(i)%w)
         
      enddo
      vapour_mass = vapour_mass + total_mass_loss


      return
   end subroutine



      ! condensation event; remove all water
   subroutine condensation(swrm, dtime, vapour_mass, rhovap, temper, vol)  
      implicit none
      type(swarm), dimension(:), allocatable, target   :: swrm        ! list of all swarms in the simulation
      real, intent(in)                                 :: dtime, rhovap, temper, vol
      real, intent(inout) :: vapour_mass
      integer  :: i
      real, dimension(:), allocatable :: grain_size
      real :: A_K, total_mass_loss, vth, total_grain_power ! follow Krijt et al. 2016 
      real :: gain_water

      vth = sqrt(8*kB*temper/pi/mH2O)
      allocate(grain_size(size(swrm)))

      grain_size(:) = (3./4./pi*swrm(:)%mass/swrm(:)%rhoi)**third
      total_grain_power = sum(swrm(:)%npar*grain_size(:)**2.)

      A_K = pi*vth/vol*total_grain_power! sum(swrm(:)%npar*grain_size(:)**2.)


      total_mass_loss = (1.-exp(-A_K*dtime))*(vapour_mass/vol-rhovap)*vol
      vapour_mass = vapour_mass - total_mass_loss

      do i=1, size(swrm)

         gain_water = total_mass_loss*grain_size(i)**2./total_grain_power ! per particle
         swrm(i)%w = swrm(i)%w + (1. + swrm(i)%w)*gain_water/swrm(i)%mass
         swrm(i)%fw = swrm(i)%w/(1.+swrm(i)%w)
         swrm(i)%mass = swrm(i)%mass + gain_water
         swrm(i)%npar = swrm(i)%mswarm / swrm(i)%mass *(1.+swrm(i)%w)
         swrm(i)%rhoi = f_dens(swrm(i)%w) 
      enddo
      !write(*,*) total_mass_loss

      deallocate(grain_size)

      return
   end subroutine

end module phase_change
#endif