!#######################################################################
!Important constants
!#######################################################################
module constants 
 
   use precision_parameters, only: dp
   implicit none 
   
   !Conversion factors
   real(dp),parameter :: atm_pa = 101325._dp                            !atm to Pa
   real(dp),parameter :: atm_bar = 1.e5_dp                              !atm to bar
   
   !Mathematical constants
   real(dp),parameter :: pi = dacos(-1._dp), twopi=2._dp*pi,&
      fourpi=4._dp*pi, pio2 = pi/2._dp, invpi = 1._dp/pi, &
      thrpio2 = 1.5_dp*pi                                               !Pi and multiples
   real(dp),parameter :: sqtwo = dsqrt(2._dp), sqthree = dsqrt(3._dp)   !Square roots
   real(dp), parameter :: &
      euler = 0.5772156649015328606065120900824024310422_dp             !Euler constant
   
   !Physical constants
   real(dp),parameter :: h_planck = 6.62607015e-34_dp                   !Planck constant [J s]
   real(dp),parameter :: c0_light = 299792458._dp                       !Speed of light [m/s]
   real(dp),parameter :: k_boltzmann = 1.380649e-23_dp                  !Boltzmann constant [J/K]
   real(dp),parameter :: c_planck1 = h_planck*c0_light**2._dp           !First Planck function constant as 
                                                                        !  defined by Siegel and Howell [W m2/sr]
   real(dp),parameter :: c_planck1_mod = twopi*c_planck1                !First Planck function constant as
                                                                        !  defined by Modest [W m2/sr]
   real(dp),parameter :: c_planck2 = h_planck*c0_light/k_boltzmann      !Second Planck function constant [m K]
   
   real(dp),parameter :: sigma = (c_planck1_mod*pi**4._dp)/&            !Stefan-Boltzmann constant [W/(m2 K4)]
                                 (15._dp*c_planck2**4._dp)
   real(dp),parameter :: avogadro = 6.02214076e23_dp                    !Avogadro number
   real(dp),parameter :: Ru = avogadro*k_boltzmann                      !Gas constant [J/(K mol)]
   real(dp),parameter :: c_los = avogadro*atm_pa*1.e-6_dp/Ru            !Loschmidt constant
   
   !Aotmic weights (reference: CRC Handbook of Chemistry and Physics)
   real(dp),parameter :: MW_C = 12.011_dp, MW_H = 1.008_dp, &
      MW_N = 14.007_dp, MW_O = 15.999_dp

   !Molecular weights of common molecules
   real(dp),parameter :: MW_CH4 = MW_C + 4._dp*MW_H
   real(dp),parameter :: MW_CO = MW_C + MW_O, MW_CO2 = MW_C + 2._dp*MW_O
   real(dp),parameter :: MW_H2 = 2._dp*MW_H, MW_H2O = MW_H2 + MW_O
   real(dp),parameter :: MW_O2 = 2._dp*MW_O, MW_OH = MW_O + MW_H
   real(dp),parameter :: MW_N2 = 2._dp*MW_N, MW_NO = MW_N + MW_O
   
   !Misc constants
   real(dp),parameter :: sigrpi = sigma*invpi
   
contains

   !==============================================================
   !Subroutine to print all constants (for verification purposes)
   !==============================================================
   subroutine print_constants

      write(*,*)
      write(*,'(a)') 'CONVERSION FACTORS'
      write(*,'(a,es23.16)') 'atm to Pa = ',atm_pa
      write(*,'(a,es23.16)') 'atm to bar = ',atm_bar
      write(*,*)
      
      write(*,*)
      write(*,'(a)') 'MATHEMATICAL CONSTANTS'
      write(*,'(a,es23.16)') 'pi = ',pi
      write(*,'(a,es23.16)') '2pi = ',twopi
      write(*,'(a,es23.16)') '4pi = ',fourpi
      write(*,'(a,es23.16)') 'pi/2 = ',pio2
      write(*,'(a,es23.16)') '3pi/2 = ',thrpio2
      write(*,'(a,es23.16)') '1/pi = ',invpi
      write(*,'(a,es23.16)') 'sqrt(2) = ',sqtwo
      write(*,'(a,es23.16)') 'sqrt(3) = ',sqthree
      write(*,'(a,es23.16)') 'Euler constant = ',euler
      write(*,*)
   
      write(*,*)
      write(*,'(a)') 'PHYSICAL CONSTANTS'
      write(*,'(a,es23.16)') 'Planck constant [J s] = ',h_planck
      write(*,'(a,es23.16)') 'Speed of light = c0_light',c0_light
      write(*,'(a,es23.16)') 'Boltzmann constant = ',k_boltzmann
      write(*,'(a,es23.16)') &
         'First Planck constant (Siegel) [W m2/sr] = ',c_planck1
      write(*,'(a,es23.16)') &
         'First Planck constant (Modest) [W m2/sr] = ',c_planck1_mod
      write(*,'(a,es23.16)') &
         'Second Planck constant [m K] = ',c_planck2
      write(*,'(a,es23.16)') &
         'Stefan-Boltzmann constant [W/(m2 K4)] = ',sigma
      write(*,'(a,es23.16)') 'Avogadro number = ',avogadro
      write(*,'(a,es23.16)') 'Gas constant [J/(K mol)] = ',Ru
      write(*,'(a,es23.16)') 'Loschmidt constant = ',c_los

      write(*,*)
      write(*,'(a)') 'ATOMIC WEIGHTS'
      write(*,'(a,es23.16)') 'Carbon:',MW_C
      write(*,'(a,es23.16)') 'Hydrogen:',MW_H
      write(*,'(a,es23.16)') 'Nitrogen:',MW_N
      write(*,'(a,es23.16)') 'Oxygen:',MW_O

      write(*,*)
      write(*,'(a)') 'MOLECULAR WEIGHTS'
      write(*,'(a,es23.16)') 'CH4:',MW_CH4
      write(*,'(a,es23.16)') 'CO:',MW_CO
      write(*,'(a,es23.16)') 'CO2:',MW_CO2
      write(*,'(a,es23.16)') 'H2:',MW_H2
      write(*,'(a,es23.16)') 'H2O:',MW_H2O
      write(*,'(a,es23.16)') 'O2:',MW_O2
      write(*,'(a,es23.16)') 'OH:',MW_OH
      write(*,'(a,es23.16)') 'N2:',MW_N2
      write(*,'(a,es23.16)') 'NO:',MW_NO

   endsubroutine print_constants
   
end module constants 
