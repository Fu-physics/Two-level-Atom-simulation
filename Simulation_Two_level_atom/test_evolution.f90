!!  gfortran CH_4_state_vector.f95 rk4.f90 -o CH_4 -lfftw3



!  --------------------  |a>, w_a = 10.0
!
!      |
!      |
!      |
!      |   \nu_pu = 1.0,   om_pu = ?
!      |
!      |
!      |
!
!  --------------------- |b>,  w_b = 0.0



module para
 implicit none

 complex*8 :: imag = sqrt((-1.0,0.0))
 real*8 :: pi = 4.0 * atan(1.0)

 integer :: i , j, p, q

 !!!!!           Unit of frequency is 10^14 HZ      !!!!!!!
 ! the atomic energy level frequencies
 real(8) :: w_a = 1.0D+00, &
            w_b = 0.0D+00

 ! the laser frequeny
 real(8) ::  nu_pu = 1.0D+00     
 real(kind = 8) :: om_pu = 0.05
            
 !!!!!!            Unit of time is    10^(-14) second !!!!!!
 ! the pulse parameters
 real(8) :: tau_pu = 100.0D+00

 real(8) :: t0_pu = 1000.0D+00

 contains
    FUNCTION Omega_pu(t)
    real(8) :: Omega_pu,  t
    !Omega_pu = 0.0
    !Omega_pu = sin(t)
    !Omega_pu = om_pu * cos(nu_pu*t)
    Omega_pu = om_pu * exp(-(t - t0_pu)*(t - t0_pu)/(2*tau_pu*tau_pu)) * cos(nu_pu*t)
    return 
   END FUNCTION

 end module para


!*****************************************************************************80

subroutine rk4vec_test ( )
 use para

  implicit none

  integer ( kind = 4 ), parameter :: n = 2
  integer ( kind = 8 )  step_number

  real ( kind = 8 ), parameter :: dt = 0.005D+00
  real ( kind = 8 ) t0
  real ( kind = 8 ) t1
  real ( kind = 8 ), parameter :: tmax = 30000.0D+00

  complex ( kind = 8 ) u0(n)
  complex ( kind = 8 ) u1(n)

  real*8 :: fre, fre_min = 0.9, fre_max =1.1  
  integer*4  index_ini, index_final, data_number
  integer ierr_1, ierr_2, ierr_3, ierr_4, ierr_5

  complex(8), dimension(:), allocatable :: rho_ab_time, rho_ab_freq
  integer ( kind = 4 ), parameter :: fftw_estimate = 128
  integer ( kind = 4 ), parameter :: fftw_forward = -1
  integer*8 plan


  step_number = int(tmax/dt)
  write(*, *)  "the Total_Number   Time_Step       Time_END"
  write(*, '(2x,i12,4x,  g14.6, 2x, g14.6)')  step_number,  dt, tmax 

  !!---------- delete the previse data ------------------------------
    open(1, file='data.csv', status='unknown', iostat=ierr_1)
    if ( ierr_1 .eq. 0) then
      !    file opened successfully, delete 
       close(1, status='delete')
    endif

    open(2, file='fre_apt.csv', status='unknown', iostat=ierr_2)
    if ( ierr_2 .eq. 0) then
      !    file opened successfully, delete 
       close(2, status='delete')
    endif

  !! -------------  find the index which gives the frequency region we care about!
   index_ini = 0
  
   do i = 1, step_number   

     fre = 2.0 * pi * i/tmax
    
     if (fre < fre_min) then

     ! -----  mark the initial index of frequency
       index_ini = index_ini + 1
       Cycle

     else if (fre > fre_min .AND. fre < fre_max) then 
        cycle
     else 

       ! -----  mark the initial index of frequency
       index_final = i
       
       EXIT
     end if
   end do


  allocate(rho_ab_time(step_number))
  allocate(rho_ab_freq(step_number))


  open (unit = 4, file = 'data.csv', ACTION="write")
  open (unit = 5, file = 'fre_apt.csv', ACTION="write")

     t0 = 0.0D+00
     u0(1) = (1.0D+00, 0.0D+00)
     u0(2) = (0.0D+00, 0.0D+00)
    

     do  i = 1 , step_number

      rho_ab_time(i) = real(u0(1) * u0(2))

       write ( 4, * ) t0, ',' , Omega_pu(t0), ',', abs(u0(1))**2, ',' , abs(u0(2))**2 , ',', real(u0(1)*u0(2))
   
       t1 = t0 + dt
       call rk4vec ( t0, n, u0, dt, rk4vec_test_f, u1 )

       !
       !  ----------  Shift the data to prepare for another step.
       !
        t0 = t1
        u0(1:n) = u1(1:n)

     end do

    call dfftw_plan_dft_1d(plan,step_number,rho_ab_time,rho_ab_freq,fftw_forward,fftw_estimate)
    call dfftw_execute_dft(plan,rho_ab_time,rho_ab_freq)

    do i = index_ini + 1, index_final - 1
        write(5,*) 2.0 * pi * i/tmax, ',', abs(rho_ab_freq(i)) 
    end do

    call dfftw_destroy_plan(plan)


  close(4)

  return

 contains

 subroutine rk4vec_test_f ( t, n, u, uprime )

   !*****************************************************************************80
   use para
   implicit none
   integer ( kind = 4 ) n
   real ( kind = 8 ) t

   real(8) ham(2,2)
   real ( kind = 8 ) vector(4)

   complex ( kind = 8 ) u(n)
   complex ( kind = 8 ) uprime(n)


   !!! ----  please note that the max number of wards is 132.
   ham(1,1) = 1.0 ;       ham(1,2) = Omega_pu(t);
   ham(2,1) = ham(1,2);                       ham(2,2) = 0 ;                         

   uprime =  -imag*matmul(ham,u)
 
   return
   
 end subroutine rk4vec_test_f


 subroutine rk4vec ( t0, m, u0, dt, f, u )


        implicit none

        integer ( kind = 4 ) m

        real ( kind = 8 ) dt
        external f
        complex ( kind = 8 ) f0(m)
        complex ( kind = 8 ) f1(m)
        complex ( kind = 8 ) f2(m)
        complex ( kind = 8 ) f3(m)

        real ( kind = 8 ) t0
        real ( kind = 8 ) t1
        real ( kind = 8 ) t2
        real ( kind = 8 ) t3

        complex ( kind = 8 ) u(m)
        complex ( kind = 8 ) u0(m)
        complex ( kind = 8 ) u1(m)
        complex ( kind = 8 ) u2(m)
        complex ( kind = 8 ) u3(m)
       !
       !  Get four sample values of the derivative.
       !
        call f ( t0, m, u0, f0 )

        t1 = t0 + dt / 2.0D+00
        u1(1:m) = u0(1:m) + dt * f0(1:m) / 2.0D+00
        call f ( t1, m, u1, f1 )

        t2 = t0 + dt / 2.0D+00
        u2(1:m) = u0(1:m) + dt * f1(1:m) / 2.0D+00
        call f ( t2, m, u2, f2 )

        t3 = t0 + dt
        u3(1:m) = u0(1:m) + dt * f2(1:m)
        call f ( t3, m, u3, f3 )
       !
       !  Combine them to estimate the solution U at time T1.
       !
        u(1:m) = u0(1:m) + ( dt / 6.0D+00 ) * ( &
                       f0(1:m) &
           + 2.0D+00 * f1(1:m) &
           + 2.0D+00 * f2(1:m) &
           +           f3(1:m) )

        return
 end subroutine rk4vec



end subroutine rk4vec_test


!*****************************************************************************80
program main

  implicit none

  call rk4vec_test ( )


  stop
end program main

