subroutine random_initialize ( seed )
!
!*******************************************************************************
!
!! RANDOM_INITIALIZE initializes the FORTRAN 90 random number seed.
!
!
!  Discussion:
!
!    If you don't initialize the random number generator, its behavior
!    is not specified.  If you initialize it simply by:
!
!      call random_seed
!
!    its behavior is not specified.  On the DEC ALPHA, if that's all you
!    do, the same random number sequence is returned.  In order to actually
!    try to scramble up the random number generator a bit, this routine
!    goes through the tedious process of getting the size of the random
!    number seed, making up values based on the current time, and setting
!    the random number seed.
!
!  Modified:
!
!    19 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer(4) SEED.
!    If SEED is zero on input, then you're asking this routine to come up
!    with a seed value, which is returned as output.
!    If SEED is nonzero on input, then you're asking this routine to
!    use the input value of SEED to initialize the random number generator,
!    and SEED is not changed on output.
!
  USE WPRECISION
  implicit none
!
  integer(4) count
  integer(4) count_max
  integer(4) count_rate
  logical, parameter :: debug = .false.
  integer(4) i
  integer(4) seed
  integer(4), allocatable :: seed_vector(:)
  integer(4) seed_size
  REAL(WP) t
!
!  Initialize the random number seed.
!
  call random_seed
!
!  Determine the size of the random number seed.
!
  call random_seed ( size = seed_size )
!
!  Allocate a seed of the right size.
!
  allocate ( seed_vector(seed_size) )

  if ( seed /= 0 ) then

    if ( debug ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RANDOM_INITIALIZE'
      write ( *, '(a,i20)' ) '  Initialize RANDOM_NUMBER, user SEED = ', seed
    end if

  else

    call system_clock ( count, count_rate, count_max )

    seed = count

    if ( debug ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RANDOM_INITIALIZE'
      write ( *, '(a,i20)' ) '  Initialize RANDOM_NUMBER, arbitrary SEED = ', &
        seed
    end if

  end if
!
!  Now set the seed.
!
  seed_vector(1:seed_size) = seed

  call random_seed ( put = seed_vector(1:seed_size) )
!
!  Free up the seed space.
!
  deallocate ( seed_vector )
!
!  Call the random number routine a bunch of times.
!
  do i = 1, 100
    call random_number ( harvest = t )
  end do

  return
end


subroutine rvec_random ( alo, ahi, n, a )
!
!*******************************************************************************
!
!! RVEC_RANDOM returns a random REAL(WP) vector in a given range.
!
!
!  Modified:
!
!    04 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, REAL(WP) ALO, AHI, the range allowed for the entries.
!
!    Input, integer(4) N, the number of entries in the vector.
!
!    Output, REAL(WP) A(N), the vector of randomly chosen values.
!
  USE WPRECISION
  implicit none
!
  integer(4) n
!
  REAL(WP) a(n)
  REAL(WP) ahi
  REAL(WP) alo
  integer(4) i
!
  do i = 1, n
    call r_random ( alo, ahi, a(i) )
  end do

  return
end


subroutine r_random ( rlo, rhi, r )
!
!*******************************************************************************
!
!! R_RANDOM returns a random REAL(WP) in a given range.
!
!
!  Modified:
!
!    06 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, REAL(WP) RLO, RHI, the minimum and maximum values.
!
!    Output, REAL(WP) R, the randomly chosen value.
!
  USE WPRECISION
  implicit none
!
  REAL(WP) r
  REAL(WP) rhi
  REAL(WP) rlo
  REAL(WP) t
!
!  Pick T, a random number in (0,1).
!
  call random_number ( harvest = t )
!
!  Set R in ( RLO, RHI ).
!
  r = ( 1.0E+00 - t ) * rlo + t * rhi

  return
end

