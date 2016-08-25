module nrtype_m

  implicit none

  public

  integer, parameter :: DP = kind(1.0d0)
  integer, parameter :: DPC = kind((1.0d0,1.0d0))
  complex(DPC), parameter :: CMPLX_0 = (0.0_dp,0.0_dp)
  complex(DPC), parameter :: CMPLX_1 = (1.0_dp,0.0_dp)
  complex(DPC), parameter :: CMPLX_I = (0.0_dp,1.0_dp)

end module nrtype_m
