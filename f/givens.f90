module givens

  use nrtype_m

  implicit none

  private

  public :: left_multiply, right_multiply, rotate

contains

subroutine left_multiply(a, i, j, c, s)
  complex(DPC), intent(inout) :: a(:, :, :)
  integer, intent(in) :: i
  integer, intent(in) :: j
  real(DP), intent(in) :: c
  complex(DPC), intent(in) :: s
  integer :: nmat
  integer :: nelem
  integer :: k

  nmat = size(a, dim=1)
  nelem = size(a, dim=2)

  do k=1,nelem
    call zrot(nmat, a(:, i, k), 1, a(:, j, k), 1, c, -s)
  enddo

end subroutine left_multiply

subroutine right_multiply(a, i, j, c, s)
  complex(DPC), intent(inout) :: a(:, :, :)
  integer, intent(in) :: i
  integer, intent(in) :: j
  real(DP), intent(in) :: c
  complex(DPC), intent(in) :: s
  integer :: nmat
  integer :: nelem
  integer :: k

  nmat = size(a, dim=1)
  nelem = size(a, dim=2)

  do k=1,nelem
    call zrot(nmat, a(:, k, j), 1, a(:, k, i), 1, c, s)
  enddo

end subroutine right_multiply

subroutine rotate(a, i, j, c, s)
  complex(DPC), intent(inout) :: a(:, :, :)
  integer, intent(in) :: i
  integer, intent(in) :: j
  real(DP), intent(in) :: c
  complex(DPC), intent(in) :: s
  integer :: nmat
  integer :: nelem
  integer :: k

  nmat = size(a, dim=1)
  nelem = size(a, dim=2)

  do k=1,nelem
    call zrot(nmat, a(:, i, k), 1, a(:, j, k), 1, c, -s)
  enddo
  do k=1,nelem
    call zrot(nmat, a(:, k, j), 1, a(:, k, i), 1, c, s)
  enddo

end subroutine rotate

end module
