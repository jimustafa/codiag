module codiag

  use nrtype_m

  implicit none

  private

  public :: diag2_ii_Abc, offdiag2_ikki_Abc

contains

subroutine diag2_ii_Abc(m, n, i, j, A, b, c)
  complex(DPC), intent(inout) :: m(:, :, :)
  integer, intent(in) :: n
  integer, intent(in) :: i
  integer, intent(in) :: j
  real(DP), intent(inout) :: A(3,3)
  real(DP), intent(inout) :: b(3)
  real(DP), intent(inout) :: c
  integer :: ii
  complex(DPC) :: z(n, 3)
  complex(DPC) :: b0(n)
  complex(DPC) :: A0(3, 3)

  z = 0d0
  z(:, 1) = 1d0/2 * (m(:, i, i)-m(:, j, j))
  z(:, 2) = 1d0/2 * (-(m(:, i, j)+m(:, j, i)))
  z(:, 3) = 1d0/2 * CMPLX_I * (m(:, i, j)-m(:, j, i))

  c = 1./4 * sum(abs(m(:, i, i) + m(:, j, j))**2)

  b0 = conjg(m(:, i, i)+m(:, j, j))
  b(1) = real(sum(b0*z(:, 1)))
  b(2) = real(sum(b0*z(:, 2)))
  b(3) = real(sum(b0*z(:, 3)))

  A0 = 0d0
  call zgemm('T', 'N', 3, 3, n, CMPLX_1, z(:, :), n, conjg(z(:, :)), n, CMPLX_0, A0(:, :), 3)
  A = real(A0)

end subroutine diag2_ii_Abc

subroutine offdiag2_ikki_Abc(m, n, i, j, kidx, A, b, c)
  complex(DPC), intent(inout) :: m(:, :, :)
  integer, intent(in) :: n
  integer, intent(in) :: i
  integer, intent(in) :: j
  integer, intent(in) :: kidx(:)
  real(DP), intent(inout) :: A(3,3)
  real(DP), intent(inout) :: b(3)
  real(DP), intent(inout) :: c
  integer :: nk
  integer :: ii
  integer :: k
  real(DP):: aik2, aki2, ajk2, akj2
  real(DP) :: b1, b2, b3

  c = 0d0
  b = 0d0
  A = 0d0
  nk = size(kidx, dim=1)

  do ii=1,nk
    k = kidx(ii)

    aik2 = sum(abs(m(:, i, k))**2)
    ajk2 = sum(abs(m(:, j, k))**2)
    aki2 = sum(abs(m(:, k, i))**2)
    akj2 = sum(abs(m(:, k, j))**2)

    b1 = b1 + 1./2 * (aik2 - ajk2 + aki2 - akj2)
    b2 = b2 - sum(real(m(:, i, k)*conjg(m(:, j, k)) + m(:, k, i)*conjg(m(:, k, j))))
    b3 = b3 + sum(&
        + imag(m(:, j, k)) * real(m(:, i, k)) &
        - imag(m(:, i, k)) * real(m(:, j, k)) &
        - imag(m(:, k, j)) * real(m(:, k, i)) &
        + imag(m(:, k, i)) * real(m(:, k, j)) &
      )

    c = c + 1./2 * (aik2 + ajk2 + aki2 + akj2)
  enddo

  b(1) = b1
  b(2) = b2
  b(3) = b3

end subroutine offdiag2_ikki_Abc

end module
