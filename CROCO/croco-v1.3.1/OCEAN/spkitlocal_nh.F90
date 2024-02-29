#include "cppdefs.h"
#ifdef NBQ
!------------------------------------------------------------------------------
!                               NHOMS
!                Non Hydrostatic Ocean Modeling System      
!------------------------------------------------------------------------------
!
!> @note <a href="http://poc.obs-mip.fr/auclair/WOcean.fr/SNH/index_snh_home.htm"> Main web documentation </a>
!
! DESCRIPTION: 
!
!> @brief <a href="http://poc.obs-mip.fr/auclair/WOcean.fr/SNH/Restricted/NH-NBQ/Sources/Images/png/Algebrique_Gestion.png">
!! Sparskit library
!! </a> : Compressed Sparse Row format (CSR).
!> @details Additional documentation:
!! - <a href="http://www-users.cs.umn.edu/~saad/software/SPARSKIT/paper.ps">Youcef Saad</a>
!! - <a href="http://people.sc.fsu.edu/~jburkardt/f_src/sparsekit/sparsekit.html">Sparskit</a>
! REVISION HISTORY:
!
!> @authors
!> @date 2015 January
!> @todo
!
!------------------------------------------------------------------------------
subroutine amubdg ( nrow, ncol, ncolb, ja, ia, jb, ib, ndegr, nnz, iw )

!*****************************************************************************80
!
!> AMUBDG gets the number of nonzero elements in each row of A * B.
!
!> Discussion:
!!
!!   The routine also computes the total number of nonzero elements in A * B.
!!
!!   Method: A' * A = sum [over i = 1, nrow]  a(i)^T a(i)
!!   where a(i) = i-th row of  A.  We must be careful not to add  the
!!   elements already accounted for.
!
!> Modified:
!!
!!   07 January 2004
!
!> Author:
!!
!!   Youcef Saad
!
!> Parameters:
!!
!!   Input, integer NROW, the row dimension of the matrix A.
!!
!!   Input, integer NCOL, the column dimension of the matrix A,
!!   (and the row dimension of B).
!!
!!   Input, integer NCOLB, the column dimension of the matrix B.
!!
!!   Input, ja, ia= row structure of input matrix A: ja = column indices of
!!   the nonzero elements of A stored by rows.
!!   ia = pointer to beginning of each row in ja.
!!
!!   Input, jb, ib, the row structure of input matrix B: jb = column indices of
!!   the nonzero elements of A stored by rows.
!!   ib is a pointer to beginning of each row in jb.
!!
!!   Output, integer NDEGR(NROW), contains the degrees (the number of
!!   nonzeros in each row of the matrix A * B.
!!
!!   Output, integer NNZ, the number of nonzero elements found in A * B.
!!
!!  Workspace, integer IW(NCOLB).
!!
  implicit none

  integer ncol
  integer ncolb
  integer nrow

  integer ia(*)
  integer ib(*)
  integer ii
  integer iw(*)
  integer j
  integer ja(*)
  integer jb(*)
  integer jc
  integer jr
  integer k
  integer last
  integer ldg
  integer ndegr(*)
  integer nnz

  iw(1:ncolb) = 0
  ndegr(1:nrow) = 0

  do ii = 1, nrow
!
!  For each row of A.
!
    ldg = 0
!
!  End-of-linked list.
!
    last = -1

    do j = ia(ii), ia(ii+1)-1
!
!  Row number to be added.
!
        jr = ja(j)

        do k = ib(jr), ib(jr+1)-1
           jc = jb(k)
!
!  Add one element to the linked list.
!
           if ( iw(jc) == 0 ) then
              ldg = ldg + 1
              iw(jc) = last
              last = jc
           end if

         end do

    end do

    ndegr(ii) = ldg
!
!  Reset IW to zero.
!
    do k = 1, ldg
      j = iw(last)
      iw(last) = 0
      last = j
     end do

  end do

  nnz = sum ( ndegr(1:nrow) )

  return
end
subroutine amux ( n, x, y, a, ja, ia )
!subroutine amux ( n, x, y, a, ja, nl )

!*****************************************************************************80
!
!! AMUX multiplies a CSR matrix A times a vector.
!
!  Discussion:
!
!    This routine multiplies a matrix by a vector using the dot product form.
!    Matrix A is stored in compressed sparse row storage.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
! 
!  Parameters:
!
!    Input, integer N, the row dimension of the matrix.
!
!    Input, real X(*), and array of length equal to the column dimension
!    of A.
!
!    Input, real A(*), integer JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!    Output, real Y(N), the product A * X.
!

!  use module_parallele                                           ! #MPI#

  use module_nh
  implicit none
! include 'mkl_spblas.h'
  integer n

  real ( kind = 8 ) a(*)
  integer i,j
  integer ia(*)
  integer ja(*)
  integer k
! integer nl
! character(len=1) :: TRANS='N'
  real ( kind = 8 ) t
  real ( kind = 8 ) x(*)
  real ( kind = 8 ) y(*)

  do i = 1, n
!
!  Compute the inner product of row I with vector X.
!
    t = 0.0D+00
     do k = ia(i), ia(i+1)-1
!    do k = 1 + nl*(i-1), nl*i
       t = t + a(k) * x(ja(k))
     end do
     y(i) = t
  end do

  return
end

subroutine amuxa ( n, x, y, a, ja, ia )
!subroutine amux ( n, x, y, a, ja, nl )

!*****************************************************************************80
!
!! AMUX multiplies a CSR matrix A times a vector.
!
!  Discussion:
!
!    This routine multiplies a matrix by a vector using the dot product form.
!    Matrix A is stored in compressed sparse row storage.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
! 
!  Parameters:
!
!    Input, integer N, the row dimension of the matrix.
!
!    Input, real X(*), and array of length equal to the column dimension
!    of A.
!
!    Input, real A(*), integer JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!    Output, real Y(N), the product A * X.
!

!  use module_parallele                                           ! #MPI#

  use module_nh
  implicit none
! include 'mkl_spblas.h'
  integer n

  real ( kind = 8 ) a(*)
  integer i,j
  integer ia(*)
  integer ja(*)
  integer k
! integer nl
! character(len=1) :: TRANS='N'
  real ( kind = 8 ) t
  real ( kind = 8 ) x(*)
  real ( kind = 8 ) y(*)

  do i = 1, n
!
!  Compute the inner product of row I with vector X.
!
    t = 0.0D+00
     do k = ia(i), ia(i+1)-1
!    do k = 1 + nl*(i-1), nl*i
       t = t + a(k) * x(ja(k))
     end do
     y(i) = y(i)+t
  end do

  return
end

!*****************************************************************************80
!                amux pour MOM
!*****************************************************************************80
subroutine amux_s ( n , n1, x, y, a, ja, ia )

!  use module_parallele                                           ! #MPI#

  use module_nh
  implicit none
! include 'mkl_spblas.h'
  integer n,n1

  real ( kind = 8 ) a(*)
  integer i
  integer ia(*)
  integer ja(*)
  integer k
  character(len=1) :: TRANS='N'
  real ( kind = 8 ) t
  real ( kind = 8 ) x(*)
  real ( kind = 8 ) y(*)
  double precision :: cput1, cput2

  !cput1=mpi_wtime()
  
  do i = 1, n1
    t = 0.0D+00
    do k = ia(i), ia(i+1)-1
      t = t + a(k) * x(ja(k))
    end do
    y(i) = t
  end do

  do i = n1+1, n
    y(i) = a(ia(i)) * x(ja(ia(i))) + a(ia(i)+1) * x(ja(ia(i)+1)) 
  end do


  return
  end
subroutine amub2 ( nrow, ncol, job, a, ja, ia, b, jb, ib, c, jc, ic, nzmax, &
  iw, ierr )

!*****************************************************************************
!
!! AMUB performs the matrix product C = A * B.
!
!  Discussion:
!
!    The column dimension of B is not needed.
!
!  Modified:
!
!    08 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer NROW, the row dimension of the matrix.
!
!    Input, integer NCOL, the column dimension of the matrix.
!
!    Input, integer JOB, job indicator.  When JOB = 0, only the structure
!    is computed, that is, the arrays JC and IC, but the real values
!    are ignored.
!
!    Input, real A(*), integer JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!    Input, b, jb, ib, matrix B in compressed sparse row format.
!
!    Input, integer NZMAX, the length of the arrays c and jc.
!    The routine will stop if the result matrix C  has a number
!    of elements that exceeds exceeds NZMAX.
!
! on return:
!
! c,
! jc,
! ic    = resulting matrix C in compressed sparse row sparse format.
!
! ierr      = integer. serving as error message.
!         ierr = 0 means normal return,
!         ierr > 0 means that amub stopped while computing the
!         i-th row  of C with i = ierr, because the number
!         of elements in C exceeds nzmax.
!
! work arrays:
!
!  iw      = integer work array of length equal to the number of
!         columns in A.
!
  implicit none

  integer ncol
  integer nrow
  integer nzmax

  real ( kind = 8 ) a(*)
  real ( kind = 8 ) b(*)
  real ( kind = 8 ) c(nzmax)
!  integer ia(nrow+1)
!  integer ib(ncol+1)
!  integer ic(ncol+1)
  integer ia(*)
  integer ib(*)
  integer ic(*)
  integer ierr
  integer ii
  integer iw(*)
  integer ja(*)
  integer jb(*)
  integer jc(nzmax)
  integer jcol
  integer jj
  integer job
  integer jpos
  integer k
  integer ka
  integer kb
  integer len
  real ( kind = 8 ) scal
  logical values

  values = ( job /= 0 )
  len = 0
  ic(1) = 1
  ierr = 0
!
!  Initialize IW.
!
   iw(1:ncol) = 0

   do ii = 1, nrow
!
!  Row I.
!
    do ka = ia(ii), ia(ii+1)-1

      if ( values ) then
        scal = a(ka)
      end if

      jj = ja(ka)

      do kb = ib(jj), ib(jj+1)-1

           jcol = jb(kb)
           jpos = iw(jcol)

           if ( jpos == 0 ) then
              len = len + 1
              if ( nzmax < len ) then
                 ierr = ii
                 return
              end if
              jc(len) = jcol
              iw(jcol)= len
              if ( values ) then
                c(len) = scal * b(kb)
              end if
           else
              if ( values ) then
                c(jpos) = c(jpos) + scal * b(kb)
              end if
           end if

         end do

    end do

    do k = ic(ii), len
      iw(jc(k)) = 0
    end do

    ic(ii+1) = len + 1

  end do

  return
  end
subroutine amub3 ( nrow, ncol, job, alpha, a, ja, ia, b, jb, ib, c, jc, ic, nzmax, &
  iw, ierr )

!*****************************************************************************
!
!! AMUB performs the matrix product C = A * B.
!
!  Discussion:
!
!    The column dimension of B is not needed.
!
!  Modified:
!
!    08 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer NROW, the row dimension of the matrix.
!
!    Input, integer NCOL, the column dimension of the matrix.
!
!    Input, integer JOB, job indicator.  When JOB = 0, only the structure
!    is computed, that is, the arrays JC and IC, but the real values
!    are ignored.
!
!    Input, real A(*), integer JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!    Input, b, jb, ib, matrix B in compressed sparse row format.
!
!    Input, integer NZMAX, the length of the arrays c and jc.
!    The routine will stop if the result matrix C  has a number
!    of elements that exceeds exceeds NZMAX.
!
! on return:
!
! c,
! jc,
! ic    = resulting matrix C in compressed sparse row sparse format.
!
! ierr      = integer. serving as error message.
!         ierr = 0 means normal return,
!         ierr > 0 means that amub stopped while computing the
!         i-th row  of C with i = ierr, because the number
!         of elements in C exceeds nzmax.
!
! work arrays:
!
!  iw      = integer work array of length equal to the number of
!         columns in A.
!
  implicit none

  integer ncol
  integer nrow
  integer nzmax

  real ( kind = 8 ) alpha
  real ( kind = 8 ) a(*)
  real ( kind = 8 ) b(*)
  real ( kind = 8 ) c(nzmax)
!  integer ia(nrow+1)
!  integer ib(ncol+1)
!  integer ic(ncol+1)
  integer ia(*)
  integer ib(*)
  integer ic(*)
  integer ierr
  integer ii
  integer iw(*)
  integer ja(*)
  integer jb(*)
  integer jc(nzmax)
  integer jcol
  integer jj
  integer job
  integer jpos
  integer k
  integer ka
  integer kb
  integer len
  real ( kind = 8 ) scal
  logical values

  values = ( job /= 0 )
  len = 0
  ic(1) = 1
  ierr = 0
!
!  Initialize IW.
!
   iw(1:ncol) = 0

   do ii = 1, nrow
!
!  Row I.
!
    do ka = ia(ii), ia(ii+1)-1

      if ( values ) then
        scal = a(ka)*alpha
      end if

      jj = ja(ka)

      do kb = ib(jj), ib(jj+1)-1

           jcol = jb(kb)
           jpos = iw(jcol)

           if ( jpos == 0 ) then
              len = len + 1
              if ( nzmax < len ) then
                 ierr = ii
                 return
              end if
              jc(len) = jcol
              iw(jcol)= len
              if ( values ) then
                c(len) = scal * b(kb)
              end if
           else
              c(jpos) = c(jpos) + scal * b(kb)
           end if

         end do

    end do

    do k = ic(ii), len
      iw(jc(k)) = 0
    end do

    ic(ii+1) = len + 1

  end do

  return
  end
subroutine aplb ( nrow, ncol, job, a, ja, ia, b, jb, ib, c, jc, ic, nzmax, &
  iw, ierr )

!*****************************************************************************80
!
!! APLB performs the CSR matrix sum C = A + B.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NROW, the row dimension of A and B.
!
!    Input, integer ( kind = 4 ) NCOL, the column dimension of A and B.
!
!    Input, integer ( kind = 4 ) JOB.  When JOB = 0, only the structure
!    (i.e. the arrays jc, ic) is computed and the
!    real values are ignored.
!
!    Input, real A(*), integer ( kind = 4 ) JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
! b,
! jb,
! ib      =  Matrix B in compressed sparse row format.
!
! nzmax      = integer ( kind = 4 ). The  length of the arrays c and jc.
!         amub will stop if the result matrix C  has a number
!         of elements that exceeds exceeds nzmax. See ierr.
!
! on return:
!
! c,
! jc,
! ic      = resulting matrix C in compressed sparse row sparse format.
!
! ierr      = integer ( kind = 4 ). serving as error message.
!         ierr = 0 means normal return,
!         ierr > 0 means that amub stopped while computing the
!         i-th row  of C with i = ierr, because the number
!         of elements in C exceeds nzmax.
!
! work arrays:
!
! iw      = integer ( kind = 4 ) work array of length equal to the number of
!         columns in A.
!
  implicit none

  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nrow

  real ( kind = 8 ) a(*)
  real ( kind = 8 ) b(*)
  real ( kind = 8 ) c(*)
  integer ( kind = 4 ) ia(nrow+1)
  integer ( kind = 4 ) ib(nrow+1)
  integer ( kind = 4 ) ic(nrow+1)
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) iw(ncol)
  integer ( kind = 4 ) ja(*)
  integer ( kind = 4 ) jb(*)
  integer ( kind = 4 ) jc(*)
  integer ( kind = 4 ) jcol
  integer ( kind = 4 ) job
  integer ( kind = 4 ) jpos
  integer ( kind = 4 ) k
  integer ( kind = 4 ) ka
  integer ( kind = 4 ) kb
  integer ( kind = 4 ) len
  integer ( kind = 4 ) nzmax
  logical values

  values = ( job /= 0 )
  ierr = 0
  len = 0
  ic(1) = 1
  iw(1:ncol) = 0

  do ii = 1, nrow
!
!  Row I.
!
     do ka = ia(ii), ia(ii+1)-1

        len = len + 1
        jcol = ja(ka)

        if ( nzmax < len ) then
          ierr = ii
          return
        end if

        jc(len) = jcol
        if ( values ) then
          c(len) = a(ka)
        end if
        iw(jcol) = len
     end do

     do kb = ib(ii), ib(ii+1)-1

        jcol = jb(kb)
        jpos = iw(jcol)

        if ( jpos == 0 ) then

           len = len + 1

           if ( nzmax < len ) then
             ierr = ii
             return
           end if

           jc(len) = jcol
           if ( values ) then
             c(len) = b(kb)
           end if
           iw(jcol)= len
        else
           if ( values ) then
             c(jpos) = c(jpos) + b(kb)
           end if
        end if

     end do

     do k = ic(ii), len
       iw(jc(k)) = 0
     end do

     ic(ii+1) = len+1
  end do

  return
end
subroutine aplb1 ( nrow, ncol, job, a, ja, ia, b, jb, ib, c, jc, ic, &
  nzmax, ierr )

!*****************************************************************************80
!
!! APLB1 performs the sum C = A + B for sorted CSR matrices.
!
!  Discussion:
!
!    The difference between this routine and APLB is that here the
!    resulting matrix is such that the elements of each row are sorted,
!    with increasing column indices in each row, provided the original
!    matrices are sorted in the same way.
!
!    This routine will not work if either of the two input matrices is
!    not sorted.
!
!  Modified:
!
!    11 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NROW, the row dimension of A and B.
!
!    Input, integer ( kind = 4 ) NCOL, the column dimension of A and B.
!
!    Input, integer ( kind = 4 ) JOB.  When JOB = 0, only the structure
!    (i.e. the arrays jc, ic) is computed and the
!    real values are ignored.
!
!    Input, real A(*), integer ( kind = 4 ) JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format with entries sorted.
!
! b,
! jb,
! ib      =  Matrix B in compressed sparse row format with entries sorted
!        ascendly in each row
!
! nzmax      = integer ( kind = 4 ). The  length of the arrays c and jc.
!         amub will stop if the result matrix C  has a number
!         of elements that exceeds exceeds nzmax. See ierr.
!
! on return:
!
! c,
! jc,
! ic      = resulting matrix C in compressed sparse row sparse format
!         with entries sorted ascendly in each row.
!
! ierr      = integer ( kind = 4 ). serving as error message.
!         ierr = 0 means normal return,
!         ierr > 0 means that amub stopped while computing the
!         i-th row  of C with i = ierr, because the number
!         of elements in C exceeds nzmax.
!
  implicit none

  integer ( kind = 4 ) nrow

  real ( kind = 8 ) a(*)
  real ( kind = 8 ) b(*)
  real ( kind = 8 ) c(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia(nrow+1)
  integer ( kind = 4 ) ib(nrow+1)
  integer ( kind = 4 ) ic(nrow+1)
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) ja(*)
  integer ( kind = 4 ) jb(*)
  integer ( kind = 4 ) jc(*)
  integer ( kind = 4 ) job
  integer ( kind = 4 ) ka
  integer ( kind = 4 ) kamax
  integer ( kind = 4 ) kb
  integer ( kind = 4 ) kbmax
  integer ( kind = 4 ) kc
  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nzmax
  logical values

  values = ( job /= 0 )
  ierr = 0
  kc = 1
  ic(1) = kc

  do i = 1, nrow

    ka = ia(i)
    kb = ib(i)
    kamax = ia(i+1) - 1
    kbmax = ib(i+1) - 1

    do

      if ( ka <= kamax ) then
        j1 = ja(ka)
      else
        j1 = ncol + 1
      end if

      if ( kb <= kbmax ) then
        j2 = jb(kb)
      else
        j2 = ncol + 1
      end if
!
!  Three cases
!
      if ( j1 == j2 ) then
        if ( values ) then
          c(kc) = a(ka) + b(kb)
        end if
        jc(kc) = j1
        ka = ka + 1
        kb = kb + 1
        kc = kc + 1
      else if ( j1 < j2 ) then
        jc(kc) = j1
        if ( values ) then
          c(kc) = a(ka)
        end if
        ka = ka + 1
        kc = kc + 1
      else if ( j2 < j1 ) then
        jc(kc) = j2
        if ( values ) then
          c(kc) = b(kb)
        end if
        kb = kb + 1
        kc = kc + 1
      end if

      if ( nzmax < kc ) then
        ierr = i
        return
      end if

      if ( kamax < ka .and. kbmax < kb ) then
        exit
      end if

     end do

     ic(i+1) = kc

  end do

  return
end
subroutine aplsca1 ( nrow, a, ja, ia, iw )

!*****************************************************************************80
!
!! APLSCA adds 1. to the diagonal entries of a sparse matrix A :=A + s I
!
!  Discussion:
!
!    The column dimension of A is not needed.
!
!    important: the matrix A may be expanded slightly to allow for
!    additions of nonzero elements to previously nonexisting diagonals.
!    The is no checking as to whether there is enough space appended
!    to the arrays a and ja. if not sure allow for n additional
!    elements.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NROW, the row dimension of the matrix.
!
!    Input, real A(*), integer ( kind = 4 ) JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
! on return:
!
!
! a,
! ja,
! ia      = matrix A with diagonal elements shifted (or created).
!
! iw    = integer ( kind = 4 ) work array of length n. On return iw will
!         contain  the positions of the diagonal entries in the
!         output matrix. (i.e., a(iw(k)), ja(iw(k)), k = 1,...n,
!         are the values/column indices of the diagonal elements
!         of the output matrix. ).
!
  implicit none

  integer ( kind = 4 ) nrow

  real ( kind = 8 ) a(*)
  integer ( kind = 4 ) ia(nrow+1)
  integer ( kind = 4 ) icount
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) iw(*)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ja(*)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) k2
  integer ( kind = 4 ) ko
  logical test

  call diapos ( nrow, ja, ia, iw )
  icount = 0

  do j = 1, nrow

     if ( iw(j) == 0 ) then
        icount = icount + 1
     else
        a(iw(j)) = a(iw(j)) + 1. 
     end if

  end do
!
!  If no diagonal elements to insert in data structure, return.
!
  if ( icount == 0 ) then
    return
  end if
!
!  Shift the nonzero elements if needed, to allow for created
!  diagonal elements.
!
  ko = ia(nrow+1) + icount
!
!  Copy rows backward.
!
  do ii = nrow, 1, -1
!
!  Go through row II.
!
     k1 = ia(ii)
     k2 = ia(ii+1) - 1
     ia(ii+1) = ko
     test = ( iw(ii) == 0 )

     do k = k2, k1, -1

        j = ja(k)

        if ( test .and. j < ii ) then
           test = .false.
           ko = ko - 1
           a(ko) = 1.
           ja(ko) = ii
           iw(ii) = ko
        end if

        ko = ko - 1
        a(ko) = a(k)
        ja(ko) = j

    end do
!
!  The diagonal element has not been added yet.
!
     if ( test ) then
        ko = ko - 1
        a(ko)  = 1. 
        ja(ko) = ii
        iw(ii) = ko
     end if

  end do

  ia(1) = ko

  return
end
subroutine diapos ( n, ja, ia, idiag )

!*****************************************************************************80
!
!! DIAPOS returns the positions of the diagonal elements of a sparse matrix.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the row dimension of the matrix.
!
!    Input, JA(*), IA(N+1), the matrix information, (but no values) 
!    in CSR Compressed Sparse Row format.
!
!    Output, integer ( kind = 4 ) IDIAG(N); the I-th entry of IDIAG points to the 
!    diagonal element A(I,I) in the arrays A and JA.  That is,
!    A(IDIAG(I)) = element A(I,I) of matrix A.  If no diagonal element 
!    is found, the entry is set to 0.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia(n+1)
  integer ( kind = 4 ) idiag(n)
  integer ( kind = 4 ) ja(*)
  integer ( kind = 4 ) k

  idiag(1:n) = 0
!
!  Sweep through the data structure.
!
  do i = 1, n
    do k = ia(i), ia(i+1) -1
      if ( ja(k) == i ) then
        idiag(i) = k
      end if
    end do
  end do

  return
end


#ifdef FOR_TESTS

subroutine amux_1 ( n, x, y, a, ja, ia )
!subroutine amux ( n, x, y, a, ja, nl )

!*****************************************************************************80
!
!! AMUX multiplies a CSR matrix A times a vector.
!
!  Discussion:
!
!    This routine multiplies a matrix by a vector using the dot product form.
!    Matrix A is stored in compressed sparse row storage.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
! 
!  Parameters:
!
!    Input, integer N, the row dimension of the matrix.
!
!    Input, real X(*), and array of length equal to the column dimension
!    of A.
!
!    Input, real A(*), integer JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!    Output, real Y(N), the product A * X.
!

!  use module_parallele                                           ! #MPI#

  use module_nh
  use module_nbq
  implicit none
! include 'mkl_spblas.h'
  integer n

  real ( kind = 8 ) a(*)
  integer i,j
  integer ia(*)
  integer ja(*)
  integer k
! integer nl
! character(len=1) :: TRANS='N'
  real ( kind = 8 ) x(*)
  real ( kind = 8 ) y(*) 

  namux1_nbq=0
  do i = 1, n
!
!  Compute the inner product of row I with vector X.
!
   ! t = 0.0D+00
     do k = ia(i), ia(i+1)-1
!    do k = 1 + nl*(i-1), nl*i
   !    t = t + a(k) * x(ja(k))
       namux1_nbq = namux1_nbq+1
       amux1a_nbq(namux1_nbq)=i
       amux1b_nbq(namux1_nbq)=k
       amux1c_nbq(namux1_nbq)=ja(k)     
     end do
   !  y(i) = t
  end do

  return
end

subroutine amux_2 ( n, x, y, a )
!subroutine amux ( n, x, y, a, ja, nl )

!*****************************************************************************80
!
!! AMUX multiplies a CSR matrix A times a vector.
!
!  Discussion:
!
!    This routine multiplies a matrix by a vector using the dot product form.
!    Matrix A is stored in compressed sparse row storage.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
! 
!  Parameters:
!
!    Input, integer N, the row dimension of the matrix.
!
!    Input, real X(*), and array of length equal to the column dimension
!    of A.
!
!    Input, real A(*), integer JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!    Output, real Y(N), the product A * X.
!

!  use module_parallele                                           ! #MPI#

  use module_nh
  use module_nbq
  implicit none
! include 'mkl_spblas.h'
  integer n

  real ( kind = 8 ) a(*)
  integer i,j
  integer k
! integer nl
! character(len=1) :: TRANS='N'
  real ( kind = 8 ) x(*)
  real ( kind = 8 ) y(*) 

  integer icall

  y(1:n)=0.
  do l_nbq=1,namux1_nbq
     y(amux1a_nbq(l_nbq)) = y(amux1a_nbq(l_nbq)) +a(amux1b_nbq(l_nbq))*x(amux1c_nbq(l_nbq))
  enddo

  return

end
#endif

#endif

