! This file is part of multicharge.
! SPDX-Identifier: Apache-2.0
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

module multicharge_output
   use mctc_env, only : wp, timer_type, format_time
   use mctc_io, only : structure_type
   use mctc_io_constants, only : pi
   use mctc_io_convert, only : autoaa
   use multicharge_model, only : mchrg_model_type
   use multicharge_version, only : get_multicharge_version
   implicit none
   private

   public :: write_ascii_model, write_ascii_properties, write_ascii_results
   public :: json_results
   public :: print_cg_header, print_cg_convergence, print_cg_iteration, &
      & print_cg_final, write_direct_solver, print_direct_final


contains


!> Write the charge model parameters in tabular form
subroutine write_ascii_model(unit, mol, model, verbosity, timer)

   !> Formatted unit
   integer, intent(in) :: unit

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Electronegativity equilibration model
   class(mchrg_model_type), intent(in) :: model

   !> Verbosity level for output
   integer, intent(in) :: verbosity

   !> Timer for performance measurement
   real(wp), intent(in) :: timer

   integer :: isp
   real(wp), parameter :: sqrt2pi = sqrt(2.0_wp/pi)

   write(unit, '(a, ":")') "Charge model parameter"
   write(unit, '(54("-"))')
   write(unit, "(a4,5x,*(1x,a10))") "Z", "chi/Eh", "kcn_chi/Eh", "eta/Eh", "rad/AA"
   write(unit, '(54("-"))')
   do isp = 1, mol%nid
      write(unit, "(i4, 1x, a4, *(1x,f10.4))") &
         & mol%num(isp), mol%sym(isp), model%chi(isp), model%kcnchi(isp), &
         & model%eta(isp) + sqrt2pi/model%rad(isp), model%rad(isp) * autoaa
   end do
   write(unit, '(54("-"),/)')

   if (verbosity > 1) then
      write(unit, "(a, 1x, a)") "Model setup time :", format_time(timer)
      write(unit, "(a)") ""
   end if

end subroutine write_ascii_model

subroutine write_ascii_properties(unit, mol, model, cn, qvec)

   !> Unit for output
   integer, intent(in) :: unit

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Electronegativity equilibration model
   class(mchrg_model_type), intent(in) :: model

   !> Coordination numbers
   real(wp), intent(in) :: cn(:)

   !> Atomic partial charges
   real(wp), intent(in) :: qvec(:)

   integer :: iat, isp

   write(unit, '(54("-"))')
   write(unit, "(24x,a)") "Results"
   write(unit, '(54("-"))')
   write(unit, "(a)") ""
   write(unit, '(a,":")') "Electrostatic properties (in atomic units)"
   write(unit, '(54("-"))')
   write(unit, "(a10,1x,a4,5x,*(1x,a10))") "#", "Z", "CN", "q", "chi"
   write(unit, '(54("-"))')
   do iat = 1, mol%nat
      isp = mol%id(iat)
      write(unit, "(i10,1x,i4,1x,a4,*(1x,f10.4))") &
      & iat, mol%num(isp), mol%sym(isp), cn(iat), qvec(iat), &
      & model%chi(isp) - model%kcnchi(isp) * sqrt(cn(iat))
   end do
   write(unit, '(54("-"))')
   write(unit, "(a7,26x,f10.4)") "Σ", sum(qvec)
   write(unit, '(54("-"),/)')

end subroutine write_ascii_properties


!> Write the electrostatic energy, gradient, and charge derivatives
subroutine write_ascii_results(unit, mol, energy, gradient, sigma, dqdr, dqdL)

   !> Unit for output
   integer, intent(in) :: unit

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   real(wp), intent(in) :: energy(:)
   real(wp), intent(in), optional :: gradient(:, :)
   real(wp), intent(in), optional :: sigma(:, :)

   !> Derivative of the partial charges w.r.t. the Cartesian coordinates
   real(wp), intent(in), optional :: dqdr(:, :, :)

   !> Derivative of the partial charges w.r.t. strain deformations
   real(wp), intent(in), optional :: dqdL(:, :, :)

   integer :: iat, jat, isp, jsp
   logical :: grad, qgrad
   character(len=1), parameter :: comp(3) = ["x", "y", "z"]

   grad = present(gradient) .and. present(sigma)
   qgrad = present(dqdr) .and. present(dqdL)

   write(unit, '(a,":", t25, es20.13, 1x, a)') &
      & "Electrostatic energy", sum(energy), "Eh"
   write(unit, "(a)")
   if (grad) then
      write(unit, '(a,":", t25, es20.13, 1x, a)') &
      & "Energy gradient norm", norm2(gradient), "Eh/a0"
      write(unit, '(54("-"))')
      write(unit, "(a10,1x,a4,5x,*(1x,a10))") "#", "Z", "dE/dx", "dE/dy", "dE/dz"
      write(unit, '(54("-"))')
      do iat = 1, mol%nat
         isp = mol%id(iat)
         write(unit, "(i10,1x,i4,1x,a4,*(es11.3))") &
         & iat, mol%num(isp), mol%sym(isp), gradient(:, iat)
      end do
      write(unit, '(54("-"))')
      write(unit, "(a)")

      write(unit, '(a,":")') &
      & "Energy virial"
      write(unit, '(50("-"))')
      write(unit, "(a15,1x,*(1x,a10))") "component", "x", "y", "z"
      write(unit, '(50("-"))')
      do iat = 1, 3
         write(unit, "(2x,4x,1x,a4,1x,4x,*(es11.3))") &
            & comp(iat), sigma(:, iat)
      end do
      write(unit, '(50("-"))')
      write(unit, "(a)")
   end if

   if (qgrad) then
      write(unit, '(a,":", t25, es20.13, 1x, a)') &
      & "Charge gradient norm", norm2(dqdr), "a.u./a0"
      write(unit, '(72("-"))')
      write(unit, "(a10,1x,a4,3x,a6,1x,a4,3x,*(1x,a12))") &
         & "#", "Z", "#", "A", "dQ(Z)/dx(A)", "dQ(Z)/dy(A)", "dQ(Z)/dz(A)"
      write(unit, '(72("-"))')
      do iat = 1, mol%nat
         isp = mol%id(iat)
         do jat = 1, mol%nat
            jsp = mol%id(jat)
            write(unit, "(i10,1x,i3,1x,a2,1x,i6,1x,i3,1x,a2,*(2x ,es11.3))") &
            & iat, mol%num(isp), mol%sym(isp), jat, mol%num(jsp), mol%sym(jsp), &
            & dqdr(:, jat, iat)
         end do
      end do
      write(unit, '(72("-"))')
      write(unit, "(a)")

      write(unit, '(a,":")') &
      & "Charge virial"
      write(unit, '(62("-"))')
      write(unit, "(a10,1x,a4,3x,a9,1x,*(1x, a10))") &
         & "#", "Z", "component", "x", "y", "z"
      write(unit, '(62("-"))')
      do iat = 1, mol%nat
         isp = mol%id(iat)
         do jat = 1, 3
            write(unit, "(i10,1x,i3,1x,a2, 2x, a4, 5x,*(es11.3))") &
            & iat, mol%num(isp), mol%sym(isp), comp(jat), dqdL(:, jat, iat)
         end do
      end do
      write(unit, '(62("-"))')
      write(unit, "(a)")
   end if

end subroutine write_ascii_results


!> Write computed properties as a JSON object
subroutine json_results(unit, indentation, energy, gradient, dqdr, charges, cn)
   integer, intent(in) :: unit
   character(len=*), intent(in), optional :: indentation
   real(wp), intent(in), optional :: energy
   real(wp), intent(in), optional :: gradient(:, :)
   !> Derivative of the partial charges w.r.t. the Cartesian coordinates
   real(wp), intent(in), optional :: dqdr(:, :, :)
   real(wp), intent(in), optional :: charges(:)
   real(wp), intent(in), optional :: cn(:)
   character(len=:), allocatable :: indent, version_string
   character(len=*), parameter :: jsonkey = "('""',a,'"":',1x)"
   real(wp), allocatable :: array(:)

   call get_multicharge_version(string=version_string)

   if (present(indentation)) then
      indent = indentation
   else
      indent = ""
   end if

   write(unit, '("{")', advance="no")
   if (allocated(indent)) write(unit, "(/,a)", advance="no") repeat(indent, 1)
   write(unit, jsonkey, advance="no") "version"
   write(unit, "(1x,a)", advance="no") '"'//version_string//'"'
   if (present(energy)) then
      write(unit, '(",")', advance="no")
      if (allocated(indent)) write(unit, "(/,a)", advance="no") repeat(indent, 1)
      write(unit, jsonkey, advance="no") "energy"
      write(unit, "(1x,es25.16)", advance="no") energy
   end if
   if (present(gradient)) then
      write(unit, '(",")', advance="no")
      if (allocated(indent)) write(unit, "(/,a)", advance="no") repeat(indent, 1)
      write(unit, jsonkey, advance="no") "gradient"
      array = reshape(gradient, [size(gradient)])
      call write_json_array(unit, array, indent)
   end if
   if (present(dqdr)) then
      write(unit, '(",")', advance="no")
      if (allocated(indent)) write(unit, "(/,a)", advance="no") repeat(indent, 1)
      write(unit, jsonkey, advance="no") "dq/dr"
      array = reshape(dqdr, [size(dqdr)])
      call write_json_array(unit, array, indent)
   end if
   if (present(charges)) then
      write(unit, '(",")', advance="no")
      if (allocated(indent)) write(unit, "(/,a)", advance="no") repeat(indent, 1)
      write(unit, jsonkey, advance="no") "charges"
      array = reshape(charges, [size(charges)])
      call write_json_array(unit, array, indent)
   end if
   if (present(cn)) then
      write(unit, '(",")', advance="no")
      if (allocated(indent)) write(unit, "(/,a)", advance="no") repeat(indent, 1)
      write(unit, jsonkey, advance="no") "coordination numbers"
      array = reshape(cn, [size(cn)])
      call write_json_array(unit, array, indent)
   end if
   if (allocated(indent)) write(unit, "(/)", advance="no")
   write(unit, '("}")')

end subroutine json_results


subroutine write_json_array(unit, array, indent)
   integer, intent(in) :: unit
   real(wp), intent(in) :: array(:)
   character(len=:), allocatable, intent(in) :: indent
   integer :: i
   write(unit, '("[")', advance="no")
   do i = 1, size(array)
      if (allocated(indent)) write(unit, "(/,a)", advance="no") repeat(indent, 2)
      write(unit, "(es23.16)", advance="no") array(i)
      if (i /= size(array)) write(unit, '(",")', advance="no")
   end do
   if (allocated(indent)) write(unit, "(/,a)", advance="no") repeat(indent, 1)
   write(unit, '("]")', advance="no")
end subroutine write_json_array


!> Print header for CG solver
subroutine print_cg_header(unit, verbosity, maxit, tol, timer)

   !> Output unit
   integer, intent(in) :: unit

   !> Verbosity level
   integer, intent(in) :: verbosity

   !> Maximum number of iterations
   integer, intent(in) :: maxit

   !> Convergence tolerance
   real(wp), intent(in) :: tol

   !> Timer holding the accumulated initialization time
   type(timer_type), intent(in), optional :: timer

   if (verbosity > 1) then
      write(unit, "(a)") "Using Conjugate Gradient Solver"
      write(unit, "(a)")
      write(unit, "(a, 1x, i6)") "Max iterations : ", maxit
      write(unit, "(a, 1x, es10.2)") "Tolerance      : ", tol
      write(unit, "(a)") "Preconditioner : Jacobi (Diagonal)"
      write(unit, "(a, 1x, a)") "Initialisation time:", &
         & format_time(timer%get("initialization"))
      write(unit, "(a)") ""
      write(unit, "(2X,A,6X,A,8X,A,6X,A,4X,A)") &
         & "iter", "|residual|", "step", "relative residual", "Time / s"
   else if (verbosity == 1) then
      write(unit, "(a)") "Using Conjugate Gradient Solver"
      write(unit, "(a)")
      write(unit, "(a, 1x, i6)") "Max iterations : ", maxit
      write(unit, "(a, 1x, es10.2)") "Tolerance      : ", tol
      write(unit, "(a)") "Preconditioner : Jacobi (Diagonal)"
      write(unit, "(a)") ""
      write(unit, "(2X,A,6X,A,8X,A,6X,A)") &
         & "iter", "|residual|", "step", "relative residual"
   end if

end subroutine print_cg_header


!> Print the CG convergence message
subroutine print_cg_convergence(unit, iter, res_norm, verbosity)

   !> Output unit
   integer, intent(in) :: unit

   !> Iteration at which convergence was reached
   integer, intent(in) :: iter

   !> Converged residual norm
   real(wp), intent(in) :: res_norm

   !> Verbosity level
   integer, intent(in) :: verbosity

   if (verbosity > 0) then
      write(unit, "(a)") ""
      write(unit, "(a, i0, a, es15.5)") &
         "CG converged in ", iter, " iterations with residual norm ", res_norm
      write(unit, "(a)") ""
   end if

end subroutine print_cg_convergence


!> Print the CG iteration progress
subroutine print_cg_iteration(unit, iter, res_norm, step, rel_resnorm, &
      & verbosity, timer)

   !> Output unit
   integer, intent(in) :: unit

   !> Current iteration number
   integer, intent(in) :: iter

   !> Current residual norm
   real(wp), intent(in) :: res_norm

   !> Current step length
   real(wp), intent(in) :: step

   !> Current relative residual norm
   real(wp), intent(in) :: rel_resnorm

   !> Verbosity level
   integer, intent(in) :: verbosity

   !> Timer holding the accumulated iteration time
   type(timer_type), intent(in), optional :: timer

   if (verbosity == 1) then
      write(unit, "(i6,*(1x, es15.5))") iter, res_norm, step, rel_resnorm
   else if (verbosity > 1) then
      write(unit, "(i6,*(1x, es15.5))") iter, res_norm, step, rel_resnorm, &
         & timer%get("iteration")
   end if

end subroutine print_cg_iteration


!> Print the final timing summary of the CG solver
subroutine print_cg_final(unit, timer, verbosity)

   !> Output unit
   integer, intent(in) :: unit

   !> Timer holding the accumulated execution time
   type(timer_type), intent(in) :: timer

   !> Verbosity level
   integer, intent(in) :: verbosity

   if (verbosity > 1) then
      write(unit, "(a, 1x, a)") "CG total time : ", &
         & format_time(timer%get("total"))
      write(unit, "(a)") ""
   end if

end subroutine print_cg_final


!> Print the direct solver banner
subroutine write_direct_solver(unit)

   !> Output unit
   integer, intent(in) :: unit

   write(unit, "(a)") "Using Direct Solver"
   write(unit, "(a)")

end subroutine write_direct_solver


!> Print the final timing summary of the direct solver
subroutine print_direct_final(unit, timer, verbosity)

   !> Output unit
   integer, intent(in) :: unit

   !> Timer holding the accumulated execution time
   type(timer_type), intent(in) :: timer

   !> Verbosity level
   integer, intent(in) :: verbosity

   if (verbosity > 1) then
      write(unit, "(a, 1x, a)") "Direct solver time : ", &
         & format_time(timer%get("total"))
      write(unit, "(a)") ""
   end if

end subroutine print_direct_final


end module multicharge_output
