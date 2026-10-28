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

!> @file multicharge/model/type.f90
!> Provides a general base class for the charge models

!> General charge model
module multicharge_model_type

   use, intrinsic :: iso_fortran_env, only : output_unit
   use mctc_cutoff, only: get_lattice_points
   use mctc_env, only : timer_type, format_time, error_type, fatal_error, wp
   use mctc_io, only: structure_type
   use mctc_io_constants, only: pi
   use mctc_io_math, only: matinv_3x3
   use mctc_ncoord, only: ncoord_type
   use multicharge_blas, only: gemv, symv, gemm
   use multicharge_model_cache, only : mchrg_cache, ensure_shape
   use multicharge_model_hessian, only : add_symmetric_product
   use multicharge_solver_type, only : mchrg_solver_type
   implicit none
   private

   public :: mchrg_model_type, get_dir_trans, get_rec_trans

   !> Abstract multicharge model type
   type, abstract :: mchrg_model_type
      !> Electronegativity
      real(wp), allocatable :: chi(:)
      !> Charge width
      real(wp), allocatable :: rad(:)
      !> Chemical hardness
      real(wp), allocatable :: eta(:)
      !> CN scaling factor for electronegativity
      real(wp), allocatable :: kcnchi(:)
      !> Local charge scaling factor for electronegativity
      real(wp), allocatable :: kqchi(:)
      !> Local charge scaling factor for chemical hardness
      real(wp), allocatable :: kqeta(:)
      !> CN scaling factor for charge width
      real(wp), allocatable :: kcnrad
      !> Coordination number
      class(ncoord_type), allocatable :: ncoord
      !> Electronegativity weighted CN for local charge
      class(ncoord_type), allocatable :: ncoord_en
   contains
      !> Solve linear equations for the charge model
      procedure :: solve
      !> Contract the second charge response with fixed atomic weights
      procedure :: add_charge_hessian
      !> Model-specific curvature u^T (x_ab - A_ab q)
      procedure :: add_model_hessian
      !> Calculate local charges from electronegativity weighted CN
      procedure :: local_charge
      !> Update cache
      procedure(update), deferred :: update
      !> Calculate capacitance matrix
      procedure(get_capacitance_matrix), deferred :: get_capacitance_matrix
      !> Calculate right-hand side (electronegativity)
      procedure(get_xvec), deferred :: get_xvec
      !> Calculate xvec Gradients
      procedure(get_xvec_derivs), deferred :: get_xvec_derivs
      !> Calculate Coulomb matrix
      procedure(get_coulomb_matrix), deferred :: get_coulomb_matrix
      !> Calculate Coulomb matrix derivatives
      procedure(get_coulomb_derivs), deferred :: get_coulomb_derivs

   end type mchrg_model_type


   abstract interface

      !> Update model-dependent quantities and cache
      subroutine update(self, mol, cache, trans, grad)
         import :: mchrg_model_type, structure_type, mchrg_cache, wp

         !> Multicharge model type
         class(mchrg_model_type), intent(in) :: self

         !> Structure type
         type(structure_type), intent(in) :: mol

         !> Multicharge cache with CN, local charges, and Wigner-Seitz cell
         type(mchrg_cache), intent(inout) :: cache

         !> Lattice vectors
         real(wp), intent(in) :: trans(:, :)

         !> Flag to compute derivatives (dcndr, dcndL, dqlocdr, dqlocdL)
         logical, intent(in) :: grad
      end subroutine update

      !> Construct a capacitance matrix using cached CN and charge data
      subroutine get_capacitance_matrix(self, mol, ndim, cache)
         import :: mchrg_model_type, structure_type, mchrg_cache, wp

         !> Multicharge model type
         class(mchrg_model_type), intent(in) :: self

         !> Structure type
         type(structure_type), intent(in) :: mol

         !> System size
         integer, intent(in) :: ndim

         !> Multicharge cache holding capacitance matrix and derivatives
         type(mchrg_cache), intent(inout) :: cache
      end subroutine get_capacitance_matrix

      !> Coulomb interaction matrix (A-matrix) construction
      subroutine get_coulomb_matrix(self, mol, ndim, cache)
         import :: mchrg_model_type, structure_type, mchrg_cache, wp

         !> Multicharge model type
         class(mchrg_model_type), intent(in) :: self

         !> Structure type
         type(structure_type), intent(in) :: mol

         !> System size
         integer, intent(in) :: ndim

         !> Multicharge cache holding the Coulomb matrix
         type(mchrg_cache), intent(inout) :: cache
      end subroutine get_coulomb_matrix

      !> Coulomb matrix derivatives contracted with charges
      subroutine get_coulomb_derivs(self, mol, ndim, cache)
         import :: mchrg_model_type, structure_type, mchrg_cache, wp

         !> Multicharge model type
         class(mchrg_model_type), intent(in) :: self

         !> Structure type
         type(structure_type), intent(in) :: mol

         !> System size
         integer, intent(in) :: ndim

         !> Multicharge cache holding Coulomb-matrix derivatives
         type(mchrg_cache), intent(inout) :: cache
      end subroutine get_coulomb_derivs

      !> Electronegativity vector construction
      subroutine get_xvec(self, mol, ndim, cache)
         import :: mchrg_model_type, mchrg_cache, structure_type

         !> Multicharge model type
         class(mchrg_model_type), intent(in) :: self

         !> Structure type
         type(structure_type), intent(in) :: mol

         !> System size
         integer, intent(in) :: ndim

         !> Multicharge cache holding the electronegativity vector
         type(mchrg_cache), intent(inout) :: cache
      end subroutine get_xvec

      !> Derivatives of electronegativity vector
      subroutine get_xvec_derivs(self, mol, ndim, cache)
         import :: mchrg_model_type, structure_type, mchrg_cache, wp

         !> Multicharge model type
         class(mchrg_model_type), intent(in) :: self

         !> Structure type
         type(structure_type), intent(in) :: mol

         !> System size
         integer, intent(in) :: ndim

         !> Multicharge cache holding electronegativity-vector derivatives
         type(mchrg_cache), intent(inout) :: cache
      end subroutine get_xvec_derivs

   end interface


   !> Twice pi
   real(wp), parameter :: twopi = 2.0_wp * pi

   !> Smallest positive working-precision number
   real(wp), parameter :: eps = tiny(1.0_wp)

contains

!> Generate direct lattice translation vectors for a periodic structure
subroutine get_dir_trans(mol, trans, cutoff)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Translation vectors
   real(wp), allocatable, intent(out) :: trans(:, :)

   !> Optional lattice-vector cutoff
   real(wp), intent(in), optional :: cutoff

   integer, parameter :: rep(3) = [2, 2, 2]

   if (present(cutoff)) then
      call get_lattice_points(mol%periodic, mol%lattice, cutoff, trans)
   else
      call get_lattice_points(mol%lattice, rep, .true., trans)
   end if

end subroutine get_dir_trans


!> Generate reciprocal lattice translation vectors for a periodic structure
subroutine get_rec_trans(mol, trans, cutoff)

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Reciprocal translation vectors
   !> Shape: (3, ntrans)
   real(wp), allocatable, intent(out) :: trans(:, :)
   !> Optional reciprocal-vector cutoff
   real(wp), intent(in), optional :: cutoff

   integer, parameter :: rep(3) = [2, 2, 2]

   real(wp) :: rec_lat(3, 3)

   rec_lat = twopi * transpose(matinv_3x3(mol%lattice))
   if (present(cutoff)) then
      call get_lattice_points(mol%periodic, rec_lat, cutoff, trans)
   else
      call get_lattice_points(rec_lat, rep, .false., trans)
   end if

end subroutine get_rec_trans


!> Top-level solve routine with optional persistent cache
subroutine solve(self, mol, solver, cache, error, &
   & energy, gradient, sigma, qvec, dqdr, dqdL, verbosity, unit)

   !> Electronegativity-equilibration model
   class(mchrg_model_type), intent(in) :: self

   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> The solver instance
   class(mchrg_solver_type), intent(in) :: solver
   !> Cache handling
   type(mchrg_cache), intent(inout) :: cache

   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Optional atomic partial charges result
   real(wp), intent(out), contiguous, optional :: qvec(:)
   !> Optional electrostatic energy result
   real(wp), intent(inout), contiguous, optional :: energy(:)
   !> Optional gradient for electrostatic energy
   real(wp), intent(inout), contiguous, optional :: gradient(:, :)
   !> Optional stress tensor for electrostatic energy
   real(wp), intent(inout), contiguous, optional :: sigma(:, :)
   !> Optional derivative of the atomic partial charges w.r.t. atomic positions
   real(wp), intent(out), contiguous, optional :: dqdr(:, :, :)
   !> Optional derivative of the atomic partial charges w.r.t. lattice vectors
   real(wp), intent(out), contiguous, optional :: dqdL(:, :, :)

   !> Optional print verbosity level
   integer, intent(in), optional :: verbosity

   !> Output unit
   integer, intent(in), optional :: unit

   integer :: iat, ndim

   real(wp), allocatable :: unitvec(:)
   real(wp), allocatable :: uvec(:)
   real(wp), allocatable :: vvec(:)
   real(wp) :: uvecsum
   real(wp) :: vvecsum
   real(wp) :: lambda
   real(wp), allocatable :: daqxdr(:, :, :)
   real(wp), allocatable :: daqxdL(:, :, :)

   logical :: grad, cpq
   logical :: add_lagr
   type(timer_type) :: timer
   integer :: print_unit, verbosity_solve

   ! Calculate gradient if the respective arrays are present
   grad = present(gradient) .and. present(sigma)
   cpq = present(dqdr) .and. present(dqdL)

   if (.not. present(verbosity)) then
      verbosity_solve = 0
   else
      verbosity_solve = verbosity
   end if

   if (present(unit)) then
      print_unit = unit
   else
      print_unit = output_unit
   end if

   ! The CG solver requires a positive-definite system
   if (solver%need_pos_def) then
      if (cpq) then
         call fatal_error(error, "Charge derivatives require the inverse "//&
            & "matrix, which is not available with an iterative solver")
         return
      end if
      ndim = mol%nat
      add_lagr = .false.
   else
      ndim = mol%nat + 1
      add_lagr = .true.
   end if

   call timer%push("total")
   call timer%push("setup")

   ! Setup the system matrices and vectors
   call self%get_capacitance_matrix(mol, ndim, cache)
   call self%get_coulomb_matrix(mol, ndim, cache)
   call self%get_xvec(mol, ndim, cache)
   call ensure_shape(cache%vrhs, mol%nat + 1)

   ! pop setup timer
   call timer%pop

   ! Print header
   call print_solve_header(print_unit, verbosity_solve, timer%get("setup"))

   if (add_lagr) then
      if (cpq) then
         ! The inverse matrix is needed for the charge derivatives
         call ensure_shape(cache%ainv, ndim, ndim)
         call solver%solve(amat=cache%amat, xvec=cache%xvec, &
            & vrhs=cache%vrhs, ainv=cache%ainv, new_unit=print_unit, &
            & error=error)
      else
         if (allocated(cache%ainv)) deallocate(cache%ainv)
         call solver%solve(amat=cache%amat, xvec=cache%xvec, &
            & vrhs=cache%vrhs, new_unit=print_unit, error=error)
      end if
      if (allocated(error)) return

   else
      allocate(unitvec(mol%nat))
      allocate(uvec(mol%nat))
      allocate(vvec(mol%nat))
      ! Initial guess
      do iat = 1, mol%nat
         uvec(iat) = 1.0_wp / (cache%amat(iat, iat) + eps)
         vvec(iat) = - cache%xvec(iat) / (cache%amat(iat, iat) + eps)
      end do

      unitvec = 1.0_wp

      call print_constrained_system_message(print_unit, verbosity_solve, "u")
      ! Constrained response: A*uvec = 1
      call solver%solve(amat=cache%amat, xvec=unitvec, &
         & vrhs=uvec, new_unit=print_unit, error=error)
      if (allocated(error)) return
      call print_constrained_system_message(print_unit, verbosity_solve, "v")

      ! Unconstrained response: A*vvec = -xvec
      call solver%solve(amat=cache%amat, xvec=-cache%xvec, &
         & vrhs=vvec, new_unit=print_unit, error=error)
      if (allocated(error)) return
      uvecsum = sum(uvec)
      vvecsum = sum(vvec)
      ! Lagrangian multiplier
      lambda = - (mol%charge + vvecsum) / (uvecsum + eps)

      ! Projection of uvec on vvec
      cache%vrhs(:mol%nat) = -vvec - lambda * uvec
      cache%vrhs(mol%nat + 1) = lambda

   end if

   ! Partial charges if present
   if (present(qvec)) then
      qvec(:) = cache%vrhs(:mol%nat)
   end if

   ! Electrostatic energy if present
   if (present(energy)) then
      call timer%push("energy")
      call symv(cache%amat, cache%vrhs(:ndim), cache%xvec(:ndim), &
         & alpha=0.5_wp, beta=-1.0_wp, uplo="l")
      if (ndim > mol%nat) then
         ! Correct xvec to exclude constraint term
         cache%xvec(:mol%nat) = cache%xvec(:mol%nat) &
            & - 0.5_wp * cache%vrhs(mol%nat + 1)
      end if
      energy(:) = energy(:) + cache%vrhs(:mol%nat) * cache%xvec(:mol%nat)
      call timer%pop
      call print_energy_time(print_unit, verbosity_solve, timer%get("energy"))
   end if

   ! Allocate and get amat and xvec derivatives
   if (grad .or. cpq) then
      call timer%push("setup_gradient")
      call timer%push("dxdr_setup")
      call self%get_xvec_derivs(mol, ndim, cache)
      call timer%pop
      if (verbosity_solve > 1) then
         write(print_unit, "(a, 1x, a)") &
            & "Electronegativity derivatives setup time : ", &
            & format_time(timer%get("dxdr_setup"))
         write(print_unit, "(a)") ""
      end if
      call timer%push("dadr_setup")
      call self%get_coulomb_derivs(mol, ndim, cache)
      call timer%pop
      if (verbosity_solve > 1) then
         write(print_unit, "(a, 1x, a)") &
            & "Coulomb matrix derivatives setup time : ", &
            & format_time(timer%get("dadr_setup"))
         write(print_unit, "(a)") ""
      end if

      ! pop gradient setup
      call timer%pop
      call print_gradient_header(print_unit, verbosity_solve, &
         & timer%get("setup_gradient"))
   end if

   ! Calculate gradients if requested
   if (grad) then
      call timer%push("gradient")
      call gemv(cache%dadr(:, :, :mol%nat), cache%vrhs(:mol%nat), gradient, &
         & beta=1.0_wp, alpha=0.5_wp)
      call gemv(cache%dxdr(:, :, :mol%nat), cache%vrhs(:mol%nat), gradient, &
         & beta=1.0_wp, alpha=-1.0_wp)
      call gemv(cache%dadL(:, :, :mol%nat), cache%vrhs(:mol%nat), sigma, &
         & beta=1.0_wp, alpha=0.5_wp)
      call gemv(cache%dxdL(:, :, :mol%nat), cache%vrhs(:mol%nat), sigma, &
         & beta=1.0_wp, alpha=-1.0_wp)
      ! pop gradient timer
      call timer%pop
      call print_gradient_time(print_unit, verbosity_solve, &
         & timer%get("gradient"))
   end if

   ! Calculate charge derivatives if requested
   if (cpq) then
      allocate(daqxdr(3, mol%nat, ndim), source=0.0_wp)
      allocate(daqxdL(3, 3, ndim), source=0.0_wp)

      call timer%push("cpq")
      do iat = 1, mol%nat
         daqxdr(:, :, iat) = cache%dxdr(:, :, iat) - cache%dadr(:, :, iat)
         daqxdL(:, :, iat) = cache%dxdL(:, :, iat) - cache%dadL(:, :, iat)
      end do
      call gemm(daqxdr, cache%ainv(:, :mol%nat), dqdr, alpha=1.0_wp)
      call gemm(daqxdL, cache%ainv(:, :mol%nat), dqdL, alpha=1.0_wp)
      ! pop cpq timer
      call timer%pop
      call print_gradient_time(print_unit, verbosity_solve, timer%get("cpq"))
   end if

   ! pop total solve timer
   call timer%pop
   call print_total_time(print_unit, verbosity_solve, timer%get("total"))

end subroutine solve

!> Add the Cartesian charge Hessian contracted with fixed weights and lattice.
!> Ay=x, A^T u=(weights,0): H_ab=u^T(x_ab-A_ab y-A_a y_b-A_b y_a).
subroutine add_charge_hessian(self, mol, error, weights, hessian, cache)
   class(mchrg_model_type), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(error_type), allocatable, intent(out) :: error
   real(wp), intent(in) :: weights(:)
   real(wp), intent(inout) :: hessian(:, :)

   !> Response cache retained by get_charges
   type(mchrg_cache), intent(inout) :: cache

   integer :: nat
   real(wp), allocatable :: adj(:), rhs(:), qvec(:)
   real(wp), allocatable :: dadr(:, :, :), dadL(:, :, :)

   nat = mol%nat
   if (any(mol%periodic) .and. .not. all(mol%periodic)) then
      call fatal_error(error, "Analytical charge Hessian requires molecular or 3D periodic boundary conditions")
      return
   end if
   if (size(weights) /= nat .or. any(shape(hessian) /= 3*nat)) then
      call fatal_error(error, "Invalid dimensions for charge Hessian contraction")
      return
   end if
   if (.not. allocated(cache%vrhs) .or. .not. allocated(cache%ainv) &
      & .or. .not. allocated(cache%dqdr) .or. .not. allocated(cache%cn)) then
      call fatal_error(error, "Charge Hessian requires a response cache from get_charges")
      return
   end if
   if (size(cache%vrhs) /= nat + 1 .or. any(shape(cache%ainv) /= nat + 1) &
      & .or. any(shape(cache%dqdr) /= [3, nat, nat])) then
      call fatal_error(error, "Charge response cache has incompatible dimensions")
      return
   end if

   allocate(adj(nat+1), rhs(nat+1))
   rhs(:nat) = weights
   rhs(nat+1) = 0.0_wp
   call symv(cache%ainv, rhs, adj, uplo="l")
   call self%add_model_hessian(mol, error, cache, adj, hessian)
   if (allocated(error)) return

   ! Reuse the first derivative action with the adjoint in place of the charges
   call move_alloc(cache%vrhs, qvec)
   call move_alloc(cache%dadr, dadr)
   call move_alloc(cache%dadL, dadL)
   cache%vrhs = adj
   call self%get_coulomb_derivs(mol, nat + 1, cache)
   call add_symmetric_product(reshape(cache%dadr(:, :, :nat), [3 * nat, nat]), &
      & reshape(cache%dqdr, [3 * nat, nat]), hessian, alpha=-1.0_wp)
   call move_alloc(qvec, cache%vrhs)
   call move_alloc(dadr, cache%dadr)
   call move_alloc(dadL, cache%dadL)
end subroutine add_charge_hessian


!> Concrete models override this curvature contraction when supported.
subroutine add_model_hessian(self, mol, error, cache, adj, hessian)
   class(mchrg_model_type), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(error_type), allocatable, intent(out) :: error

   !> Response cache retained by get_charges
   type(mchrg_cache), intent(inout) :: cache
   real(wp), intent(in) :: adj(:)
   real(wp), intent(inout) :: hessian(:, :)

   call fatal_error(error, "Analytical charge Hessian is not implemented for this model")
end subroutine add_model_hessian


subroutine local_charge(self, mol, trans, qloc, dqlocdr, dqlocdL)
   !> Electronegativity equilibration model
   class(mchrg_model_type), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Lattice points
   real(wp), intent(in) :: trans(:, :)
   !> Local atomic partial charges
   real(wp), intent(out) :: qloc(:)
   !> Optional derivative of local atomic partial charges w.r.t. atomic positions
   real(wp), intent(out), optional :: dqlocdr(3, mol%nat, mol%nat)
   !> Optional derivative of local atomic partial charges w.r.t. lattice vectors
   real(wp), intent(out), optional :: dqlocdL(3, 3, mol%nat)

   qloc = 0.0_wp
   if (present(dqlocdr) .and. present(dqlocdL)) then
      dqlocdr = 0.0_wp
      dqlocdL = 0.0_wp
   end if
   ! Get the electronegativity weighted CN for local charge
   ! Derivatives depend only in this CN
   if (allocated(self%ncoord_en)) then
      call self%ncoord_en%get_coordination_number(mol, trans, qloc, dqlocdr, dqlocdL)
   end if

   ! Distribute the total charge equally
   qloc = qloc + mol%charge / real(mol%nat, wp)

end subroutine local_charge


!> Print header for charge equilibration solver
subroutine print_solve_header(unit, verbosity, timer)

   !> Output unit
   integer, intent(in) :: unit

   !> Verbosity level
   integer, intent(in) :: verbosity

   !> Elapsed setup time
   real(wp), intent(in) :: timer

   if (verbosity > 0) then
      write(unit, '(54("-"))')
      write(unit, "(13x, a)") "Charge equilibration solver"
      write(unit, '(54("-"))')
      write(unit, "(a)") ""
      if (verbosity > 1) then
         write(unit, "(a, 1x, a)") "Setup time : ", format_time(timer)
         write(unit, "(a)") ""
      end if
   end if
end subroutine print_solve_header


!> Print header for gradient calculations
subroutine print_gradient_header(unit, verbosity, timer)

   !> Output unit
   integer, intent(in) :: unit

   !> Verbosity level
   integer, intent(in) :: verbosity

   !> Elapsed gradient setup time
   real(wp), intent(in) :: timer

   if (verbosity > 0) then
      write(unit, '(54("-"))')
      write(unit, "(17x, a)") "Gradient Calculations"
      write(unit, '(54("-"))')
      write(unit, "(a)") ""
      if (verbosity > 1) then
         write(unit, "(a, 1x, a)") "Gradient setup time : ", format_time(timer)
         write(unit, "(a)") ""
      end if
   end if
end subroutine print_gradient_header


!> Print message for constrained system solves
subroutine print_constrained_system_message(unit, verbosity, vector)

   !> Output unit
   integer, intent(in) :: unit

   !> Verbosity level
   integer, intent(in) :: verbosity

   !> Constrained-system identifier
   character, intent(in) :: vector

   if (verbosity > 0) then
      if (vector == "u") then
         write(unit, "(a)") "Solving constrained system: J*u = 1"
      else if (vector == "v") then
         write(unit, "(a)") "Solving unconstrained system: J*v = chi"
      end if
      write(unit, "(a)") ""
   end if
end subroutine print_constrained_system_message


!> Print gradient calculation time
subroutine print_gradient_time(unit, verbosity, timer)

   !> Output unit
   integer, intent(in) :: unit

   !> Verbosity level
   integer, intent(in) :: verbosity

   !> Elapsed gradient calculation time
   real(wp), intent(in) :: timer

   if (verbosity > 1) then
      write(unit, "(a, 1x, a)") &
         & "Gradient calculation time : ", format_time(timer)
      write(unit, "(a)") ""
   end if
end subroutine print_gradient_time


!> Print energy calculation time
subroutine print_energy_time(unit, verbosity, timer)

   !> Output unit
   integer, intent(in) :: unit

   !> Verbosity level
   integer, intent(in) :: verbosity

   !> Elapsed energy calculation time
   real(wp), intent(in) :: timer

   if (verbosity > 1) then
      write(unit, "(a, 1x, a)") "Energy calculation time : ", format_time(timer)
      write(unit, "(a)") ""
   end if
end subroutine print_energy_time


!> Print total solve time
subroutine print_total_time(unit, verbosity, timer)

   !> Output unit
   integer, intent(in) :: unit

   !> Verbosity level
   integer, intent(in) :: verbosity

   !> Elapsed total solve time
   real(wp), intent(in) :: timer

   if (verbosity > 1) then
      write(unit, "(a, 1x, a)") "Total solve time : ", format_time(timer)
      write(unit, "(a)") ""
   end if
end subroutine print_total_time

end module multicharge_model_type
