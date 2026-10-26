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

program main
   use, intrinsic :: iso_fortran_env, only: output_unit, error_unit, input_unit
   use mctc_cutoff, only: get_lattice_points
   use mctc_env, only : error_type, fatal_error, get_argument, wp, timer_type, &
      & format_time
   use mctc_io, only: structure_type, read_structure, filetype, get_filetype
   use multicharge, only : mchrg_model_type, mchrg_model, mchrg_cache, &
      & new_eeq2019_model, new_eeqbc2025_model, get_multicharge_version, &
      & write_ascii_model, write_ascii_properties, write_ascii_results
   use multicharge_output, only: json_results
   use multicharge_solver, only : new_mchrg_solver, mchrg_solver_type, &
      & mchrg_solver_input, cg_input, direct_input
   implicit none
   character(len=*), parameter :: prog_name = "multicharge"
   character(len=*), parameter :: json_output = "multicharge.json"

   character(len=:), allocatable :: input, chargeinput
   integer, allocatable :: input_format
   integer :: stat, unit, model_id
   type(error_type), allocatable :: error
   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model
   type(mchrg_cache), allocatable :: cache
   class(mchrg_solver_type), allocatable :: solver
   class(mchrg_solver_input), allocatable :: solver_input
   logical :: grad, egrad, qgrad, json, exist
   real(wp), allocatable :: trans(:, :)
   real(wp), allocatable :: energy(:), gradient(:, :), sigma(:, :)
   real(wp), allocatable :: qvec(:)
   real(wp), allocatable :: dqdr(:, :, :), dqdL(:, :, :)
   real(wp), allocatable :: charge
   integer, allocatable :: verbosity
   type(timer_type) :: timer

   call timer%push("total")

   call get_arguments(input, model_id, input_format, egrad, qgrad, charge, &
      & json, solver_input, verbosity, error)
   if (allocated(error)) then
      write(error_unit, "(a)") error%message
      error stop
   end if

   call new_mchrg_solver(solver, solver_input, error)
   if (allocated(error)) then
      write(error_unit, "(a)") error%message
      error stop
   end if

   if (input == "-") then
      if (.not. allocated(input_format)) input_format = filetype%xyz
      call read_structure(mol, input_unit, input_format, error)
   else
      call read_structure(mol, input, error, input_format)
   end if
   if (allocated(error)) then
      write(error_unit, "(a)") error%message
      error stop
   end if

   if (allocated(charge)) then
      mol%charge = charge
   else
      chargeinput = ".CHRG"
      inquire(file=chargeinput, exist=exist)
      if (exist) then
         open(file=chargeinput, newunit=unit)
         allocate(charge)
         read(unit, *, iostat=stat) charge
         if (stat == 0) then
            mol%charge = charge
            write(output_unit, "(a,/)") &
               "[Info] Molecular charge read from '"//chargeinput//"'"
         else
            write(output_unit, "(a,/)") &
               "[Warn] Could not read molecular charge from '"//chargeinput//"'"
         end if
         close(unit)
      end if
   end if

   allocate(cache)

   call timer%push("model_setup")

   if (model_id == mchrg_model%eeq2019) then
      call new_eeq2019_model(mol, model, error)
   else if (model_id == mchrg_model%eeqbc2025) then
      call new_eeqbc2025_model(mol, model, error)
   else
      call fatal_error(error, "Invalid model was chosen.")
   end if
   if (allocated(error)) then
      write(error_unit, "(a)") error%message
      error stop
   end if

   call get_lattice_points(mol%periodic, mol%lattice, model%ncoord%cutoff, trans)

   call timer%pop

   call write_ascii_model(output_unit, mol, model, verbosity, &
      & timer%get("model_setup"))

   allocate(energy(mol%nat), qvec(mol%nat))
   energy(:) = 0.0_wp
   qvec(:) = 0.0_wp

   if (egrad) then
      allocate(gradient(3, mol%nat), sigma(3, 3))
      gradient(:, :) = 0.0_wp
      sigma(:, :) = 0.0_wp
   end if

   if (qgrad) then
      allocate(dqdr(3, mol%nat, mol%nat), dqdL(3, 3, mol%nat))
      dqdr(:, :, :) = 0.0_wp
      dqdL(:, :, :) = 0.0_wp
   end if

   grad = egrad .or. qgrad

   call timer%push("update")
   call model%update(mol, cache, trans, grad)
   call timer%pop
   if (verbosity > 1) then
      write(output_unit, "(a, 1x, a)") "Get coordination number time : ", &
         & format_time(timer%get("update"))
   end if
   call model%solve(mol, solver, cache, error, &
      & energy, gradient, sigma, qvec, dqdr, dqdL, &
      & verbosity=verbosity, unit=output_unit)

   if (allocated(error)) then
      write(error_unit, "(a)") error%message
      error stop
   end if

   call write_ascii_properties(output_unit, mol, model, cache%cn, qvec)
   call write_ascii_results(output_unit, mol, energy, gradient, sigma, dqdr, dqdL)

   call timer%pop
   if (verbosity > 1) then
      write(output_unit, "(a, 1x, a)") "Total execution time : ", &
         & format_time(timer%get("total"))
   end if

   if (json) then
      open(file=json_output, newunit=unit)
      call json_results(unit, "  ", energy=sum(energy), gradient=gradient, &
         & dqdr=dqdr, charges=qvec, cn=cache%cn)
      close(unit)
      write(output_unit, "(a)") &
         "[Info] JSON dump of results written to '"//json_output//"'"
   end if

contains

subroutine help(unit)
   integer, intent(in) :: unit

   write(unit, "(a, *(1x, a))") &
      "Usage: "//prog_name//" [options] <input>"

   write(unit, "(a)") &
      "", &
      "Electronegativity equilibration model for atomic charges and", &
      "higher multipole moments", &
      ""

   write(unit, "(2x, a, t45, a)") &
      "-m, -model, --model <model>", "Choose the charge model (eeq or eeqbc)", &
      "-i, -input, --input <format>", "Hint for the format of the input file", &
      "-c, -charge, --charge <value>", "Provide the molecular charge", &
      "-solver, --solver <type>", &
      & "Provide the partial charge solver: 'cg' or 'direct' (default)", &
      "-it, -maxiter, --maxiter <int>", "Provide the maximal number of CG iterations", &
      "-tol, -tolerance, --tolerance <real>", &
      & "Provide the convergence tolerance of the CG solver", &
      "-g, -eg, -grad, --grad, -egrad, --egrad", &
      & "Evaluate molecular energy gradient and virial.", &
      "-qg, -qgrad, --qgrad", "Evaluate molecular charge gradient and virial.", &
      "-v, -verbose, --verbose", "Show more", &
      "-s, -silent, --silent", "Show less", &
      "-j, -json, --json", &
      & "Provide output in JSON format to the file 'multicharge.json'", &
      "-version, --version", "Print program version and exit", &
      "-h, -help, --help", "Show this help message"

   write(unit, "(a)")

end subroutine help

subroutine version(unit)
   integer, intent(in) :: unit
   character(len=:), allocatable :: version_string

   call get_multicharge_version(string=version_string)
   write(unit, "(a, *(1x, a))") &
      & prog_name, "version", version_string

end subroutine version

subroutine get_arguments(input, model_id, input_format, egrad, qgrad, charge, &
   & json, solver_input, verbosity, error)

   !> Input file name
   character(len=:), allocatable, intent(out) :: input

   !> ID of choosen model type
   integer, intent(out) :: model_id

   !> Input file format
   integer, allocatable, intent(out) :: input_format

   !> Evaluate energy gradient
   logical, intent(out) :: egrad

   !> Evaluate charge gradient
   logical, intent(out) :: qgrad

   !> Charge
   real(wp), allocatable, intent(out) :: charge

   !> Provide JSON output
   logical, intent(out) :: json

   !> Solver input configuration
   class(mchrg_solver_input), allocatable, intent(out) :: solver_input

   !> Verbosity level
   integer, allocatable, intent(out) :: verbosity

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: iarg, narg, iostat
   character(len=:), allocatable :: arg, solver_name
   integer, allocatable :: maxiter
   real(wp), allocatable :: tol

   model_id = mchrg_model%eeq2019
   egrad = .false.
   qgrad = .false.
   json = .false.
   iarg = 0
   verbosity = 1
   narg = command_argument_count()
   do while(iarg < narg)
      iarg = iarg + 1
      call get_argument(iarg, arg)
      select case(arg)
      case("-h", "-help", "--help")
         call help(output_unit)
         stop
      case("-version", "--version")
         call version(output_unit)
         stop
      case("-v", "-verbose", "--verbose")
         verbosity = verbosity + 1
      case("-s", "-silent", "--silent")
         verbosity = verbosity - 1
      case default
         if (.not. allocated(input)) then
            call move_alloc(arg, input)
            cycle
         end if
         call fatal_error(error, "Too many positional arguments present")
         exit
      case("-m", "-model", "--model")
         iarg = iarg + 1
         call get_argument(iarg, arg)
         if (.not. allocated(arg)) then
            call fatal_error(error, "Missing argument for model")
            exit
         end if
         if (arg == "eeq2019" .or. arg == "eeq") then
            model_id = mchrg_model%eeq2019
         else if (arg == "eeqbc2025" .or. arg == "eeqbc") then
            model_id = mchrg_model%eeqbc2025
         else
            call fatal_error(error, "Invalid model")
            exit
         end if
      case("-i", "-input", "--input")
         iarg = iarg + 1
         call get_argument(iarg, arg)
         if (.not. allocated(arg)) then
            call fatal_error(error, "Missing argument for input format")
            exit
         end if
         input_format = get_filetype("."//arg)
      case("-c", "-charge", "--charge")
         if (allocated(charge)) then
            call fatal_error(error, "Cannot use multiple charges")
            exit
         end if
         iarg = iarg + 1
         call get_argument(iarg, arg)
         if (.not. allocated(arg)) then
            call fatal_error(error, "Missing argument for charge")
            exit
         end if
         allocate(charge)
         read(arg, *, iostat=iostat) charge
         if (iostat /= 0) then
            call fatal_error(error, "Invalid charge value")
            exit
         end if
      case("-g", "-eg", "-grad", "--grad", "-egrad", "--egrad")
         egrad = .true.
      case("-qg", "-qgrad", "--qgrad")
         qgrad = .true.
      case("-j", "-json", "--json")
         json = .true.
      case("-solver", "--solver")
         if (allocated(solver_name)) then
            call fatal_error(error, "Cannot use multiple solvers")
            exit
         end if
         iarg = iarg + 1
         call get_argument(iarg, solver_name)
         if (.not. allocated(solver_name)) then
            call fatal_error(error, "Missing argument for solver")
            exit
         end if
         if (solver_name == "DIRECT" .or. solver_name == "direct") then
            allocate(direct_input :: solver_input)
         else if (solver_name == "CG" .or. solver_name == "cg") then
            allocate(cg_input :: solver_input)
         else
            call fatal_error(error, "Invalid solver '"//solver_name//"'")
            exit
         end if
      case("-it", "-maxiter", "--maxiter")
         if (.not. allocated(maxiter)) allocate(maxiter)
         iarg = iarg + 1
         call get_argument(iarg, arg)
         if (.not. allocated(arg)) then
            call fatal_error(error, "Missing argument for maximal number of iterations")
            exit
         end if
         read(arg, *, iostat=iostat) maxiter
         if (iostat /= 0 .or. maxiter <= 0) then
            call fatal_error(error, "Invalid maximal number of iterations")
            exit
         end if
      case("-tol", "-tolerance", "--tolerance")
         if (.not. allocated(tol)) allocate(tol)
         iarg = iarg + 1
         call get_argument(iarg, arg)
         if (.not. allocated(arg)) then
            call fatal_error(error, "Missing argument for tolerance")
            exit
         end if
         read(arg, *, iostat=iostat) tol
         if (iostat /= 0 .or. .not. tol > 0.0_wp) then
            call fatal_error(error, "Invalid tolerance")
            exit
         end if
      end select
   end do
   if (allocated(error)) return

   ! Default solver is direct
   if (.not. allocated(solver_input)) then
      allocate(direct_input :: solver_input)
   end if

   select type (solver_input)
   type is (cg_input)
      if (qgrad) then
         call fatal_error(error, &
            & "Charge gradient cannot be evaluated using cg solver.")
         return
      end if
      if (allocated(maxiter)) then
         solver_input%cgmiter = maxiter
      end if
      if (allocated(tol)) then
         solver_input%cgtol = tol
      end if
      if (allocated(verbosity)) then
         solver_input%verbosity = verbosity
      end if
   type is (direct_input)
      if (allocated(maxiter) .or. allocated(tol)) then
         call fatal_error(error, "Maximal number of iterations and tolerance "//&
            & "can only be used with the cg solver.")
         return
      end if
      if (allocated(verbosity)) then
         solver_input%verbosity = verbosity
      end if
   end select

   if (.not. allocated(input)) then
      if (.not. allocated(error)) then
         call help(output_unit)
         error stop
      end if
   end if

end subroutine get_arguments

end program main
