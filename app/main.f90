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
   use, intrinsic :: iso_fortran_env, only : output_unit, error_unit, input_unit
   use mctc_env, only : error_type, fatal_error, get_argument, wp
   use mctc_io, only : structure_type, read_structure, filetype, get_filetype
   use multicharge, only : mchrg_model_type, new_eeq2019_model, &
      & write_ascii_model, write_ascii_properties, write_ascii_results, &
      & get_coordination_number, get_covalent_rad, get_lattice_points, &
      & get_multicharge_version
   use multicharge_output, only : json_results
   implicit none
   character(len=*), parameter :: prog_name = "multicharge"
   character(len=*), parameter :: json_output = "multicharge.json"

   character(len=:), allocatable :: input, chargeinput
   integer, allocatable :: input_format
   integer :: stat, unit
   type(error_type), allocatable :: error
   type(structure_type) :: mol
   type(mchrg_model_type) :: model
   logical :: grad, json, exist
   real(wp), parameter :: cn_max = 8.0_wp, cutoff = 25.0_wp
   real(wp), allocatable :: cn(:), dcndr(:, :, :), dcndL(:, :, :), rcov(:), trans(:, :)
   real(wp), allocatable :: energy(:), gradient(:, :), sigma(:, :)
   real(wp), allocatable :: qvec(:), dqdr(:, :, :), dqdL(:, :, :)
   real(wp), allocatable :: charge

   call get_arguments(input, input_format, grad, charge, json, error)
   if (allocated(error)) then
      write(error_unit, '(a)') error%message
      error stop
   end if

   if (input == "-") then
      if (.not.allocated(input_format)) input_format = filetype%xyz
      call read_structure(mol, input_unit, input_format, error)
   else
      call read_structure(mol, input, error, input_format)
   end if
   if (allocated(error)) then
      write(error_unit, '(a)') error%message
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
            write(output_unit, '(a,/)') &
               "[Info] Molecular charge read from '"//chargeinput//"'"
         else
            write(output_unit, '(a,/)') &
               "[Warn] Could not read molecular charge read from '"//chargeinput//"'"
         end if
         close(unit)
      end if
   end if

   call new_eeq2019_model(mol, model)
   call get_lattice_points(mol%periodic, mol%lattice, cutoff, trans)

   call write_ascii_model(output_unit, mol, model)

   allocate(cn(mol%nat))
   if (grad) then
      allocate(dcndr(3, mol%nat, mol%nat), dcndL(3, 3, mol%nat))
   end if

   rcov = get_covalent_rad(mol%num)
   call get_coordination_number(mol, trans, cutoff, rcov, cn, dcndr, dcndL, cut=cn_max)

   allocate(energy(mol%nat), qvec(mol%nat))
   energy(:) = 0.0_wp
   if (grad) then
      allocate(gradient(3, mol%nat), sigma(3, 3))
      gradient(:, :) = 0.0_wp
      sigma(:, :) = 0.0_wp
   end if

   call model%solve(mol, cn, dcndr, dcndL, energy, gradient, sigma, &
      & qvec, dqdr, dqdL)

   call write_ascii_properties(output_unit, mol, model, cn, qvec)
   call write_ascii_results(output_unit, mol, energy, gradient, sigma)

   if (json) then
      open(file=json_output, newunit=unit)
      call json_results(unit, "  ", energy=sum(energy), gradient=gradient, charges=qvec, cn=cn)
      close(unit)
      write(output_unit, '(a)') &
         "[Info] JSON dump of results written to '"// json_output //"'"
   end if

contains


subroutine help(unit)
   integer, intent(in) :: unit

   write(unit, '(a, *(1x, a))') &
      "Usage: "//prog_name//" [options] <input>"

   write(unit, '(a)') &
      "", &
      "Electronegativity equilibration model for atomic charges and", &
      "higher multipole moments", &
      ""

   write(unit, '(2x, a, t25, a)') &
      "-i, -input, --input <format>", "Hint for the format of the input file", &
      "-c, -charge, --charge <value>", "Set the molecular charge", &
      "-g, -grad, --grad", "Evaluate molecular gradient and virial", &
      "-j, -json, --json", "Provide output in JSON format to the file 'multicharge.json'", &
      "-v, -version, --version", "Print program version and exit", &
      "-h, -help, --help", "Show this help message"

   write(unit, '(a)')

end subroutine help


subroutine version(unit)
   integer, intent(in) :: unit
   character(len=:), allocatable :: version_string

   call get_multicharge_version(string=version_string)
   write(unit, '(a, *(1x, a))') &
      & prog_name, "version", version_string

end subroutine version


subroutine get_arguments(input, input_format, grad, charge, json, error)

   !> Input file name
   character(len=:), allocatable :: input

   !> Input file format
   integer, allocatable, intent(out) :: input_format

   !> Evaluate gradient
   logical, intent(out) :: grad

   !> Provide JSON output
   logical, intent(out) :: json

   !> Charge
   real(wp), allocatable, intent(out) :: charge

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: iarg, narg, iostat
   character(len=:), allocatable :: arg

   grad = .false.
   json = .false.
   iarg = 0
   narg = command_argument_count()
   do while(iarg < narg)
      iarg = iarg + 1
      call get_argument(iarg, arg)
      select case(arg)
      case("-h", "-help", "--help")
         call help(output_unit)
         stop
      case("-v", "-version", "--version")
         call version(output_unit)
         stop
      case default
         if (.not.allocated(input)) then
            call move_alloc(arg, input)
            cycle
         end if
         call fatal_error(error, "Too many positional arguments present")
         exit
      case("-i", "-input", "--input")
         iarg = iarg + 1
         call get_argument(iarg, arg)
         if (.not.allocated(arg)) then
            call fatal_error(error, "Missing argument for input format")
            exit
         end if
         input_format = get_filetype("."//arg)
      case("-c", "-charge", "--charge")
         iarg = iarg + 1
         call get_argument(iarg, arg)
         if (.not.allocated(arg)) then
            call fatal_error(error, "Missing argument for charge")
            exit
         end if
         allocate(charge)
         read(arg, *, iostat=iostat) charge
         if (iostat /= 0) then
            call fatal_error(error, "Invalid charge value")
            exit
         end if
      case("-g", "-grad", "--grad")
         grad = .true.
      case("-j", "-json", "--json")
         json = .true.
      end select
   end do

   if (.not.allocated(input)) then
      if (.not.allocated(error)) then
         call help(output_unit)
         error stop
      end if
   end if

end subroutine get_arguments

end program main
