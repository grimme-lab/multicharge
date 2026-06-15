program tester
  use multicharge_version, only : get_multicharge_version
  implicit none
  character(len=:), allocatable :: version
  call get_multicharge_version(string=version)
  print *, version
end program tester
