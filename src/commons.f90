module get_commons
  use xtb_mctc_accuracy, only: wp
  implicit none

  !common/mobilcal/b2max(100),parab2max(100),temp1,temp2
  !common/mpicon/imyrank,inprocs,imp_per_node,inp_per_node,s_time
  !common/runparams/itn,inp,imp,igas

  private

  public  :: mpi_parameters

  type  :: mpi_parameters


    integer :: ierr, imyrank,inprocs,imp_per_node,inp_per_node,s_time
    integer :: itn,inp,imp,igas
    integer :: i1,i2,i3,i4,i5,i6

    real(wp) :: b2max(100),parab2max(100),temp1,temp2


  end type


end module get_commons

