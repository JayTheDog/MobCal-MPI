module get_commons
  use xtb_mctc_accuracy, only: wp
  implicit none

  !common/mobilcal/b2max(100),parab2max(100),temp1,temp2
  !common/mpicon/imyrank,inprocs,imp_per_node,inp_per_node,s_time
  !common/runparams/itn,inp,imp,igas

  private

  public  :: mpi_parameters, lj_parameters2, ff_parameters, constants, &
    & coordinates, angles, trajectory, printswitch, read_inp

  type  :: mpi_parameters
    integer :: ierr, imyrank,inprocs,imp_per_node,inp_per_node,s_time
    integer :: itn,inp,imp,igas
    integer :: i1,i2,i3,i4,i5,i6

    real(wp) :: b2max(100),parab2max(100),temp1,temp2
  end type

  type  :: lj_parameters2
    real(wp) :: RijStar(100000),eij(100000)
  end type

  type  :: ff_parameters
    real(wp) :: alphaN2, NiN2, GiN2, NiHe, GiHe
    real(wp) :: RN2Star, alphaHe, RHeStar
    real(wp) :: MMFF_B, MMFF_beta
  end type

  type :: constants
    real(wp) :: mu,ro,eo,ro2,dipol,emax,m1,m2
    real(wp) :: mconst,romax
  end type

  type :: coordinates
    real(wp) :: fx(100000),fy(100000),fz(100000)
    real(wp) :: ox(100000),oy(100000),oz(100000)
  end type

  type :: angles
    real(wp) :: theta,phi,gamma
  end type

  type :: trajectory
    real(wp) :: sw1,sw2,dtsf1,dtsf2,cmin,ifail,ifailc,inwr
  end type

  type :: printswitch
    real(wp) ::ip,it,iu1,iu2,iu3,iv,im2,im4,igs 
  end type

  type :: read_inp
    integer :: inatom, icoord
  end type

end module get_commons

