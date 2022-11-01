!     PROGRAM MOBCAL_MPI v1.2
!
!     Mobility Calculation program published in ...
!**************************************************************************
!**************************************************************************
!     Adapted from original code with reference:
!     PROGRAM MOBCAL
!
!     Program to Calculate Mobilities
!
!     See: M. F. Mesleh, J. M. Hunter, A. A. Shvartsburg, G. C. Schatz,
!     and M. F. Jarrold, "Structural Information from Ion Mobility
!     Measurements: Effects of the Long Range Potential" J. Phys. Chem.
!     1996, 100, 16082-16086; A. A. Shvartsburg and M. F. Jarrold,
!     "An Exact Hard Spheres Scattering Model for the Mobilities of
!     Polyatomic Ions" Chem. Phys. Lett. 1996, 261, 86.
!
!     Version of 06/15/00 junko.f
!     Modified to include Na+ parameters and how masses are read
!     Version of 09/14/98
!     Corrected center of mass calculation in fcoord and ncoord and
!     changed how structural asymmetry parameter is calculated
!     Version of 08/15/98
!     Modified to calculate structural asymmetry parameter
!     Version of 02/20/98
!     Modified GSANG: removed bug that caused trajectory start position
!     to go into endless loop for long structures. Also fixed TESTXRAND
!     and summary print out.
!     Version of 01/29/98
!     Modified MOBIL2: removed a bug that messed up b2max calculation
!     for large molecules
!     Version of 12/24/97
!     Modified to permit averaging over multiple conformations
!     Version of 12/08/97
!     Extensively modified to allow for a molecule consisting of
!     an unlimited variety of different atoms. Exact hard sphere
!     scattering incorporated. RANLUX random number generator
!     incorporated. Several obsolete subroutines removed.
!     Version of 10/16/97
!     Modified to allow uniform and non-uniform charge distributions
!     Version of 08/23/97
!     Modified to allow calculations with a non-uniform charge distribution
!     Version of 21/10/19
!
!     Version 17/01/21
!     Adding multiple temp runs in same input file
!***  Explicit seed and correction to calculation of ion mobility (K); does not affect previous CCS calculations
!
!     MOBIL2 calculates the mobility using a Lennard-Jones plus
!     ion-induced dipole potential. MOBIL2 uses Monte Carlo integrations
!     over orientation and impact parameter and a numerical integration
!     over g* (the reduced velocity). MOBIL2 includes second order
!     corrections, though these don't appear to be important.
!
!     GSANG/DERIV/DIFFEQ are subroutines derived from code provided by
!     George Schatz. They use Runge-Kutta-Gill and Adams-Moulton predictor
!     -corrector integration methods. These subroutines are now set-up to
!     work in three dimensions.
!
!     DLJPOT calculates the potential and the derivatives of the potential.
!     The potential is given by a sum of 6-12 two body interactions and
!     ion-induced dipole interactions. The charge can be distributed over
!     all atoms equally, a specific charge distribution can be employed, or
!     the charge can be set to zero (only the 6-12 part of the potential
!     is used). Each atom in the cluster or molecule can have different
!     Lennard-Jones parameters.
!
!     MOBIL4 calculates the exact hard sphere scattering mobility and
!     the projection approximation mobility. Adapted from code written
!     by Alexandre Shvartsburg.
!
!     XRAND/RANLUX are random number generators. RANLUX is the standard
!     ranlux number generator of M. Luscher (code from F. James). XRAND
!     is a function that calls either RANLUX or RAND (the standard F77
!     random number generator).
!
!     FCOORD reads in the coordinates, integer masses, and partial charges.
!
!     RANTATE/ROTATE rotates the coordinates of the cluster.
!
!     TRAJ calculates a series of trajectories for closer examination.
!
!     TRAJONE calculates one trajectory for a given velocity, impact
!     parameter, and angles (ang%theta, ang%phi, and ang%gamma in degrees).
!
!     POTENT determines the average potential. (Only for near spherical
!     molecules or clusters).
!
!
!     ***************************************************************
!******************************************************************************
!******************************************************************************
!     Major changes include inclusion of "softer" exp-6 Lennard-Jones
!     potential and more exhaustive paramater set which needs to be
!     included as input
!
!

! Version updated by Jay -- Jeroen Koopman

program MobCal
  use constants
  use mpi
  use get_commons, only: mpi_parameters, lj_parameters2, ff_parameters, constants, &
    & coordinates, trajectory, read_inp
  use xtb_mctc_accuracy
  implicit none
   
  !include 'mpif.h'

  integer :: i, i1
  integer :: iadd
  real(wp) :: t, mob, cs, sdevpc
  
  real(wp) :: asymp
  real(wp) :: xrand
  
  integer, parameter :: ixlen=100000

  integer :: itcnt
  integer :: ip, it, iu1, iu2, iu3, iv, im2, im4, igs
  integer :: itloop

  real(wp) :: AiHe, AiN2
  real(wp) :: pcharge(100000)
  
  character(len=50) :: filen1,filen2,unit,dchar,xlabel,infile
  character(len=2) :: flend(32)
  

  real(wp) dens

  !****NEED COMMON BLOCK FOR MPI ELEMENTS
  !common/mpicon/imyrank,inprocs,imp_per_node,inp_per_node,s_time
  !common/runparams/itn,inp,imp,igas
  !****Add ability for temparray.. up to max of 50 temps
  !dimension temparray(50)
  !dimension tmmarray(50),csarray(50),sdevarray(50)
  real(wp) :: temparray(50)
  real(wp) ::  tmmarray(50),csarray(50),sdevarray(50)

  !******
  data flend/'_1','_2','_3','_4','_5','_6','_7','_8',&
  &'_9','10','11','12','13','14','15','16','17','18',&
  &'19','20','21','22','23','24','25','26','27','28',&
  &'29','30','31','32'/


  type(mpi_parameters) :: mpiP
  type(lj_parameters2) :: lj
  type(ff_parameters) :: mmff 
  type(constants) :: const
  type(coordinates) :: coord
  type(trajectory) :: trj
  type(read_inp) :: inp


  !****Start MPI Stuff
  call MPI_INIT(mpiP%ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,mpiP%inprocs,mpiP%ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,mpiP%imyrank,mpiP%ierr)
  mpiP%s_time=MPI_WTIME()
  !
  !***Read input file and output file name
  if(mpiP%imyrank.eq.0)then
     call getarg(1,infile)
     open(unit=20,file=infile)
     read(20,'(a50)') filen1
     read(20,'(a50)') filen2
     mpiP%i2=-85790349
     close(20)
  endif
  !
  !************
  !     Define parameters for MOBIL2
  !
  !     Number of complete cycles for average mobility calculation in
  !     MOBIL2. Default value is 10.
  !
  mpiP%itn=10
  !
  !     Number of points in velocity integration in MOBIL2. Default
  !     value is 20.
  !
  mpiP%inp=40
  !
  !     Number of points in Monte Carlo integrations of impact parameter
  !     and orientation in MOBIL2. Default value is 500.
  !
  mpiP%imp=500
  mpiP%igas=2
  !****Only node 0 opens output and input file
  if(mpiP%imyrank.eq.0)then
     open(unit=8,file=filen2)
  else
     open(unit=1000+mpiP%imyrank,file='output'//flend(mpiP%imyrank)//'.out',&
     &status='unknown')
  endif
  
  !******
  !
  !     print switches ip=1  print scattering angles
  !                    it=1  print trajectory
  !                    iu1=1 print initial coordinates from FCOORD
  !                    iu2=1 print angles and coordinates from ROTATE
  !                    iu3=1 print angles from ROTATE
  !                    iv=1  print all potentials from POTENT
  !                    im2=1 print brief information from MOBIL2
  !                    im4=1 print brief information from MOBIL4
  !                    igs=1 print out angles, v, and b from MOBIL4 into
  !                          a temporary file called hold
  !
  ip=0
  it=0
  iu1=0
  iu2=0
  iu3=0
  iv=0
  im2=0
  im4=0
  igs=0
  !
  !
  !     Temperature
  !      t=301.d0
  !      t=298.d0
  !      xmv=8.20573660809596d-5*t
  
  !****
  !
  !*** Define MMFF parameters for Diatomic Nitrogen
  mmff%alphaN2=1.739688487d0
  mmff%NiN2=5.918d0
  AiN2=3.1781544d0
  mmff%GiN2=1.16175013d0
  mmff%RN2Star=AiN2*(mmff%alphaN2**(0.25d0))
  !*** Define MMFF parameters for Helium
  mmff%alphaHe=0.205236d0
  mmff%NiHe=1.42d0
  AiHe=4.40d0
  mmff%GiHe=1.209d0
  mmff%RHeStar=AiHe*(mmff%alphaHe**(0.25d0))
  !****************
  mmff%MMFF_B=0.2d0
  mmff%MMFF_beta=12.0d0
  !*****
  if(mpiP%imyrank.eq.0)then
     call fcoord(filen1,unit,dchar,xlabel,asymp,temparray,itcnt)
     close(9)
  endif
  !*******************
  !****Broadcast******
  call MPI_BCAST(itcnt,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpiP%ierr)
  call MPI_BCAST(temparray,itcnt,MPI_DOUBLE_PRECISION,0,&
  &MPI_COMM_WORLD,mpiP%ierr)
  call MPI_BCAST(mpiP%itn,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpiP%ierr)
  call MPI_BCAST(mpiP%i2,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpiP%ierr)
  call MPI_BCAST(mpiP%inp,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpiP%ierr)
  call MPI_BCAST(mpiP%imp,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpiP%ierr)
  call MPI_BCAST(mpiP%igas,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpiP%ierr)
  call MPI_BCAST(inp%inatom,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpiP%ierr)
  call MPI_BCAST(inp%icoord,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpiP%ierr)
  call MPI_BCAST(const%m2,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiP%ierr)
  call MPI_BCAST(const%romax,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiP%ierr)
  call MPI_BCAST(coord%fx,ixlen,MPI_DOUBLE_PRECISION,0,&
  &MPI_COMM_WORLD,mpiP%ierr)
  call MPI_BCAST(coord%fy,ixlen,MPI_DOUBLE_PRECISION,0,&
  &MPI_COMM_WORLD,mpiP%ierr)
  call MPI_BCAST(coord%fz,ixlen,MPI_DOUBLE_PRECISION,0,&
  &MPI_COMM_WORLD,mpiP%ierr)
  call MPI_BCAST(pcharge,ixlen,MPI_DOUBLE_PRECISION,0,&
  &MPI_COMM_WORLD,mpiP%ierr)
  call MPI_BCAST(lj%eij,ixlen,MPI_DOUBLE_PRECISION,0,&
  &MPI_COMM_WORLD,mpiP%ierr)
  call MPI_BCAST(lj%RijStar,ixlen,MPI_DOUBLE_PRECISION,0,&
  &MPI_COMM_WORLD,mpiP%ierr)
  call MPI_BCAST(coord%ox,ixlen,MPI_DOUBLE_PRECISION,0,&
  &MPI_COMM_WORLD,mpiP%ierr)
  call MPI_BCAST(coord%oy,ixlen,MPI_DOUBLE_PRECISION,0,&
  &MPI_COMM_WORLD,mpiP%ierr)
  call MPI_BCAST(coord%oz,ixlen,MPI_DOUBLE_PRECISION,0,&
  &MPI_COMM_WORLD,mpiP%ierr)
  !******k*************
  !***Added loop here over temperature array
  do itloop=1,itcnt
  !*******************
  !     Temperature
  !
  !      t=301.d0
  !      t=298.d0
     t=temparray(itloop)
  !     Lennard-Jones scaling parameters
  !
     const%eo=1.34d-03*xe
     const%ro=3.043d0*1.0d-10
     const%ro2=const%ro*const%ro
     if(mpiP%imyrank.eq.0)write(8,600) const%eo/xe,const%ro/1.0d-10
  !
  !     Constant for ion-induced dipole potential
  !
  !     Polarizability of helium = 0.204956d-30 m3
  !     xeo is permitivity of vacuum, 8.854187817d-12 F.m-1
  !
     if(mpiP%igas.eq.2)then
        const%dipol=1.740d-30/(2.d0*4.d0*pi*xeo)
        const%dipol=const%dipol*xe*xe
        if(mpiP%imyrank.eq.0)write(8,601) const%dipol
  !
  !     Mass constants
  !
        const%m1=28.0134d0
        const%mu=((const%m1*const%m2)/(const%m1+const%m2))/(xn*1.0d3)
     elseif(mpiP%igas.eq.1)then
        const%dipol=0.204956d-30/(2.0d0*4.0d0*pi*xeo)
        const%dipol=const%dipol*xe*xe
        if(mpiP%imyrank.eq.0)write(8,601) const%dipol
  !
  !     Mass constants
  !
        const%m1=4.0026d0
        const%mu=((const%m1*const%m2)/(const%m1+const%m2))/(xn*1.0d3)
     endif
  !
  !     Mobility constant
  !
     const%mconst=dsqrt(18.d0*pi)/16.d0
     const%mconst=const%mconst*dsqrt(xn*1.0d3)*dsqrt((1.d0/const%m1)+(1.d0/const%m2))
     const%mconst=const%mconst*xe/dsqrt(xk)
     dens=xn/(8.20573660809596d-5*t)
     const%mconst=const%mconst/dens
  !      if(mpiP%imyrank.eq.0)write(8,602) const%mconst
  !
  
     if(mpiP%imyrank.eq.0)write(8,603) t
  !
  !     Define parameters for random number generator
  !
  !     If mpiP%i5=1 RANLUX is used otherwise RAND is used. If RAND is used
  !     mpiP%i2 is the seed integer. If RANLUX is used mpiP%i2, mpiP%i3, and i4 are seed
  !     integers. mpiP%i3 and i4 (which are used to start RANLUX in a particular
  !     place) are usually set to zero. i1 contain the "luxury level" - how
  !     good the random number generator is. Values of 0 to 4 can be used
  !     and the default value in RANLUX is 3. At this level "any
  !     theoretically possible correlations have a very small chance of
  !     being observed".
  !
     i1=3
     mpiP%i3=0
     mpiP%i4=0
     mpiP%i5=1
     if(mpiP%i5.ne.1) then
        if(mpiP%imyrank.eq.0)write(8,604) mpiP%i2
     else
        if(mpiP%imyrank.eq.0)write(8,619)
     endif
  !     initialize the random number generators
  !****Broadcast mpiP%i2 and initialize several times to mix things up...
     call MPI_BCAST(mpiP%i2,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpiP%ierr)
     do i=1,mpiP%imyrank+1
        call rluxgo(i1,mpiP%i2,mpiP%i3,mpiP%i4)
        call srand(mpiP%i2)
        mpiP%i2=int(xrand()*1d8+xrand()*1d7+&
        &xrand()*1d6+xrand()*1d5+&
        &xrand()*1d4+xrand()*1d3+&
        &xrand()*1d2)
     enddo
  !************
  !***Calc mpiP%inp per node
     if(mod(mpiP%inp,mpiP%inprocs).ne.0)then
        iadd=mpiP%inprocs-mod(mpiP%inp,mpiP%inprocs)
        mpiP%inp=mpiP%inp+iadd
     endif
     mpiP%inp_per_node=mpiP%inp/mpiP%inprocs
  !*********************
  !***Test mpiP%imp to make sure it fits evenly over mpiP%inprocs
     if(mod(mpiP%imp,mpiP%inprocs).ne.0)then
        iadd=mpiP%inprocs-mod(mpiP%imp,mpiP%inprocs)
        mpiP%imp=mpiP%imp+iadd
     endif
     mpiP%imp_per_node=mpiP%imp/mpiP%inprocs
  
  !      mpiP%imp=500
  !
  !     Minimum value of (1-cosX). This quantity determines the maximum
  !     impact parameter at each velocity. Default value is 0.0005.
  !
     trj%cmin=0.0005d0
  !
  !     Define some parameters for trajectories: trj%sw1 defines the potential
  !     energy where the trajectory starts and trj%dtsf1 is related to the
  !     time step at the start of the trajectory. As the trajectory comes
  !     close to a collision the time step is reduced. trj%sw2 defines the
  !     potential energy where the time step is reduced and trj%dtsf2 defines
  !     the reduced time step. Default values are:
  !     trj%sw1 = 0.00005   trj%dtsf1=0.5
  !     trj%sw2 = 0.0025    trj%dtsf2=0.05
  !     trj%inwr is the number of integration steps before the program tests
  !     to see if the trajectory is done or lost. trj%ifail is the number of
  !     failed trajectories that are permitted (a failed trajectory is one
  !     that does not conserve energy to within 1%. Default values are:
  !     trj%inwr = 1        trj%ifail = 100
  !
     trj%sw1=0.00005d0
     trj%sw2=0.005d0
     trj%dtsf1=0.5d0
     trj%dtsf2=0.1d0
     trj%inwr=1
     trj%ifail=100
     trj%ifailc=0
  !******
     mob=0.d0
     cs=0.0d0
  !***********************************************************
  !***********************************************************
     call mobil2(t,mob,cs,sdevpc)
  !***********************************************************
  !***********************************************************
     tmmarray(itloop)=mob
     csarray(itloop)=cs
     sdevarray(itloop)=sdevpc
  !     print out summary
  !********************************
  !***Only print out for 0 node
  !********************************
     if(mpiP%imyrank.eq.0)then
  !
        write(8,605) filen1,xlabel
        write(8,618) 1.d0/mob,cs*1.d20,sdevpc,trj%ifailc
        write(8,615)mpiP%itn,mpiP%inp,mpiP%imp,mpiP%itn*mpiP%inp*mpiP%imp
        if(mpiP%igas.eq.2)then
           write(8,*)'Mobility Calculated under N2 gas'
        elseif(mpiP%igas.eq.1)then
           write(8,*)'Mobility Calculated under He gas'
        endif
     endif
  !****End temperature loop here
  enddo
  !****Collective output:
  write(8,'(t1,a)')'*********Mobility Summary*********'
  write(8,'(t1,a)')' temp  ,average CCS ,percent error '
  do itloop=1,itcnt
     if(mpiP%imyrank.eq.0)then
        write(8,'(t1,f7.3,",",e13.6,",",e13.6)')&
        &temparray(itloop),csarray(itloop)*1.d20,sdevarray(itloop)
     endif
  enddo
  !********************************
  !***Only print out for 0 node
  !********************************
  if(mpiP%imyrank.eq.0)then
     close(8)
  else
     close(1000+mpiP%imyrank)
  endif
  !********************************
  !********************************
  call MPI_FINALIZE(mpiP%ierr)
  !****
  601 format(1x,'dipole constant =',1pe11.4)
  602 format(1x,'mobility constant =',1pe11.4)
  603 format(1x,'temperature =',1pe11.4)
  619 format(/1x,'using RANLUX',/)
  604 format(1x,'using RAND with seed integer =',i8)
  635 format(/)
  621 format(/1x,'coordinate set =',i5)
  637 format(/1x,'structural asymmetry parameter =',f8.4)
  608 format(1x,'using no charge - only LJ interactions')
  612 format(/1x,'inverse average EHS mobility =',1pe11.4,&
  &/1x,'average EHS cross section =',e11.4)
  611 format(/1x,'inverse average PA mobility =',1pe11.4,&
  &/1x,'average PA cross section =',e11.4)
  609 format(/1x,'mobility calculation by MOBIL4 (HS scattering)')
  !610 format(/1x,'number of Monte Carlo trajectories =',i7,&
  !&/1x,'maximum number of reflections encountered =',mpiP%i3)
  616 format(1x,'temperature =',1pe11.4)
  617 format(1x,'using RAND with seed integer =',i8)
  620 format(1x,'using RANLUX with seed integer =',i8)
  636 format(1x,'structural asymmetry parameter =',1pe11.4)
  607 format(1x,'using a calculated (non-uniform) charge distribution')
  606 format(1x,'using a uniform charge distribution')
  605 format(///1x,'SUMMARY',//1x,'program version = junkn.f',&
  &/1x,'input file name = ',a50,/1x,'input file label = ',a50)
  615 format(/1x,'number of complete cycles (mpiP%itn) =',i6,/1x,&
  &'number of velocity points (mpiP%inp) =',i6,/1x,&
  &'number of random points (mpiP%imp) =',i6,/1x,&
  &'total number of points =',i7)
  614 format(/1x,'trajectory parameters',/1x,'trj%sw1 =',1pe11.4,7x,&
  &'trj%sw2 =',e11.4,/1x,'trj%dtsf1 =',e11.4,5x,'trj%dtsf2 =',e11.4)
  613 format(/1x,'mobility calculation by MOBIL2 (trajectory method)')
  618 format(/1x,'inverse average (second order) TM mobility =',&
  &1pe11.4,/1x,'average TM cross section =',e11.4,/1x,&
  &'standard deviation (percent) =',e11.4,/1x,&
  &'number of failed trajectories =',i4)
  634 format(/1x,'minimum and maximum number of reflections =',i3,2x,i3)
  622 format(//4x,'set',5x,'PA CS',7x,'PA MOB^-1',6x,'EHSS CS',&
  &6x,'EHSS MOB^-1',4x,'ASYMP')
  626 format(/1x,'number of Monte Carlo trajectories =',i7)
  625 format(/1x,'mobility calculation by MOBIL4 (HS scattering)')
  623 format(1x,i5,3x,1pe11.4,3x,e11.4,3x,e11.4,3x,e11.4,3x,f8.4)
  624 format(/3x,'AVGE',2x,1pe11.4,3x,e11.4,3x,e11.4,3x,e11.4,3x,f8.4)
  628 format(/1x,'trajectory parameters',/1x,'trj%sw1 =',1pe11.4,7x,&
  &'trj%sw2 =',e11.4,/1x,'trj%dtsf1 =',e11.4,5x,'trj%dtsf2 =',e11.4)
  631 format(1x,i5,3x,1pe11.4,3x,e11.4)
  629 format(/1x,'number of complete cycles (mpiP%itn) =',i6,/1x,&
  &'number of velocity points (mpiP%inp) =',i6,/1x,&
  &'number of random points (mpiP%imp) =',i6,/1x,&
  &'total number of points =',i7)
  630 format(//4x,'set',5x,'TM CS',7x,'TM MOB^-1')
  633 format(/1x,'total number of failed trajectories =',i4)
  632 format(/3x,'AVGE',2x,1pe11.4,3x,e11.4)
  627 format(/1x,'mobility calculation by MOBIL2 (trajectory method)')
  600 format(1x,'Lennard-Jones scaling parameters: const%eo=',&
  &1pe11.4,1x,'const%ro=',e11.4)
  750 format(1x,a7,1x,f5.2)
  stop
  !end
  end program MobCal
  !
  !     ***************************************************************
  !
  subroutine fcoord(filen1,unit,dchar,xlabel,asymp,temparray,itcnt)
  !
  !     Reads in coordinates and other parameters.
      use constants
      use mpi
      use get_commons, only: mpi_parameters, lj_parameters2, ff_parameters, constants,&
        & coordinates, angles, read_inp
      use xtb_mctc_accuracy
      implicit none
  !
     integer :: i, iatom, igamma, iphi, ilen, isw, k
     integer :: itcnt, istart, ifin
     integer :: correct

     real(wp) :: alphai(100000), Ni(100000),Ai(100000),Gi(100000)
     real(wp) :: RiiStar, gammaij, Rsum
     real(wp) :: asymp
     real(wp) :: temparray(50)
     real(wp) :: enersca, distsca
     real(wp) :: pcharge(100000)
     real(wp) :: xmass(100000)
     real(wp) :: acharge, tcharge
     real(wp) :: coeff1, coeff2, coeff3
     real(wp) :: fxo, fyo, fzo
     real(wp) :: hold
     real(wp) :: xyz, yz, xyzsum, yzsum

     character*50 filen1,unit,dchar,xlabel
     integer, parameter :: ixlen=100000

     character*1000 readline

     type(mpi_parameters) :: mpiP
     type(lj_parameters2) :: lj
     type(ff_parameters) :: mmff
     type(constants) :: const
     type(coordinates) :: coord
     type(angles) :: ang
     type(read_inp) :: inp
  !
     write(8,603) filen1
     open (9,file=filen1)
     read(9,'(a30)') xlabel !mol name 
     write(8,601) xlabel
     read(9,*) inp%icoord ! ? = 1
     write(8,650) inp%icoord
  !
     read(9,*) inp%inatom !no of atoms
     write(8,612) inp%inatom
     read(9,'(a30)') unit !units (angström)
     read(9,'(a30)') dchar !task
  !
     if(unit.eq.'au') write(8,611)
     if(unit.eq.'ang') write(8,614)
     if(unit.ne.'au'.and.unit.ne.'ang') then
        write(8,610)
        close(8)
        call MPI_FINALIZE(mpiP%ierr)
        stop
     endif
  !
  read(9,*) correct !dont know, should be set to = 1
     write(8,613) correct
  !
     read(9,'(t1,a1000)')readline !read entire options line in .mfj
     ilen=len(trim(readline))
     isw=0
     k=0
     itcnt=0
     do i=1,ilen
        if(isw.eq.1)then
           if(readline(i:i).eq.' '.or.i.eq.ilen)then ! go through options
              ifin=i-1
              if(i.eq.ilen)ifin=i
              isw=0
              k=k+1
              if(k.eq.1)then
                 read(readline(istart:ifin),*)mpiP%itn !no of 
              elseif(k.eq.2)then
                 read(readline(istart:ifin),*)mpiP%inp !no of 
              elseif(k.eq.3)then
                 read(readline(istart:ifin),*)mpiP%imp !no of impact iterations
              elseif(k.eq.4)then
                 read(readline(istart:ifin),*)mpiP%igas !collision gas
              elseif(k.eq.5)then
                 read(readline(istart:ifin),*)mpiP%i2 !random seed
              elseif(k.ge.6)then
                 itcnt=itcnt+1 
                 read(readline(istart:ifin),*)temparray(itcnt) ! iterate over temperatures
              endif
           endif
        endif
        if(isw.eq.0)then
           if(readline(i:i).ne.' ')then
              istart=i
              isw=1
           endif
        endif
     enddo
  !***PATCH PATCH
  !***PATCH PATCH
  !***PATCH PATCH
  !***PATCH PATCH
  !****node 0 can read in coord. then broadcast
     if(mpiP%igas.eq.2)then
        enersca=0.80d0
        distsca=0.78d0
     else
        enersca=0.81d0
        distsca=0.98d0
     endif
  !
     if(dchar.eq.'equal') write(8,630)
     if(dchar.eq.'calc') write(8,631)
     if(dchar.eq.'none') write(8,633)
     if(dchar.ne.'equal'.and.dchar.ne.'calc'.and.dchar&
     &.ne.'none') then
        write (8,632)
        close(8)
        call MPI_FINALIZE(mpiP%ierr)
        stop
     endif
  !
     tcharge=0.d0
     acharge=0.d0

  !> Read in xyz, mass, charge and vdw parameters from input
  do iatom=1,inp%inatom
        read(9,*) coord%fx(iatom),coord%fy(iatom),coord%fz(iatom),&
        & xmass(iatom), pcharge(iatom),&
        & alphai(iatom), Ni(iatom), Ai(iatom), Gi(iatom)

  !*****Generate pair paramters (lj%RijStar and lj%eij) following combination
  !*****rules
        RiiStar = Ai(iatom)*(alphai(iatom)**0.25d0)
        if(mpiP%igas.eq.2)then
           Rsum=RiiStar+ mmff%RN2Star
           gammaij=(RiiStar - mmff%RN2Star)/Rsum
           coeff1=-mmff%MMFF_beta*gammaij*gammaij
           lj%RijStar(iatom)=0.5d0*Rsum*(1.0d0+mmff%MMFF_B*(1.0d0-dexp(coeff1)))
  
           coeff2=181.16d0 * Gi(iatom) * mmff%GiN2 * alphai(iatom) * mmff%alphaN2
           coeff3=dsqrt(alphai(iatom)/Ni(iatom)) + dsqrt(mmff%alphaN2/mmff%NiN2)
        elseif(mpiP%igas.eq.1)then
           Rsum=RiiStar+mmff%RHeStar
           gammaij=(RiiStar-mmff%RHeStar)/Rsum
           coeff1=-mmff%MMFF_beta*gammaij*gammaij
           lj%RijStar(iatom)=0.5d0*Rsum*(1.0d0+mmff%MMFF_B*(1.0d0-dexp(coeff1)))
  
           coeff2=181.16d0*Gi(iatom)*mmff%GiHe*alphai(iatom)*mmff%alphaHe
           coeff3=dsqrt(alphai(iatom)/Ni(iatom))+dsqrt(mmff%alphaHe/mmff%NiHe)
        endif
        lj%eij(iatom)=coeff2/(coeff3*(lj%RijStar(iatom)**6.0d0))
  !***Unit conversion from ang to meters and kcal/mol to J
        lj%RijStar(iatom)=lj%RijStar(iatom)*1.0d-10*distsca
        lj%eij(iatom)=lj%eij(iatom)*(4.184d3/xn)*enersca
  !****well depth in MM3 is 1.1195eps
        lj%eij(iatom)=lj%eij(iatom)/1.1195d0
  
        tcharge=tcharge+pcharge(iatom)
        acharge=acharge+dabs(pcharge(iatom))
        if(unit.eq.'au') then
           coord%fx(iatom)=coord%fx(iatom)*0.52917706d0
           coord%fy(iatom)=coord%fy(iatom)*0.52917706d0
           coord%fz(iatom)=coord%fz(iatom)*0.52917706d0
        endif
  end do
  !
     if(dchar.eq.'equal') then
        do 2011 iatom=1,inp%inatom
  ! 2011 pcharge(iatom)=1.d0/dble(inp%inatom)
           pcharge(iatom)=tcharge/dble(inp%inatom)
  2011  continue
     endif
  !
     if(dchar.eq.'none') then
        do 2012 iatom=1,inp%inatom
           pcharge(iatom)=0.d0
  2012  continue
     endif
  !
  !      if(dchar.ne.'none') write(8,615) tcharge,acharge
  !**********
     const%m2=0.d0
     do 2021 iatom=1,inp%inatom
        const%m2=const%m2+xmass(iatom)
  2021 continue
     write(8,604) const%m2
  !
  !
  !      if(iu1.eq.1) write(8,620)
  !
     fxo=0.d0
     fyo=0.d0
     fzo=0.d0
     do 2009 iatom=1,inp%inatom
        fxo=fxo+(coord%fx(iatom)*xmass(iatom))
        fyo=fyo+(coord%fy(iatom)*xmass(iatom))
        fzo=fzo+(coord%fz(iatom)*xmass(iatom))
  2009 continue
     fxo=fxo/const%m2
     fyo=fyo/const%m2
     fzo=fzo/const%m2
     do 2010 iatom=1,inp%inatom
        coord%fx(iatom)=(coord%fx(iatom)-fxo)*1.d-10*correct
        coord%fy(iatom)=(coord%fy(iatom)-fyo)*1.d-10*correct
        coord%fz(iatom)=(coord%fz(iatom)-fzo)*1.d-10*correct
  2010 continue
  !
     if(inp%icoord.eq.1) close(9)
  !
     do 3000 iatom=1,inp%inatom
        coord%ox(iatom)=coord%fx(iatom)
        coord%oy(iatom)=coord%fy(iatom)
        coord%oz(iatom)=coord%fz(iatom)
  3000 continue
  !
     const%romax=0.d0
     do 3001 iatom=1,inp%inatom
        if(lj%RijStar(iatom).gt.const%romax) const%romax=lj%RijStar(iatom)
  3001 continue
  !
     const%romax=const%romax+(1.1055d-10/2.0)
  !
  !     determine structural asymmetry parameter
  !
     ang%theta=0.d0
     asymp=0.d0
     do igamma=0,360,2
        do iphi=0,180,2
           ang%gamma=dble(igamma)/cang
           ang%phi=dble(iphi)/cang
           call rotate
           xyzsum=0.d0
           yzsum=0.d0
           do iatom=1,inp%inatom
              xyz=dsqrt(coord%fx(iatom)**2+coord%fy(iatom)**2+coord%fz(iatom)**2)
              yz=dsqrt(coord%fy(iatom)**2+coord%fz(iatom)**2)
              xyzsum=xyzsum+xyz
              yzsum=yzsum+yz
           end do
           hold=((pi/4.d0)*xyzsum)/yzsum
           if(hold.gt.asymp) asymp=hold
         end do
      end do
  !
  603 format(1x,'input file name = ',a50)
  601 format(1x,'input file label = ',a50)
  650 format(1x,'number of coordinate sets =',i5)
  612 format(1x,'number of atoms =',i4)
  611 format(1x,'coordinates in atomic units')
  614 format(1x,'coordinates in angstroms')
  610 format(1x,'units not specified')
  613 format(1x,'correction factor for coordinates =',1pe11.4)
  630 format(1x,'using a uniform charge distribution')
  631 format(1x,'using a calculated (non-uniform) charge distribution')
  633 format(1x,'using no charge - only LJ interactions')
  632 format(1x,'charge distribution not specified')
  615 format(1x,'total charge =',1pe11.4,/1x,&
     &'total absolute charge =',e11.4)
  620 format(/9x,'initial coordinates',9x,'mass',3x,'charge',&
     &9x,'LJ parameters',/)
  602 format(1x,'type not defined for atom number',i3)
  604 format(1x,'mass of ion =',1pd11.4)
  623 format(1x,'center of mass coordinates = ',1pe11.4,',',e11.4,&
     &',',e11.4)
  600 format(1x,1pe11.4,1x,e11.4,1x,e11.4,1x,i3,1x,e11.4,1x,e11.4,1x,e11.4)
  621 format(/)
     return
  end subroutine fcoord
  !
  !     ***************************************************************
  !
  subroutine rantate
  !
  !     Rotates the cluster/molecule to a random orientation.
      use mpi
      use get_commons, only: constants, angles
      implicit none
  !
     !implicit double precision (a-h,m-z)
     !include 'mpif.h'
  !
    type(constants) :: const
    type(angles) :: ang

     rnt=xrand()
     rnp=xrand()
     rng=xrand()
     ang%theta=rnt*2.d0*pi
     ang%phi=dasin((rnp*2.d0)-1.d0)+(pi/2.d0)
     ang%gamma=rng*2.d0*pi
     call rotate
  !
     return
  end subroutine rantate

  !
  !     ***************************************************************
  !
  subroutine rotate
  !
  !     Rotates the cluster/molecule.
      use constants
      use mpi
      use get_commons, only: constants, coordinates, angles, read_inp
      use xtb_mctc_accuracy
      implicit none
  !
     !include 'mpif.h'

     integer :: iatom

     real(wp) :: ogamma, ngamma, ophi, nphi, otheta, ntheta
     real(wp) :: rxy, rzy

    type(constants) :: const
    type(coordinates) :: coord
    type(angles) :: ang
    type(read_inp) :: inp
  !
  !      if(mpiP%imyrank.eq.0)then
  !       if(iu2.eq.1.or.iu3.eq.1) write(8,610) ang%theta*cang,ang%phi*cang,
  !     1  ang%gamma*cang
  !      else
  !       if(iu2.eq.1.or.iu3.eq.1) write(1000+mpiP%imyrank,610)ang%theta*cang,
  !     1  ang%phi*cang,
  !     1  ang%gamma*cang
  !      endif
  !
     do iatom=1,inp%inatom
        rxy=dsqrt(coord%ox(iatom)*coord%ox(iatom)+(coord%oy(iatom)*coord%oy(iatom)))
        if(rxy.gt.0.0d0)then
           otheta=dacos(coord%ox(iatom)/rxy)
           if(coord%oy(iatom).lt.0.d0) otheta=(2.d0*pi)-otheta
           ntheta=otheta+ang%theta
        endif
        coord%fx(iatom)=dcos(ntheta)*rxy
        coord%fy(iatom)=dsin(ntheta)*rxy
    end do
  !
     do iatom=1,inp%inatom
        rzy=dsqrt(coord%oz(iatom)*coord%oz(iatom)+(coord%fy(iatom)*coord%fy(iatom)))
        if(rzy.gt.0.d0)then
           ophi=dacos(coord%oz(iatom)/rzy)
           if(coord%fy(iatom).lt.0.d0) ophi=(2.d0*pi)-ophi
           nphi=ophi+ang%phi
        endif
        coord%fz(iatom)=dcos(nphi)*rzy
        coord%fy(iatom)=dsin(nphi)*rzy
      end do
  !
     do iatom=1,inp%inatom
        rxy=dsqrt(coord%fx(iatom)*coord%fx(iatom)+(coord%fy(iatom)*coord%fy(iatom)))
        if(rxy.gt.0.d0)then
           ogamma=dacos(coord%fx(iatom)/rxy)
           if(coord%fy(iatom).lt.0.d0) ogamma=(2.d0*pi)-ogamma
           ngamma=ogamma+ang%gamma
        endif
        coord%fx(iatom)=dcos(ngamma)*rxy
        coord%fy(iatom)=dsin(ngamma)*rxy
      end do
  !
  !      if(iu2.ne.0)then
  !       if(mpiP%imyrank.eq.0)then
  !        write(8,620)
  !       else
  !        write(1000+mpiP%imyrank,620)
  !       endif
  !       do iatom=1,inp%inatom
  !        if(mpiP%imyrank.eq.0)then
  !         write(8,600) coord%ox(iatom),coord%oy(iatom),coord%oz(iatom),coord%fx(iatom),
  !     1    coord%fy(iatom),coord%fz(iatom)
  !        else
  !         write(1000+mpiP%imyrank,600) coord%ox(iatom),coord%oy(iatom),
  !     1    coord%oz(iatom),coord%fx(iatom),
  !     1    coord%fy(iatom),coord%fz(iatom)
  !        endif
  !       enddo
  !      endif
  !
  600 format(1x,1pe11.4,2(1x,e11.4),5x,3(1x,e11.4))
  610 format(//1x,'coordinates rotated by ROTATE',//1x,&
     &'ang%theta=',1pe11.4,1x,'ang%phi=',e11.4,1x,'ang%gamma=',1pe11.4,/)
  620 format(9x,'initial coordinates',24x,'new coordinates',/)
     return
  end
  !
  !     *****************************************************
  !
  subroutine dljpotHe(x,y,z,pot,dpotx,dpoty,dpotz,dmax,dchar)
  !
  !     Subroutine to calculate L-J + ion-dipole potential.
      use get_commons, only: lj_parameters2, constants, coordinates, read_inp
      use xtb_mctc_accuracy
      implicit none
  !
     integer, parameter :: ixlen=100000

     real(wp) :: x, y, z
     real(wp) :: pot, dpotx, dpoty, dpotz, dmax

     character(len=4) :: dchar

     dimension pottry(2,3)
     dimension pot_mol(3)
     dimension dpotxtry(2,3)
     dimension dpotx_mol(3)
     dimension dpotytry(2,3)
     dimension dpoty_mol(3)
     dimension dpotztry(2,3)
     dimension dpotz_mol(3)

     type(lj_parameters2) :: lj
     type(constants) :: const
     type(coordinates) :: coord
     type(read_inp) :: inp
  !
  !     nitrogen : five charge model (Allen Tildesley, 14 page)
  !     Every data from B3LYP//aug-cc-pVDZ
     dmax=2.d0*const%romax
     rx=0.d0
     ry=0.d0
     rz=0.d0
     e00=0.d0
     de00x=0.d0
     de00y=0.d0
     de00z=0.d0
     sum1=0.d0
     sum2=0.d0
     sum3=0.d0
     sum4=0.d0
     sum5=0.d0
     sum6=0.d0
     do 1100 iatom=1,inp%inatom
  !
        xx=x-coord%fx(iatom)
        xx2=xx*xx
  !
        yy=y-coord%fy(iatom)
        yy2=yy*yy
  !
        zz=z-coord%fz(iatom)
        zz2=zz*zz
  !
        rxyz2=xx2+yy2+zz2
  !
        rxyz=dsqrt(rxyz2)
  !*****
  !***** Exp6 potential calculation
        R=rxyz/lj%RijStar(iatom)
        R2=R*R
        R6=R2*R2*R2
        preexp=1.84d5
        expterm=preexp*dexp(-12.0d0*R)
        predisp=2.25d0
        dispterm=predisp/R6
  !*****
  
        if(rxyz.lt.dmax) dmax=rxyz
        rxyz5=rxyz3*rxyz2
        !**** Exp6 potential
        rxyz3=rxyz2*rxyz
        e00=e00+lj%eij(iatom)*(expterm-dispterm)
  !**** Exp6 derivative
        de00=lj%eij(iatom)*(expterm*(-12.0d0)+dispterm*6.0d0/R)&
        &/rxyz/lj%RijStar(iatom)
        de00x=de00x+(de00*xx)
        de00y=de00y+(de00*yy)
        de00z=de00z+(de00*zz)
  !     ion-induced dipole potential
        if(dchar.ne.'none')then
           if(pcharge(iatom).ne.0.d0)then
              rxyz3i=pcharge(iatom)/rxyz3
              rxyz5i=-3.d0*pcharge(iatom)/rxyz5
              rx=rx+(xx*rxyz3i)
              ry=ry+(yy*rxyz3i)
              rz=rz+(zz*rxyz3i)
  !     ion-induced dipole derivative
              sum1=sum1+(rxyz3i+(xx2*rxyz5i))
              sum2=sum2+(xx*yy*rxyz5i)
              sum3=sum3+(xx*zz*rxyz5i)
              sum4=sum4+(rxyz3i+(yy2*rxyz5i))
              sum5=sum5+(yy*zz*rxyz5i)
              sum6=sum6+(rxyz3i+(zz2*rxyz5i))
           endif
        endif
  !
  1100 continue
  !
     pot=e00-(const%dipol*((rx*rx)+(ry*ry)+(rz*rz)))
     dpotx=de00x-(const%dipol*((2.0d0*rx*sum1)+(2.0d0*ry*sum2)&
     &+(2.0d0*rz*sum3)))
     dpoty=de00y-(const%dipol*((2.0d0*rx*sum2)+(2.0d0*ry*sum4)&
     &+(2.0d0*rz*sum5)))
     dpotz=de00z-(const%dipol*((2.0d0*rx*sum3)+(2.0d0*ry*sum5)&
     &+(2.0d0*rz*sum6)))
  !
     return
  end
  !
  !     *******************************************************
  !     *****************************************************
  !
  subroutine dljpotN2(x,y,z,pot,dpotx,dpoty,dpotz,dmax,dchar)
  !
  !     Subroutine to calculate L-J + ion-dipole potential.
  !
     use mpi
     use xtb_mctc_accuracy
     use get_commons, only: mpi_parameters, lj_parameters2, constants, coordinates, &
       & read_inp
     use constants
     implicit none
     !include 'mpif.h'

     integer :: isamp, ibatom, iatom

     real(wp) :: x, y, z
     real(wp) :: pot, dpotx, dpoty, dpotz, dmax
     real(wp) :: Ptfn, xkT, pc, pc_center
     real(wp) :: dpolx, dpoly, dpolz, dipolzz, dipolxx, pot_min
     real(wp) :: qpol, dqpolx, dqpoly, dqpolz

     real(wp) :: rx, ry, rz
     real(wp) :: de00,e00, de00x, de00y, de00z, sum1, sum2, sum3, sum4, sum5, sum6
     real(wp) :: xc, yc, zc
     real(wp) :: xx_center, xx, xx_center2, xx2
     real(wp) :: yy_center, yy, yy_center2, yy2
     real(wp) :: zz_center, zz, zz_center2, zz2
     real(wp) :: rxyz_center, rxyz, rxyz_center2, rxyz2
     real(wp) :: R, R2, R6, preexp, expterm, predisp, dispterm
     real(wp) :: const_k

     real(wp) :: pottry(2,3)
     real(wp) :: pot_mol(3)
     real(wp) :: dpotxtry(2,3)
     real(wp) :: dpotx_mol(3)
     real(wp) :: dpotytry(2,3)
     real(wp) :: dpoty_mol(3)
     real(wp) :: dpotztry(2,3)
     real(wp) :: dpotz_mol(3)

     real(wp), parameter :: bond=1.0976d-10

     character(len=4) :: dchar

     type(mpi_parameters) :: mpiP
     type(lj_parameters2) :: lj
     type(constants) :: const
     type(coordinates) :: coord
     type(read_inp) :: inp

  !
  !     nitrogen : five charge model (Allen Tildesley, 14 page)
  !     Every data from B3LYP//aug-cc-pVDZ
     dmax=2.d0*const%romax
     
     Ptfn=0.d0
     xkT=500.d0*xk
     pc=-0.4825d0
     predisp=2.25d0
     preexp=1.84d5
  !     -1.447 D A
     pc_center=-(pc)
  !     B3LYP/aug-cc-pVDZ
  !      dipolzz=2.263d-30/(2.d0*4.d0*pi*xeo)*1.0
     dipolzz=1.710d-30/(2.d0*4.d0*pi*xeo)
  !      dipolzz=1.7543d-30/(2.d0*4.d0*pi*xeo)
     dipolzz=dipolzz*xe*xe
     dipolxx=1.710d-30/(2.d0*4.d0*pi*xeo)
  !      dipolxx=1.500d-30/(2.d0*4.d0*pi*xeo)*1.0
  !      dipolxx=1.7543d-30/(2.d0*4.d0*pi*xeo)
     dipolxx=dipolxx*xe*xe
     pot_min=1.0d8
  !
     do 1300 isamp=1,3
        do 1200 ibatom=1,2
           rx=0.d0
           ry=0.d0
           rz=0.d0
           e00=0.d0
           de00x=0.d0
           de00y=0.d0
           de00z=0.d0
           sum1=0.d0
           sum2=0.d0
           sum3=0.d0
           sum4=0.d0
           sum5=0.d0
           sum6=0.d0
           qpol=0.d0
           dqpolx=0.d0
           dqpoly=0.d0
           dqpolz=0.d0
           xc=0.0d0
           yc=0.0d0
           zc=0.0d0
           if(isamp.eq.1)then
              xc=(bond/2.d0)*(2.d0*dble(ibatom)-3.d0)
              dpolx=dipolzz
              dpoly=dipolxx
              dpolz=dipolxx
           elseif(isamp.eq.2)then
              yc=(bond/2.d0)*(2.d0*dble(ibatom)-3.d0)
              dpolx=dipolxx
              dpoly=dipolzz
              dpolz=dipolxx
           elseif(isamp.eq.3)then
              zc=(bond/2.d0)*(2.d0*dble(ibatom)-3.d0)
              dpolx=dipolxx
              dpoly=dipolxx
              dpolz=dipolzz
           endif
  !
           do 1100 iatom=1,inp%inatom
  !
              xx_center=x-coord%fx(iatom)
              xx=xx_center+xc
              xx_center2=xx_center*xx_center
              xx2=xx*xx
  !
              yy_center=y-coord%fy(iatom)
              yy=yy_center+yc
              yy_center2=yy_center*yy_center
              yy2=yy*yy
  !
              zz_center=z-coord%fz(iatom)
              zz=zz_center+zc
              zz_center2=zz_center*zz_center
              zz2=zz*zz
              !
              rxyz_center2=xx_center2+yy_center2+zz_center2
              rxyz2=xx2+yy2+zz2
  !
              rxyz_center=dsqrt(rxyz_center2)
              rxyz=dsqrt(rxyz2)
  !*****
  !***** Exp6 potential calculation
              R=rxyz/lj%RijStar(iatom)
              R2=R*R
              R6=R2*R2*R2
              expterm=preexp*dexp(-12.0d0*R)
              dispterm=predisp/R6
  !*****
  
              if(rxyz.lt.dmax) dmax=rxyz
              rxyz3=rxyz2*rxyz
              rxyz5=rxyz3*rxyz2
              rxyz6=rxyz5*rxyz
              rxyz8=rxyz5*rxyz3
              rxyz12=rxyz6*rxyz6
              rxyz14=rxyz12*rxyz2
              rxyz_center3=rxyz_center2*rxyz_center
              rxyz_center5=rxyz_center3*rxyz_center2
  !**** Exp6 potential
              e00=e00+lj%eij(iatom)*(expterm-dispterm)
  !**** Exp6 derivative
              de00=lj%eij(iatom)*(expterm*(-12.0d0)+dispterm*6.0d0/R)&
              &/rxyz/lj%RijStar(iatom)
              de00x=de00x+(de00*xx)
              de00y=de00y+(de00*yy)
              de00z=de00z+(de00*zz)
  !     ion-induced dipole potential
              if(dchar.ne.'none')then
                 if(pcharge(iatom).ne.0.d0)then
                    rxyz3i=pcharge(iatom)/rxyz_center3
                    rxyz5i=-3.d0*pcharge(iatom)/rxyz_center5
                    rx=rx+(xx_center*rxyz3i)
                    ry=ry+(yy_center*rxyz3i)
                    rz=rz+(zz_center*rxyz3i)
  !     ion-induced dipole derivative
                    sum1=sum1+(rxyz3i+(xx_center2*rxyz5i))
                    sum2=sum2+(xx_center*yy_center*rxyz5i)
                    sum3=sum3+(xx_center*zz_center*rxyz5i)
                    sum4=sum4+(rxyz3i+(yy_center2*rxyz5i))
                    sum5=sum5+(yy_center*zz_center*rxyz5i)
                    sum6=sum6+(rxyz3i+(zz_center2*rxyz5i))
  !     ion-partial charge coulomb potential(quadrupole)
                    const_k=pcharge(iatom)*(xe*xe)/(4.d0*pi*xeo)
                    qpol=qpol+(pc_center*const_k/rxyz_center)
                    qpol=qpol+(pc*const_k/rxyz)
  !     ion-partial charge coulomb derivative(quadrupole)
                    dqpolx=dqpolx-((pc_center*const_k/rxyz_center3)*(xx_center))
                    dqpoly=dqpoly-((pc_center*const_k/rxyz_center3)*(yy_center))
                    dqpolz=dqpolz-((pc_center*const_k/rxyz_center3)*(zz_center))
                    dqpolx=dqpolx-((pc*const_k/rxyz3)*(xx))
                    dqpoly=dqpoly-((pc*const_k/rxyz3)*(yy))
                    dqpolz=dqpolz-((pc*const_k/rxyz3)*(zz))
                 endif
              endif
  !
  1100     continue
  !
           pottry(ibatom,isamp)=e00-0.5d0*(((dpolx*rx*rx)+(dpoly*ry*ry)&
           &+(dpolz*rz*rz)))+qpol
           dpotxtry(ibatom,isamp)=de00x-0.5d0*((dpolx*2.0d0*rx*sum1)&
           &+(dpoly*2.0d0*ry*sum2)+(dpolz*2.d0*rz*sum3))+dqpolx
           dpotytry(ibatom,isamp)=de00y-0.5d0*((dpolx*2.0d0*rx*sum2)&
           &+(dpoly*2.0d0*ry*sum4)+(dpolz*2.d0*rz*sum5))+dqpoly
           dpotztry(ibatom,isamp)=de00z-0.5d0*((dpolx*2.0d0*rx*sum3)&
           &+(dpoly*2.0d0*ry*sum5)+(dpolz*2.0d0*rz*sum6))+dqpolz
  !
  !
  1200  continue
        pot_mol(isamp)=pottry(1,isamp)+pottry(2,isamp)
        tpot=pot_mol(isamp)
        if(pot_min.ge.pot_mol(isamp)) pot_min=pot_mol(isamp)
        dpotx_mol(isamp)=dpotxtry(1,isamp)+dpotxtry(2,isamp)
        dpoty_mol(isamp)=dpotytry(1,isamp)+dpotytry(2,isamp)
        dpotz_mol(isamp)=dpotztry(1,isamp)+dpotztry(2,isamp)
  1300 continue
  !
     do 1410 isamp=1,3
        temp_pot=pot_mol(isamp)-pot_min
        Ptfn=Ptfn+dexp(-temp_pot/xkT)
  1410 continue
  !
     pot=0.d0
     dpotx=0.d0
     dpoty=0.d0
     dpotz=0.d0
     do 1400 isamp=1,3
        temp_pot=pot_mol(isamp)-pot_min
        weight=dexp(-temp_pot/xkT)/Ptfn
        pot=pot+weight*pot_mol(isamp)
        dpotx=dpotx+weight*dpotx_mol(isamp)
        dpoty=dpoty+weight*dpoty_mol(isamp)
        dpotz=dpotz+weight*dpotz_mol(isamp)
  1400 continue
  !
     return
  end
  !
  !     *******************************************************
  !
  subroutine gsang(v,b,erat,ang,d1,istep)
  !
  !     Calculates trajectory. Adapted from code written by George Schatz
  !     A fixed step integrator is used which combines a runge-kutta-gill
  !     initiator with an adams-moulton predictor-corrector propagator.
  !
     use mpi
     use get_commons, only: mpi_parameters, constants, trajectory
     use xtb_mctc_accuracy
     implicit none
     !include 'mpif.h'

     integer :: ns,nw,l
     integer :: istep

     real(wp) :: x, b, v
     real(wp) :: vy, vx, vxyz
     real(wp) :: erat, ang, d1

     dimension w(6),dw(6)

  type(mpi_parameters) :: mpiP
  type(constants) :: const
  type(trajectory) :: trj

     vy=-v
     vx=0.d0
     vz=0.d0
     vxyz=dabs(vy)
  !
  !     determine time step
  !
     top=(v/95.2381d0)-0.5d0
     if(v.ge.1000.d0) top=10.d0
     if(v.ge.2000.d0) top=10.0d0-((v-2000.d0)*7.5d-3)
     if(v.ge.3000.d0) top=2.5d0
  !
     dt1=top*trj%dtsf1*1.0d-11/v
     dt2=dt1*trj%dtsf2
     dt=dt1
  !
  !     determine trajectory start position
  !
     e0=0.5d0*const%mu*v*v
     x=b
     z=0.d0
  !
     ymin=0.d0
     ymax=0.d0
  !      do 200 i=1,inp%inatom
  !       if(coord%fy(i).gt.ymax) ymax=coord%fy(i)
  !       if(coord%fy(i).lt.ymin) ymin=coord%fy(i)
  !200   continue
     ymax=ymax/1.0d-10
     ymin=ymin/1.0d-10
     iymin=nint(ymin)-1
     iymax=nint(ymax)+1
  !***********************************
  !***********************************
     id2=iymax
     istop=0
     iswitch=0
     ireturn=0
     icheck=0
     pot=0.0d0
     y=dble(id2)*1.0d-10
     if(mpiP%igas.eq.2)then
        call dljpotN2(x,y,z,pot,dpotx,dpoty,dpotz,dmax,'calc')
     elseif(mpiP%igas.eq.1)then
        call dljpotHe(x,y,z,pot,dpotx,dpoty,dpotz,dmax,'calc')
     endif
     if(dabs(pot/e0).gt.trj%sw1)then
  !****
        do while(dabs(pot/e0).gt.trj%sw1)
           id2=id2+10
           y=dble(id2)*1.0d-10
           if(mpiP%igas.eq.2)then
              call dljpotN2(x,y,z,pot,dpotx,dpoty,dpotz,dmax,'calc')
           elseif(mpiP%igas.eq.1)then
              call dljpotHe(x,y,z,pot,dpotx,dpoty,dpotz,dmax,'calc')
           endif
        enddo
        do while(dabs(pot/e0).lt.trj%sw1)
           id2=id2-1
           y=dble(id2)*1.0d-10
           if(mpiP%igas.eq.2)then
              call dljpotN2(x,y,z,pot,dpotx,dpoty,dpotz,dmax,'calc')
           elseif(mpiP%igas.eq.1)then
              call dljpotHe(x,y,z,pot,dpotx,dpoty,dpotz,dmax,'calc')
           endif
        enddo
     else
        do while(dabs(pot/e0).lt.trj%sw1)
           id2=id2-1
           y=dble(id2)*1.0d-10
           if(mpiP%igas.eq.2)then
              call dljpotN2(x,y,z,pot,dpotx,dpoty,dpotz,dmax,'calc')
           elseif(mpiP%igas.eq.1)then
              call dljpotHe(x,y,z,pot,dpotx,dpoty,dpotz,dmax,'calc')
           endif
           if((id2.lt.iymin).and.(icheck.eq.0)) then
  !        if(it.eq.1) write(8,621)
              ang=0.d0
              erat=1.000d0
              ireturn=1
           endif
        enddo
  !****
     endif
     if(ireturn.eq.1)return
     y=dble(id2)*1.0d-10
     etot=e0+pot
  !      if(it.eq.1) write(8,622) y*1.0d10
     d1=y
  !***********************************
  !     initial coordinates and momenta
  !
     w(1)=x
     w(2)=vx*const%mu
     w(3)=y
     w(4)=vy*const%mu
     w(5)=z
     w(6)=vz*const%mu
     tim=0.d0
  
  !      if(it.eq.1) write(8,623)
  !
  !     initialize the time derivatives of the coordinates and momenta
  !
     call deriv(w,dw,pot,dpotx,dpoty,dpotz,dmax)
     ns=0
     nw=0
     l=0
  !**********************************
     istop=0
     ireturn=0
     do while(istop.eq.0)
        call diffeq(l,tim,dt,w,dw,pot,dmax)
        nw=nw+1
        if (nw.eq.trj%inwr)then
           ns=ns+nw
           nw=0
  !     check if trajectory has become "lost" (too many steps)
  !
           if(ns.gt.30000) then
              if(mpiP%imyrank.eq.0)then
                 write(8,105) b,v
              else
                 write(1000+mpiP%imyrank,105) b,v
              endif
              ang=pi/2.d0
              e=0.5d0*const%mu*(dw(1)**2+dw(3)**2+dw(5)**2)
              erat=(e+pot)/etot
              istep=ns
              ireturn=1
              istop=1
           else
  !
  !     check if the trajectory is finished
  !
              istop=1
              if(dmax.lt.const%romax)istop=0
              if(dabs(pot/e0).gt.trj%sw2.and.dt.eq.dt1.and.istop.eq.1) then
                 dt=dt2
                 l=0
              endif
              if(dabs(pot/e0).lt.trj%sw2.and.dt.eq.dt2.and.istop.eq.1) then
                 dt=dt1
                 l=0
              endif
              if(dabs(pot/e0).gt.trj%sw1)istop=0
              if(ns.lt.50)istop=0
           endif
  !
        endif
     enddo
     if(ireturn.eq.1)return
  !**********************************
     istep=ns
  !
  !     determine scattering angle
  !
     if(dw(1).gt.0.d0) then
        num=dw(3)*(-v)
        den=v*dsqrt(dw(1)**2+dw(3)**2+dw(5)**2)
        ang=dacos(num/den)
     endif
     if(dw(1).lt.0.d0) then
        num=dw(3)*(-v)
        den=v*dsqrt(dw(1)**2+dw(3)**2+dw(5)**2)
        ang=(-dacos(num/den))
     endif
  !
  !     check for energy conservation
  !
     e=0.5d0*const%mu*(dw(1)**2+dw(3)**2+dw(5)**2)
     erat=(e+pot)/etot
  622 format(1x,'trajectory start position =',1pe11.4)
  621 format(1x,'trajectory not started - potential too small')
  623 format(//1x,'trajectory ns, x,  y,  z,  kin e, dt,    tot e',/1x,&
     &'               vx, vy, vz, pot e, pot/e0'/)
  105 format(1x,'trajectory lost: b =',1pe11.4,' v =',e11.4)
  108 format(/1x,'specific trajectory parameters',//1x,&
     &'v =',1pe11.4,4x,'b =',e11.4)
  107 format(1x,'time steps, dt1 =',1pe11.4,' dt2 =',e11.4)
     return
  !
  !
  end
  !
  !     *******************************************************
  !
  subroutine diffeq(l,tim,dt,w,dw,pot,dmax)
  !
  !     Integration subroutine - uses 5th order runge-kutta-gill to
  !     initiate and 5th order adams-moulton predictor-corrector to
  !     propagate. Parameter l is initially set to zero and then
  !     incremented to tell the subroutine when to switch between
  !     integration methods. DIFFEQ calls subroutine DERIV to define
  !     the equations of motion to be integrated.
  !
     use mpi
     use xtb_mctc_accuracy
     implicit none
     !include 'mpif.h'

     integer :: ns,nw,l
     integer :: istep

     !real(wp) :: x, b, v
     real(wp) :: tim, dt, pot, dmax

     real(wp) :: w(6),dw(6)

     dimension a(4),b(4),c(4),ampc(5),amcc(4),array(6,40),savw(40),&
     &savdw(40),q(40)

     data a/0.50d0,0.292893218814d0,1.70710678118d0,0.1666666666667d0/
     data b/2.0d0,1.0d0,1.0d0,2.0d0/
     data c/-0.5d0,-0.292893218814d0,-1.70710678118d0,-0.5d0/
     data ampc/-0.111059153612d0,0.672667757774d0,-1.70633621697d0,&
     &2.33387888707d0,-1.8524668225d0/
     data amcc/0.0189208128941d0,-0.121233356692d0,0.337771548703d0,&
     &-0.55921513665d0/
     data var,cvar,acst/2.97013888888d0,0.990972222222d0,&
     &0.332866152768d0/
     save hvar,hcvar
  !
     if(l.ge.0)then
  !      if (l) 4,1,3
        if(l.eq.0)then
           do 2 j=1,6
              q(j)=0.0
  2        continue
           hvar=dt*var
           hcvar=dt*cvar
           dt=0.5*dt
        endif
        l=l+1
  !
  !     This is the runge-kutta-gill part...the steps are broken up into
  !     half stekps to improve accuracy.
  !
        k=0
        istop=0
        do while(istop.eq.0)
           do 17 j=1,4
              if (((-1)**j).gt.0) tim=tim+0.5*dt
              call deriv(w,dw,pot,dpotx,dpoty,dpotz,dmax)
              do 7 i=1,6
                 dw(i)=dt*dw(i)
                 r=a(j)*(dw(i)-b(j)*q(i))
                 w(i)=w(i)+r
                 q(i)=q(i)+3.0*r+c(j)*dw(i)
  7           continue
  17       continue
           call deriv(w,dw,pot,dpotx,dpoty,dpotz,dmax)
           if(k.gt.0)istop=1
           if(k.le.0)k=1
        enddo
  !********
        if((l-6).ge.0)then
           l=-1
           dt=2.0*dt
        else
           do 8 j=1,6
              array(l,j)=dw(j)
  8        continue
        endif
  !********
     else
  !     This is the adams-moulton predictor-corrector part.
  !
        do 10 j=1,6
           savw(j)=w(j)
           savdw(j)=dw(j)
           array(6,j)=savdw(j)
           do 9 i=1,5
              array(6,j)=array(6,j)+ampc(i)*array(i,j)
  9        continue
           w(j)=array(6,j)*hvar+w(j)
  10    continue
        tim=tim+dt
        call deriv(w,dw,pot,dpotx,dpoty,dpotz,dmax)
        do 12 j=1,6
           array(6,j)=acst*dw(j)
           do 11 i=1,4
              array(i,j)=array(i+1,j)
              array(6,j)=array(i,j)*amcc(i)+array(6,j)
  11       continue
           array(5,j)=savdw(j)
           w(j)=savw(j)+hcvar*(array(5,j)+array(6,j))
  12    continue
        call deriv(w,dw,pot,dpotx,dpoty,dpotz,dmax)
     endif
  !********************
     return
  end
  !
  !     ***************************************************************
  !
  subroutine deriv(w,dw,pot,dpotx,dpoty,dpotz,dmax)
  !
  !     Defines Hamilton's equations of motion as the time derivatives
  !     of the coordinates and momenta.
  !
      use mpi
      use get_commons, only: mpi_parameters, constants
      use xtb_mctc_accuracy
     implicit none
     !include 'mpif.h'
     real(wp) :: w(6),dw(6)
     real(wp) :: pot, dpotx, dpoty, dpotz, dmax


  type(mpi_parameters) :: mpiP
  type(constants) :: const
  !
  !     From Hamilton's equations, the time derivatives of the coordinates
  !     are the conjugates divided by the mass.
  !
     dw(1)=w(2)/const%mu
     dw(3)=w(4)/const%mu
     dw(5)=w(6)/const%mu
  !
  !     Hamilton's equations for the time derivatives of the momenta
  !     evaluated by using the coordinate derivatives together with the
  !     chain rule.
  !
     x=w(1)
     y=w(3)
     z=w(5)
  !
  !    These are analytical derivatives.
  !
  !
     if(mpiP%igas.eq.2)then
        call dljpotN2(x,y,z,pot,dpotx,dpoty,dpotz,dmax,'calc')
     elseif(mpiP%igas.eq.1)then
        call dljpotHe(x,y,z,pot,dpotx,dpoty,dpotz,dmax,'calc')
     endif
     dw(2)=-(dpotx)
     dw(4)=-(dpoty)
     dw(6)=-(dpotz)
  !
  !
     return
  end
  !
  !     ***************************************************************
  !
  subroutine mobil2 (t,mob,cs,sdevpc)
  !
  !     Subroutine to determine average mobility by trajectory method.
  !     All integrations Monte Carlo (except over velocity).
  !
   use xtb_mctc_accuracy
   use constants
   use mpi
   use get_commons, only: mpi_parameters, constants, coordinates, angles, trajectory, &
     & read_inp
   implicit none
   !  implicit double precision (a-h,m-z)
   !include 'mpif.h'

   integer :: i, iatom, ibst, ic, icc, ifinish, ig, ihold, im, im2, ip, ir, irn
   integer :: istart, istep,  istop, it, iu2, iu3

   real(wp) :: t, mob, cs, sdevpc
   real(wp) :: b2max(100),parab2max(100),temp1,temp2
   real(wp) :: abc_time
   real(wp) ::  d1, dbst2, dbst22, ddd, delta, delta_time, dgst, dmax
   real(wp) :: dpotx, dpoty, dpotz
   real(wp) :: e_time, emaxx, erat, f, gst, gst2, gstt, hold, hold1, hold2, hold3
   real(wp) :: mom11st, mom12st, mom13st, mom22st, pot, r, r00x, rmax, rmaxx
   real(wp) :: rnb, rxy, rzy, sdom11st, sterr, sum, sum1, sum2, temp, term, tst, tst3
   real(wp) :: u2, v, w, x, y, z 

   real(wp) :: pgst(100),wgst(100)
   real(wp) :: q1st(100),q2st(100),cosx(0:500)
   real(wp) :: om11st(100),om12st(100),om13st(100),om22st(100)
   real(wp) :: q1stt(100),q2stt(100)

   real(wp) :: ayst, b, best, cest, bst2

  type(mpi_parameters) :: mpiP
  type(constants) :: const
  type(coordinates) :: coord
  type(angles) :: ang
  type(trajectory) :: trj
  type(read_inp) :: inp
  !
  !****PATCH PATCH
     if(mpiP%imyrank.eq.0)then
        abc_time=MPI_WTIME()
        delta_time=abc_time-mpiP%s_time
        write(8,*)'Starting mobil2 in ',delta_time,' s from start'
     endif
  !****PATCH PATCH
     if(im2.eq.0)then
        if(mpiP%imyrank.eq.0)then
           write(8,631)
           write(8,603) trj%sw1,trj%sw2,trj%dtsf1,trj%dtsf2,trj%inwr,trj%ifail
        else
           write(1000+mpiP%imyrank,631)
           write(1000+mpiP%imyrank,603) trj%sw1,trj%sw2,trj%dtsf1,trj%dtsf2,trj%inwr,trj%ifail
        endif
     endif
     it=0
     iu2=0
  !
  !     determine maximum extent and orientate along x axis
  !
     if(im2.eq.0)then
        if(mpiP%imyrank.eq.0)then
           write(8,632)
        else
           write(1000+mpiP%imyrank,632)
        endif
     endif
     rmax=0.d0
     do 1000 iatom=1,inp%inatom
        r=dsqrt((coord%ox(iatom)*coord%ox(iatom))+(coord%oy(iatom)*coord%oy(iatom))+&
        &(coord%oz(iatom)*coord%oz(iatom)))
        if(r.gt.rmax) then
           rmax=r
           ihold=iatom
        endif
  1000 continue
  !
     rzy=dsqrt((coord%oz(ihold)*coord%oz(ihold))+(coord%oy(ihold)*coord%oy(ihold)))
     ang%phi=dacos(coord%oz(ihold)/rzy)
     ang%phi=ang%phi+(pi/2.d0)
     if(coord%oy(ihold).lt.0.d0) ang%phi=(2.d0*pi)-ang%phi
     ang%phi=(2.d0*pi)-ang%phi
     ang%theta=0.d0
     ang%gamma=0.d0
     call rotate
     rxy=dsqrt((coord%fx(ihold)*coord%fx(ihold))+(coord%fy(ihold)*coord%fy(ihold)))
     ang%gamma=dacos(coord%fx(ihold)/rxy)
     if(coord%fy(ihold).lt.0.d0) ang%gamma=(2.d0*pi)-ang%gamma
     ang%gamma=(2.d0*pi)-ang%gamma
     if(im2.eq.0) iu3=1
     if(ip.eq.1) iu2=1
     call rotate
     iu3=0
     iu2=0
     hold=coord%fx(ihold)/rmax
     if(hold.lt.0.9999999999d0.or.hold.gt.1.0000000001d0.or.&
     &coord%fy(ihold).gt.1.0d-20.or.coord%fz(ihold).gt.1.0d-20.or.&
     &coord%fy(ihold).lt.-1.0d-20.or.coord%fz(ihold).lt.-1.0d-20) then
        if(mpiP%imyrank.eq.0)write(8,601)
        do 1001 iatom=1,inp%inatom
           hold=dsqrt(coord%fx(iatom)**2+coord%fy(iatom)**2+coord%fz(iatom)**2)
           if(mpiP%imyrank.eq.0)then
              write(8,602) iatom,coord%fx(iatom),coord%fy(iatom),coord%fz(iatom),hold
           else
              write(1000+mpiP%imyrank,602) iatom,coord%fx(iatom),coord%fy(iatom),coord%fz(iatom),&
              &hold
           endif
  1001  continue
        if(mpiP%imyrank.eq.0)close(8)
        call MPI_FINALIZE(mpiP%ierr)
        stop
     endif
  !****PATCH PATCH
     if(mpiP%imyrank.eq.0)then
        abc_time=MPI_WTIME()
        delta_time=abc_time-mpiP%s_time
        write(8,*)'Got through hold check in ',delta_time
     endif
  !****PATCH PATCH
  !
  !     determine rmax,  and r00 along x, y, and z directions
  !
  !     if(ip.eq.1) write(8,689)
     irn=1000
     ddd=(rmax+const%romax)/dble(irn)
  !
     y=0.d0
     z=0.d0
     emaxx=0.d0
     istop=0
     do while(istop.eq.0)
        do 1101 ir=1,irn
           x=rmax+const%romax-(dble(ir)*ddd)
           if(mpiP%igas.eq.2)then
              call dljpotN2(x,y,z,pot,dpotx,dpoty,dpotz,dmax,'calc')
           elseif(mpiP%igas.eq.1)then
              call dljpotHe(x,y,z,pot,dpotx,dpoty,dpotz,dmax,'calc')
           endif
           if(pot.gt.0.d0)then
              istop=1
           else
              r00x=x
              if(pot.lt.emaxx) then
                 rmaxx=x
                 emaxx=pot
              endif
           endif
  1101  continue
     enddo
     if(im2.eq.0)then
        if(mpiP%imyrank.eq.0)then
           write(8,614) emaxx/xe,rmaxx*1.0d10,r00x*1.0d10
        else
           write(1000+mpiP%imyrank,614) emaxx/xe,rmaxx*1.0d10,r00x*1.0d10
        endif
     endif
  !****PATCH PATCH
     if(mpiP%imyrank.eq.0)then
        abc_time=MPI_WTIME()
        delta_time=abc_time-mpiP%s_time
        write(8,*)'determined rmaxx in ',delta_time,' s from start'
        write(8,*)'rmaxx=',rmaxx,', const%ro=',const%ro
     endif
  !****PATCH PATCH
  !
  !      x=0.d0
  !      z=0.d0
  !      emaxy=0.d0
  !      istop=0
  !      do while(istop.eq.0)
  !       do 1100 ir=1,irn
  !        y=rmax+const%romax-(dble(ir)*ddd)
  !        call dljpot(x,y,z,pot,dpotx,dpoty,dpotz,dmax,'calc')
  !        if(pot.gt.0.d0)then
  !         istop=1
  !        else
  !         r00y=y
  !         if(pot.lt.emaxy) then
  !          rmaxy=y
  !          emaxy=pot
  !         endif
  !        endif
  !1100   continue
  !      enddo
  !      if(im2.eq.0)then
  !       if(mpiP%imyrank.eq.0)then
  !        write(8,613) emaxy/xe,rmaxy*1.0d10,r00y*1.0d10
  !       else
  !        write(1000+mpiP%imyrank,613) emaxy/xe,rmaxy*1.0d10,r00y*1.0d10
  !       endif
  !      endif
  !
  !      x=0.d0
  !      y=0.d0
  !      emaxz=0.d0
  !      istop=0
  !      do while(istop.eq.0)
  !       do 1102 ir=1,irn
  !        z=rmax+const%romax-(dble(ir)*ddd)
  !        call dljpot(x,y,z,pot,dpotx,dpoty,dpotz,dmax,'calc')
  !        if(pot.gt.0.d0)then
  !         istop=1
  !        else
  !         r00z=z
  !         if(pot.lt.emaxz) then
  !          rmaxz=z
  !          emaxz=pot
  !         endif
  !        endif
  !1102   continue
  !      enddo
  !      if(im2.eq.0)then
  !       if(mpiP%imyrank.eq.0)then
  !        write(8,615) emaxz/xe,rmaxz*1.0d10,r00z*1.0d10
  !       else
  !        write(1000+mpiP%imyrank,615) emaxz/xe,rmaxz*1.0d10,r00z*1.0d10
  !       endif
  !      endif
  !
  !     set-up integration over gst
  !
     tst=xk*t/const%eo
     if(im2.eq.0)then
        if(mpiP%imyrank.eq.0)then
           write(8,600) tst
        else
           write(1000+mpiP%imyrank,600) tst
        endif
     endif
     tst3=tst*tst*tst
  !
     dgst=5.0d-7*6.d0*dsqrt(tst)
     gst=dgst
     sum=0.d0
     sum1=0.d0
     sum2=0.d0
  !
     do 2020 i=1,mpiP%inp
        sum1=sum1+dsqrt(dble(i))
  2020 continue
  !
     if(im2.eq.0)then
        if(mpiP%imyrank.eq.0)then
           write(8,611)
        else
           write(1000+mpiP%imyrank,611)
        endif
     endif
     do 2000 i=1,mpiP%inp
        hold1=dsqrt(dble(i))
        hold2=dsqrt(dble(i-1))
        sum2=sum2+hold2
        wgst(i)=hold1/sum1
        gstt=tst3*(sum2+(hold1/2.d0))/sum1
  !
        istop=0
        do while(istop.eq.0)
           sum=sum+(dexp(-gst*gst/tst)*gst*gst*gst*gst*gst*dgst)
           gst=gst+dgst
           if(sum.gt.gstt) pgst(i)=gst-(dgst/2.d0)
           istop=1
           if(sum.lt.gstt)istop=0
        enddo
  !
        hold1=dsqrt((pgst(i)*pgst(i)*const%eo)/(0.5d0*const%mu))
        hold2=0.5d0*const%mu*hold1*hold1/(xk*t)
        hold3=dexp(-pgst(i)*pgst(i)/tst)*pgst(i)**5.d0
        if(im2.eq.0)then
           if(mpiP%imyrank.eq.0)then
              write(8,610) pgst(i),wgst(i),hold1,hold2,&
              &hold3,sum/tst3
           else
              write(1000+mpiP%imyrank,610) pgst(i),wgst(i),hold1,hold2,&
              &hold3,sum/tst3
           endif
        endif
  2000 continue
  !****PATCH PATCH
     if(mpiP%imyrank.eq.0)then
        abc_time=MPI_WTIME()
        delta_time=abc_time-mpiP%s_time
        write(8,*)'determined gst in ',delta_time,' s from start'
     endif
  !****PATCH PATCH
  !
  !     determine b2max
  !*****
     do ig=1,100
        parab2max(ig)=0.0d0
        b2max(ig)=0.0d0
     enddo
  !*****
  !
  !      dbst2=1.d0
     dbst2=rmaxx/const%ro
     if(dbst2.lt.1.0d0)dbst2=1.0d0
     dbst22=dbst2/10.d0
     trj%cmin=0.0005d0
     if(im2.eq.0)then
        if(mpiP%imyrank.eq.0)then
           write(8,652) trj%cmin
        else
           write(1000+mpiP%imyrank,652) trj%cmin
        endif
     endif
  !********************************************
     istart=mpiP%inp-mpiP%inp_per_node*mpiP%imyrank
     ifinish=mpiP%inp-mpiP%inp_per_node*(mpiP%imyrank+1)+1
  !********************************************
     do 3030 ig=istart,ifinish,-1
        gst2=pgst(ig)*pgst(ig)
        v=dsqrt((gst2*const%eo)/(0.5d0*const%mu))
        ibst=nint(rmaxx/const%ro)-6
        if(ig.lt.istart)ibst=nint(b2max(ig+1)/dbst2)-6
        if(ibst.lt.0) ibst=0
  !*****************
        istop=0
        do while(istop.eq.0)
           bst2=dbst2*dble(ibst)
           b=const%ro*dsqrt(bst2)
           call gsang(v,b,erat,ang,d1,istep)
           cosx(ibst)=1.d0-dcos(ang)
  !       if(ip.eq.1) write(8,651) b,bst2,ang,cosx(ibst),erat
           if(ibst.ge.5)then
              if(cosx(ibst).lt.trj%cmin.and.cosx(ibst-1).lt.trj%cmin.and.&
              &cosx(ibst-2).lt.trj%cmin.and.cosx(ibst-3).lt.trj%cmin.and.&
              &cosx(ibst-4).lt.trj%cmin) istop=1
           endif
           if(istop.eq.0)then
              ibst=ibst+1
  !***Error checking
              if(ibst.gt.750) then
                 if(mpiP%imyrank.eq.0)then
                    write(8,653)
                    call flush(8)
                 else
                    write(1000+mpiP%imyrank,653)
                    call flush(1000+mpiP%imyrank)
                 endif
                 if(mpiP%imyrank.eq.0)close(8)
                 call MPI_FINALIZE(mpiP%ierr)
                 stop
              endif
           endif
  !***Printing
        enddo
        b2max(ig)=dble(ibst-5)*dbst2
        istop=0
        do while(istop.eq.0)
           b2max(ig)=b2max(ig)+dbst22
           b=const%ro*dsqrt(b2max(ig))
           call gsang(v,b,erat,ang,d1,istep)
           if((1.d0-dcos(ang)).le.trj%cmin)istop=1
        enddo
  !*****************
  3030 continue
  !********************************************
  !********Collect b2max values and broadcast them
  !******
     call MPI_REDUCE(b2max,parab2max,100,&
     &MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,mpiP%ierr)
     do ig=1,100
        b2max(ig)=parab2max(ig)
     enddo
     call MPI_BCAST(b2max,100,MPI_DOUBLE_PRECISION,0,&
     &MPI_COMM_WORLD,mpiP%ierr)
  !*******
  !****PATCH PATCH
     if(mpiP%imyrank.eq.0)then
        abc_time=MPI_WTIME()
        delta_time=abc_time-mpiP%s_time
        write(8,*)'determined b2max in ',delta_time,' s from start'
     endif
  !****PATCH PATCH
     if(im2.eq.0) then
        if(mpiP%imyrank.eq.0)then
           write(8,637)
           do ig=1,mpiP%inp
              write(8,630) pgst(ig),b2max(ig),const%ro*dsqrt(b2max(ig))*1.0d10
           enddo
        else
           write(1000+mpiP%imyrank,637)
           do ig=1,mpiP%inp
              write(1000+mpiP%imyrank,630)&
              &pgst(ig),b2max(ig),const%ro*dsqrt(b2max(ig))*1.0d10
           enddo
        endif
     endif
  !****
  !     Calculate Omega(1,1)*, Omega(1,2)*, Omega(1,3)*, and Omega(2,2)*
  !     by integrating Q(1)* or Q(2)* over all orientations, and initial
  !     relative velocities.
  !
  
     if(im2.eq.0)then
        if(mpiP%imyrank.eq.0)then
           write(8,672) mpiP%itn,mpiP%inp,mpiP%imp,mpiP%itn*mpiP%inp*mpiP%imp
        else
           write(1000+mpiP%imyrank,672) mpiP%itn,mpiP%inp,mpiP%imp,mpiP%itn*mpiP%inp*mpiP%imp
        endif
     endif
  !
  !      if(mpiP%imyrank.eq.0)then
  !       if(ip.eq.1) write(8,680)
  !      endif
     do 4011 ig=1,mpiP%inp
        q1st(ig)=0.d0
        q2st(ig)=0.d0
        q1stt(ig)=0.d0
        q2stt(ig)=0.d0
  4011 continue
  !
     do 4040 ic=1,mpiP%itn
        if(mpiP%imyrank.eq.0)then
           if(ip.eq.1) write(8,681) ic
           e_time=MPI_WTIME()
           write(8,681) ic
           write(8,*)'in',e_time-mpiP%s_time,'s'
           call flush(8)
        else
           if(ip.eq.1) write(1000+mpiP%imyrank,681) ic
           e_time=MPI_WTIME()
           write(1000+mpiP%imyrank,681) ic
           write(1000+mpiP%imyrank,*)'in',e_time-mpiP%s_time,'s'
           call flush(1000+mpiP%imyrank)
        endif
        om11st(ic)=0.d0
        om12st(ic)=0.d0
        om13st(ic)=0.d0
        om22st(ic)=0.d0
        om11stt(ic)=0.d0
        om12stt(ic)=0.d0
        om13stt(ic)=0.d0
        om22stt(ic)=0.d0
  !******
        do 4010 ig=1,mpiP%inp
           gst2=pgst(ig)*pgst(ig)
           v=dsqrt((gst2*const%eo)/(0.5d0*const%mu))
  !        if(mpiP%imyrank.eq.0)then
  !         if(ip.eq.1) write(8,682) ic,ig,gst2,v
  !        endif
           temp1=0.d0
           temp2=0.d0
  !****************************************
  !***************************************
           do 4000 im=1,mpiP%imp_per_node
              rnb=xrand()
              call rantate
              bst2=rnb*b2max(ig)
              b=const%ro*dsqrt(bst2)
              call gsang(v,b,erat,ang,d1,istep)
              hold1=1.d0-dcos(ang)
              hold2=dsin(ang)*dsin(ang)
              temp1=temp1+(hold1*b2max(ig)/dble(mpiP%imp))
              temp2=temp2+(1.5d0*hold2*b2max(ig)/dble(mpiP%imp))
  4000     continue
           om11st(ic)=om11st(ic)+(temp1*wgst(ig))
           om12st(ic)=om12st(ic)+(temp1*pgst(ig)*pgst(ig)*wgst(ig)*&
           &(1.d0/(3.d0*tst)))
           om13st(ic)=om13st(ic)+(temp1*(pgst(ig)**4)*wgst(ig)*&
           &(1.d0/(12.d0*tst*tst)))
           om22st(ic)=om22st(ic)+(temp2*pgst(ig)*pgst(ig)*wgst(ig)*&
           &(1.d0/(3.d0*tst)))
           q1st(ig)=q1st(ig)+temp1
           q2st(ig)=q2st(ig)+temp2
  4010  continue
  !*
  4040 continue
  !***********************************
  !***From here node 0 can take it..
  !***********************************
     call MPI_REDUCE(om11st,om11stt,100,MPI_DOUBLE_PRECISION,&
     &MPI_SUM,0,MPI_COMM_WORLD,mpiP%ierr)
     call MPI_REDUCE(om12st,om12stt,100,MPI_DOUBLE_PRECISION,&
     &MPI_SUM,0,MPI_COMM_WORLD,mpiP%ierr)
     call MPI_REDUCE(om13st,om13stt,100,MPI_DOUBLE_PRECISION,&
     &MPI_SUM,0,MPI_COMM_WORLD,mpiP%ierr)
     call MPI_REDUCE(om22st,om22stt,100,MPI_DOUBLE_PRECISION,&
     &MPI_SUM,0,MPI_COMM_WORLD,mpiP%ierr)
     call MPI_REDUCE(q1st,q1stt,100,MPI_DOUBLE_PRECISION,&
     &MPI_SUM,0,MPI_COMM_WORLD,mpiP%ierr)
     call MPI_REDUCE(q2st,q2stt,100,MPI_DOUBLE_PRECISION,&
     &MPI_SUM,0,MPI_COMM_WORLD,mpiP%ierr)
  !
  !     calculate running averages
  !
  !*******************************
     if(mpiP%imyrank.eq.0)then
        do ic=1,mpiP%itn
           om11st(ic)=om11stt(ic)
           om12st(ic)=om12stt(ic)
           om13st(ic)=om13stt(ic)
           om22st(ic)=om22stt(ic)
        enddo
        do ig=1,mpiP%inp
           q1st(ig)=q1stt(ig)
           q2st(ig)=q2stt(ig)
        enddo
  !***************************
        hold1=0.d0
        hold2=0.d0
        if(im2.eq.0) write(8,685)
        do 4041 icc=1,mpiP%itn
           temp=1.d0/(const%mconst/(dsqrt(t)*om11st(icc)*pi*const%ro*const%ro))
           hold1=hold1+om11st(icc)
           hold2=hold2+temp
           if(im2.eq.0) write(8,622) icc,om11st(icc)*pi*const%ro*const%ro*1.d20,&
           &hold1*pi*const%ro*const%ro*1.d20/dble(icc),temp,hold2/dble(icc)
  4041  continue
  !
        if(im2.eq.0) then
           write(8,675)
           do 4012 ig=1,mpiP%inp
              write(8,676) pgst(ig)*pgst(ig),wgst(ig),q1st(ig)/dble(mpiP%inp)
  4012     continue
        endif
  !
        mom11st=0.d0
        mom12st=0.d0
        mom13st=0.d0
        mom22st=0.d0
        do 4050 ic=1,mpiP%itn
           mom11st=mom11st+om11st(ic)
           mom12st=mom12st+om12st(ic)
           mom13st=mom13st+om13st(ic)
           mom22st=mom22st+om22st(ic)
  4050  continue
        mom11st=mom11st/dble(mpiP%itn)
        mom12st=mom12st/dble(mpiP%itn)
        mom13st=mom13st/dble(mpiP%itn)
        mom22st=mom22st/dble(mpiP%itn)
        sdom11st=0.d0
        do 4060 ic=1,mpiP%itn
           hold=mom11st-om11st(ic)
           sdom11st=sdom11st+(hold*hold)
  4060  continue
        sdom11st=dsqrt(sdom11st/dble(mpiP%itn))
        sterr=sdom11st/dsqrt(dble(mpiP%itn))
        if(im2.eq.0) write(8,674) mom11st,sdom11st,sterr
        cs=mom11st*pi*const%ro*const%ro
        sdevpc=100.d0*sdom11st/mom11st
  !
  !     Use omegas to obtain higher order correction factor to mobility
  !
        ayst=mom22st/mom11st
        best=((5.d0*mom12st)-(4.d0*mom13st))/mom11st
        cest=mom12st/mom11st
        term=((4.d0*ayst)/(15.d0))+(.5d0*((const%m2-const%m1)**2.d0)/(const%m1*const%m2))
        u2=term-(.08333d0*(2.4d0*best+1.d0)*(const%m1/const%m2))
        w=(const%m1/const%m2)
        delta=((((6.d0*cest)-5.d0)**2.d0)*w)/(60.d0*(1.d0+u2))
        f=1.d0/(1.d0-delta)
        if(im2.eq.0) write(8,673) f
        if(im2.eq.0) write(8,677) mom12st,mom13st,mom22st,u2,w,delta
        mob=(const%mconst*f)/(dsqrt(t)*cs)
        write(8,671) mob,1.d0/mob,cs*1.d20
        if(mpiP%imyrank.eq.0)then
           e_time=MPI_WTIME()
           write(8,*)'Job Completed in',e_time-mpiP%s_time,'s'
           call flush(8)
        endif
  !***********************************
        call flush(8)
  !***********************************
     endif
  !***********************************
  !***********************************
  671 format(//1x,'average (second order) TM mobility =',1pe11.4,&
     &/1x,'inverse average (second order) TM mobility =',e11.4,&
     &/1x,'average TM cross section =',e11.4)
  677 format(/1x,'omega*12 =',1pe11.4,2x,'omega*13 =',e11.4,2x,&
     &'omega*22 =',e11.4,/1x,'      u2 =',e11.4,2x,'       w =',&
     &e11.4,2x,'   delta =',e11.4)
  673 format(//1x,'f value for second order correction=',1pe11.4,/1x,&
     &'(integrations for second order correction are not',/1x,&
     &'accurate, check if correction becomes significant)')
  674 format(//1x,'mean OMEGA*(1,1) =',1pe11.4,/1x,&
     &'standard deviation =',&
     &e11.4,/1x,'standard error of mean =',e11.4)
  676 format(1x,1pe11.4,1x,e11.4,1x,e11.4)
  675 format(//1x,'average values for q1st',//5x,&
     &'gst2',8x,'wgst',8x,'q1st')
  622 format(1x,i3,4x,1pe11.4,4x,e11.4,4x,e11.4,4x,e11.4)
  685 format(/1x,'summary of mobility calculations',//1x,'cycle',&
     &5x,'cs/A^2',6x,'avge cs/A^2',8x,'Ko^-1',7x,'avge Ko^-1')
  620 format(/1x,'OMEGA(1,1)*=',1pe11.4,/)
  670 format(/1x,'v =',1pe11.4,5x,'q1st =',e11.4,/)
  684 format(1x,1pe11.4,7(e11.4))
  683 format(/5x,'b/A',8x,'ang',6x,'(1-cosX)',4x,'e ratio',4x,'ang%theta',&
     &7x,'ang%phi',7x,'ang%gamma')
  682 format(/1x,'ic =',i3,1x,'ig =',i4,1x,'gst2 =',1pe11.4,&
     &1x,'v =',e11.4)
  680 format(1x,'start mobility calculation')
  651 format(1x,1pe11.4,6(1x,e11.4))
  653 format(1x,'ibst greater than 750')
  637 format(/5x,'gst',11x,'b2max/const%ro2',9x,'b/A',/)
  630 format(1x,1pe11.4,5x,e11.4,5x,e11.4)
  672 format(//1x,'number of complete cycles (mpiP%itn) =',i6,/1x,&
     &'number of velocity points (mpiP%inp) =',i6,/1x,&
     &'number of random points (mpiP%imp) =',i6,/1x,&
     &'total number of points =',i7,/)
  681 format(/1x,'cycle number, ic =',i3)
  650 format(/1x,'gst2 =',1pe11.4,1x,'v =',e11.4,/6x,'b',&
     &10x,'bst2',7x,'X ang',7x,'cos(X)',6x,'e ratio')
  652 format(//1x,'set up b2 integration - integration over',&
     &' impact parameter',//1x,&
     &'minimum value of (1-cosX) =',1pe11.4,/)
  610 format(1x,1pe11.4,5(1x,e11.4))
  611 format(//1x,'set-up gst integration - integration over velocity',&
     &//5x,'pgst',8x,'wgst',9x,'v',9x,'ke/kt',7x,'gst^5*',5x,'frac of',&
     &/48x,'exp(gst^2/tst)',3x,'sum',/)
  600 format(/1x,'t*=',1pe11.4)
  615 format(1x,'along z axis const%emax =',1pe11.4,'eV rmax =',&
     &e11.4,'A r00 =',e11.4,'A',/)
  613 format(1x,'along y axis const%emax =',1pe11.4,'eV rmax =',&
     &e11.4,'A r00 =',e11.4,'A')
  614 format(1x,'along x axis const%emax =',1pe11.4,'eV rmax =',&
     &e11.4,'A r00 =',e11.4,'A')
  601 format(/1x,'Problem orientating along x axis',/)
  632 format(/1x,'maximum extent orientated along x axis')
  631 format(/1x,'mobility calculation by MOBIL2 (trajectory method)',/)
  603 format(1x,'global trajectory parameters',//1x,'trj%sw1 =',1pe11.4,7x,&
     &'trj%sw2 =',e11.4,/1x,'trj%dtsf1 =',e11.4,5x,'trj%dtsf2 =',e11.4,/1x,&
     &'trj%inwr =',i3,14x,'trj%ifail =',i5)
  602 format(1x,i4,5x,1pe11.4,3(5x,e11.4))
  689 format(/)
     return
  end
  !
  !     ***************************************************************
  !
  double precision function xrand() result(res)
    use xtb_mctc_accuracy
    use get_commons
    implicit none
     real(wp) :: rvec(10)
  !
  !     XRAND is a random number generator that uses RANLUX if mpiP%i5=1
  !     otherwise it uses the standard RAND subroutine available in
  !     FORTRAN 77. If RANLUX is used i1 contains the luxury level
  !     (1-4, 4 is the highest, 3 is default in RANLUX). mpiP%i2, mpiP%i3, and
  !     mpiP%i4 are seed integers. mpiP%i3 and mpiP%i4 will normally be zero. If the
  !     standard RAND subroutine is to be employed, mpiP%i2 contains the
  !     seed integer. RANLUX was downloaded from http://kumo.swcp.
  !     com/fortran/random2.f90.

  type(mpi_parameters) :: mpiP
  !
     mpiP%i6=mpiP%i6+1
  !
     if (mpiP%i5.eq.1) then
        call ranlux(rvec,1)
        !xrand=rvec(1)
        res=rvec(1)
        return
     else
        !xrand=rand()
        res=rand()
        return
     endif
  !
  end
  !
  !     ***************************************************************
  !
  
  SUBROUTINE RANLUX(RVEC,LENV)
  !         Subtract-and-borrow random number generator proposed by
  !         Marsaglia and Zaman, implemented by F. James with the name
  !         RCARRY in 1991, and later improved by Martin Luescher
  !         in 1993 to produce "Luxury Pseudorandom Numbers".
  !     Fortran 77 coded by F. James, 1993
  !
  !       references:
  !  M. Luscher, Computer Physics Communications  79 (1994) 100
  !  F. James, Computer Physics Communications 79 (1994) 111
  !
  !   LUXURY LEVELS.
  !   ------ ------      The available luxury levels are:
  !
  !  level 0  (p=24): equivalent to the original RCARRY of Marsaglia
  !           and Zaman, very long period, but fails many tests.
  !  level 1  (p=48): considerable improvement in quality over level 0,
  !           now passes the gap test, but still fails spectral test.
  !  level 2  (p=97): passes all known tests, but theoretically still
  !           defective.
  !  level 3  (p=223): DEFAULT VALUE.  Any theoretically possible
  !           correlations have very small chance of being observed.
  !  level 4  (p=389): highest possible luxury, all 24 bits chaotic.
  !
  !!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !!!!  Calling sequences for RANLUX:                                  ++
  !!!!      CALL RANLUX (RVEC, LEN)   returns a vector RVEC of LEN     ++
  !!!!                   32-bit random floating point numbers between  ++
  !!!!                   zero (not included) and one (also not incl.). ++
  !!!!      CALL RLUXGO(LUX,INT,K1,K2) initializes the generator from  ++
  !!!!               one 32-bit integer INT and sets Luxury Level LUX  ++
  !!!!               which is integer between zero and MAXLEV, or if   ++
  !!!!               LUX .GT. 24, it sets p=LUX directly.  K1 and K2   ++
  !!!!               should be set to zero unless restarting at a break++
  !!!!               point given by output of RLUXAT (see RLUXAT).     ++
  !!!!      CALL RLUXAT(LUX,INT,K1,K2) gets the values of four integers++
  !!!!               which can be used to restart the RANLUX generator ++
  !!!!               at the current point by calling RLUXGO.  K1 and K2++
  !!!!               specify how many numbers were generated since the ++
  !!!!               initialization with LUX and INT.  The restarting  ++
  !!!!               skips over  K1+K2*E9   numbers, so it can be long.++
  !!!!   A more efficient but less convenient way of restarting is by: ++
  !!!!      CALL RLUXIN(ISVEC)    restarts the generator from vector   ++
  !!!!                   ISVEC of 25 32-bit integers (see RLUXUT)      ++
  !!!!      CALL RLUXUT(ISVEC)    outputs the current values of the 25 ++
  !!!!                 32-bit integer seeds, to be used for restarting ++
  !!!!      ISVEC must be dimensioned 25 in the calling program        ++
  !!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     implicit double precision (a-h,o-z)
     DIMENSION RVEC(LENV)
     DIMENSION SEEDS(24), ISEEDS(24), ISDEXT(25)
     PARAMETER (MAXLEV=4, LXDFLT=3)
     DIMENSION NDSKIP(0:MAXLEV)
     DIMENSION NEXT(24)
     PARAMETER (TWOP12=4096., IGIGA=1000000000,JSDFLT=314159265)
     PARAMETER (ITWO24=2**24, ICONS=2147483563)
     SAVE NOTYET, I24, J24, CARRY, SEEDS, TWOM24, TWOM12, LUXLEV
     SAVE NSKIP, NDSKIP, IN24, NEXT, KOUNT, MKOUNT, INSEED
     INTEGER LUXLEV
     LOGICAL NOTYET
     DATA NOTYET, LUXLEV, IN24, KOUNT, MKOUNT /.TRUE., LXDFLT, 0,0,0/
     DATA I24,J24,CARRY/24,10,0./
  !                               default
  !  Luxury Level   0     1     2   *3*    4
     DATA NDSKIP/0,   24,   73,  199,  365 /
  !orresponds to p=24    48    97   223   389
  !     time factor 1     2     3     6    10   on slow workstation
  !                 1    1.5    2     3     5   on fast mainframe
  !
  !  NOTYET is .TRUE. if no initialization has been performed yet.
  !              Default Initialization by Multiplicative Congruential
     IF (NOTYET) THEN
        NOTYET = .FALSE.
        JSEED = JSDFLT
        INSEED = JSEED
  !         WRITE(8,'(A,I12)') ' RANLUX DEFAULT INITIALIZATION: ',JSEED
        LUXLEV = LXDFLT
        NSKIP = NDSKIP(LUXLEV)
        LP = NSKIP + 24
        IN24 = 0
        KOUNT = 0
        MKOUNT = 0
  !         WRITE(8,'(A,mpiP%i2,A,mpiP%i4)')  ' RANLUX DEFAULT LUXURY LEVEL =  ',
  !     +        LUXLEV,'      p =',LP
        TWOM24 = 1.
        DO 25 I= 1, 24
           TWOM24 = TWOM24 * 0.5
           K = JSEED/53668
           JSEED = 40014*(JSEED-K*53668) -K*12211
           IF (JSEED .LT. 0)  JSEED = JSEED+ICONS
           ISEEDS(I) = MOD(JSEED,ITWO24)
  25    CONTINUE
        TWOM12 = TWOM24 * 4096.
        DO 50 I= 1,24
           SEEDS(I) = REAL(ISEEDS(I))*TWOM24
           NEXT(I) = I-1
  50    CONTINUE
        NEXT(1) = 24
        I24 = 24
        J24 = 10
        CARRY = 0.
        IF (SEEDS(24) .EQ. 0.) CARRY = TWOM24
     ENDIF
  !
  !          The Generator proper: "Subtract-with-borrow",
  !          as proposed by Marsaglia and Zaman,
  !          Florida State University, March, 1989
  !
     DO 100 IVEC= 1, LENV
        UNI = SEEDS(J24) - SEEDS(I24) - CARRY
        IF (UNI .LT. 0.)  THEN
           UNI = UNI + 1.0
           CARRY = TWOM24
        ELSE
           CARRY = 0.
        ENDIF
        SEEDS(I24) = UNI
        I24 = NEXT(I24)
        J24 = NEXT(J24)
        RVEC(IVEC) = UNI
  !  small numbers (with less than 12 "significant" bits) are "padded".
        IF (UNI .LT. TWOM12)  THEN
           RVEC(IVEC) = RVEC(IVEC) + TWOM24*SEEDS(J24)
  !        and zero is forbidden in case someone takes a logarithm
           IF (RVEC(IVEC) .EQ. 0.)  RVEC(IVEC) = TWOM24*TWOM24
        ENDIF
  !        Skipping to luxury.  As proposed by Martin Luscher.
        IN24 = IN24 + 1
        IF (IN24 .EQ. 24)  THEN
           IN24 = 0
           KOUNT = KOUNT + NSKIP
           DO 90 ISK= 1, NSKIP
              UNI = SEEDS(J24) - SEEDS(I24) - CARRY
              IF (UNI .LT. 0.)  THEN
                 UNI = UNI + 1.0
                 CARRY = TWOM24
              ELSE
                 CARRY = 0.
              ENDIF
              SEEDS(I24) = UNI
              I24 = NEXT(I24)
              J24 = NEXT(J24)
  90       CONTINUE
        ENDIF
  100 CONTINUE
     KOUNT = KOUNT + LENV
     IF (KOUNT .GE. IGIGA)  THEN
        MKOUNT = MKOUNT + 1
        KOUNT = KOUNT - IGIGA
     ENDIF
     RETURN
  !
  !           Entry to input and float integer seeds from previous run
   ENTRY RLUXIN(ISDEXT)
     TWOM24 = 1.
     DO 195 I= 1, 24
        NEXT(I) = I-1
  195 TWOM24 = TWOM24 * 0.5
     NEXT(1) = 24
     TWOM12 = TWOM24 * 4096.
  !      WRITE(8,'(A)') ' FULL INITIALIZATION OF RANLUX WITH 25 INTEGERS:'
  !      WRITE(8,'(5X,5I12)') ISDEXT
     DO 200 I= 1, 24
        SEEDS(I) = REAL(ISDEXT(I))*TWOM24
  200 CONTINUE
     CARRY = 0.
     IF (ISDEXT(25) .LT. 0)  CARRY = TWOM24
     ISD = IABS(ISDEXT(25))
     I24 = MOD(ISD,100)
     ISD = ISD/100
     J24 = MOD(ISD,100)
     ISD = ISD/100
     IN24 = MOD(ISD,100)
     ISD = ISD/100
     LUXLEV = ISD
     IF (LUXLEV .LE. MAXLEV) THEN
        NSKIP = NDSKIP(LUXLEV)
  !          WRITE (8,'(A,mpiP%i2)') ' RANLUX LUXURY LEVEL SET BY RLUXIN TO: ',
  !     +                         LUXLEV
     ELSE  IF (LUXLEV .GE. 24) THEN
        NSKIP = LUXLEV - 24
  !          WRITE (8,'(A,mpiP%i5)') ' RANLUX P-VALUE SET BY RLUXIN TO:',LUXLEV
     ELSE
        NSKIP = NDSKIP(MAXLEV)
  !          WRITE (8,'(A,mpiP%i5)') ' RANLUX ILLEGAL LUXURY RLUXIN: ',LUXLEV
        LUXLEV = MAXLEV
     ENDIF
     INSEED = -1
     RETURN
  !
  !                    Entry to ouput seeds as integers
   ENTRY RLUXUT(ISDEXT)
     DO 300 I= 1, 24
        ISDEXT(I) = INT(SEEDS(I)*TWOP12*TWOP12)
  300 CONTINUE
     ISDEXT(25) = I24 + 100*J24 + 10000*IN24 + 1000000*LUXLEV
     IF (CARRY .GT. 0.)  ISDEXT(25) = -ISDEXT(25)
     RETURN
  !
  !                    Entry to output the "convenient" restart point
   ENTRY RLUXAT(LOUT,INOUT,K1,K2)
     LOUT = LUXLEV
     INOUT = INSEED
     K1 = KOUNT
     K2 = MKOUNT
     RETURN
  !
  !                    Entry to initialize from one or three integers
   ENTRY RLUXGO(LUX,INS,K1,K2)
     IF (LUX .LT. 0) THEN
        LUXLEV = LXDFLT
     ELSE IF (LUX .LE. MAXLEV) THEN
        LUXLEV = LUX
     ELSE IF (LUX .LT. 24 .OR. LUX .GT. 2000) THEN
        LUXLEV = MAXLEV
  !            WRITE (8,'(A,I7)') ' RANLUX ILLEGAL LUXURY RLUXGO: ',LUX
     ELSE
        LUXLEV = LUX
        DO 310 ILX= 0, MAXLEV
           IF (LUX .EQ. NDSKIP(ILX)+24)  LUXLEV = ILX
  310   CONTINUE
     ENDIF
     IF (LUXLEV .LE. MAXLEV)  THEN
        NSKIP = NDSKIP(LUXLEV)
  !         WRITE(8,'(A,mpiP%i2,A,mpiP%i4)') ' RANLUX LUXURY LEVEL SET BY RLUXGO :',
  !     +        LUXLEV,'     P=', NSKIP+24
     ELSE
        NSKIP = LUXLEV - 24
  !          WRITE (8,'(A,mpiP%i5)') ' RANLUX P-VALUE SET BY RLUXGO TO:',LUXLEV
     ENDIF
     IN24 = 0
     IF(INS.LT.0)INS=-INS
  !      IF (INS .LT. 0)  WRITE (6,'(A)')
  !     +   ' Illegal initialization by RLUXGO, negative input seed'
     IF (INS .GT. 0)  THEN
        JSEED = INS
  !        WRITE(8,'(A,3I12)') ' RANLUX INITIALIZED BY RLUXGO FROM SEEDS',
  !     +      JSEED, K1,K2
     ELSE
        JSEED = JSDFLT
  !        WRITE(8,'(A)')' RANLUX INITIALIZED BY RLUXGO FROM DEFAULT SEED'
     ENDIF
     INSEED = JSEED
     NOTYET = .FALSE.
     TWOM24 = 1.
     DO 325 I= 1, 24
        TWOM24 = TWOM24 * 0.5
        K = JSEED/53668
        JSEED = 40014*(JSEED-K*53668) -K*12211
        IF (JSEED .LT. 0)  JSEED = JSEED+ICONS
        ISEEDS(I) = MOD(JSEED,ITWO24)
  325 CONTINUE
     TWOM12 = TWOM24 * 4096.
     DO 350 I= 1,24
        SEEDS(I) = REAL(ISEEDS(I))*TWOM24
        NEXT(I) = I-1
  350 CONTINUE
     NEXT(1) = 24
     I24 = 24
     J24 = 10
     CARRY = 0.
     IF (SEEDS(24) .EQ. 0.) CARRY = TWOM24
  !        If restarting at a break point, skip K1 + IGIGA*K2
  !        Note that this is the number of numbers delivered to
  !        the user PLUS the number skipped (if luxury .GT. 0).
     KOUNT = K1
     MKOUNT = K2
     IF (K1+K2 .NE. 0)  THEN
        DO 500 IOUTER= 1, K2+1
           INNER = IGIGA
           IF (IOUTER .EQ. K2+1)  INNER = K1
           DO 450 ISK= 1, INNER
              UNI = SEEDS(J24) - SEEDS(I24) - CARRY
              IF (UNI .LT. 0.)  THEN
                 UNI = UNI + 1.0
                 CARRY = TWOM24
              ELSE
                 CARRY = 0.
              ENDIF
              SEEDS(I24) = UNI
              I24 = NEXT(I24)
              J24 = NEXT(J24)
  450      CONTINUE
  500   CONTINUE
  !         Get the right value of IN24 by direct calculation
        IN24 = MOD(KOUNT, NSKIP+24)
        IF (MKOUNT .GT. 0)  THEN
           IZIP = MOD(IGIGA, NSKIP+24)
           IZIP2 = MKOUNT*IZIP + IN24
           IN24 = MOD(IZIP2, NSKIP+24)
        ENDIF
  !       Now IN24 had better be between zero and 23 inclusive
        IF (IN24 .GT. 23) THEN
  !           WRITE (8,'(A/A,3I11,A,mpiP%i5)')
  !     +    '  Error in RESTARTING with RLUXGO:','  The values', INS,
  !     +     K1, K2, ' cannot occur at luxury level', LUXLEV
  !           IN24 = 0
        ENDIF
     ENDIF
     RETURN
  END
  ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !
  !     ***************************************************************
  !
  subroutine ncoord(unit,dchar,asymp)
  !
  !     Reads in a new set of coordinates
  !
     use mpi
     use get_commons, only: mpi_parameters, constants, coordinates, angles, & 
       trajectory, read_inp
     use xtb_mctc_accuracy
     implicit none
     !include 'mpif.h'

     integer :: igamma
     real(wp) :: asymp
     real(wp) :: pcharge(100000)

     character*30 unit,dchar,dummy
     integer, parameter :: ixlen=100000
     integer :: imass(100000)
     real(wp) :: xmass(100000)
  !
  type(mpi_parameters) :: mpiP
  type(constants) :: const
  type(coordinates) :: coord
  type(angles) :: ang
  type(trajectory) :: trj
  type(read_inp) :: inp

     read(9,'(a30)',end=100) dummy
  100 continue
  !
     do 2000 iatom=1,inp%inatom
        if(dchar.eq.'calc') then
           read(9,*) coord%fx(iatom),coord%fy(iatom),coord%fz(iatom),&
           &ximass,pcharge(iatom)
        else
           read(9,*) coord%fx(iatom),coord%fy(iatom),coord%fz(iatom),ximass
        endif
        imass(iatom)=nint(ximass)
        if(unit.eq.'au') then
           coord%fx(iatom)=coord%fx(iatom)*0.52917706d0
           coord%fy(iatom)=coord%fy(iatom)*0.52917706d0
           coord%fz(iatom)=coord%fz(iatom)*0.52917706d0
        endif
  2000 continue
  !
     mx=0.d0
     do 2021 iatom=1,inp%inatom
        mx=mx+xmass(iatom)
  2021 continue
  !
     if(mx.ne.const%m2) then
        write(8,624)
        if(mpiP%imyrank.eq.0)close(8)
        call MPI_FINALIZE(mpiP%ierr)
        stop
     endif
  !
  !
     fxo=0.d0
     fyo=0.d0
     fzo=0.d0
     do 2009 iatom=1,inp%inatom
        fxo=fxo+(coord%fx(iatom)*xmass(iatom))
        fyo=fyo+(coord%fy(iatom)*xmass(iatom))
        fzo=fzo+(coord%fz(iatom)*xmass(iatom))
  2009 continue
     fxo=fxo/const%m2
     fyo=fyo/const%m2
     fzo=fzo/const%m2
     do 2010 iatom=1,inp%inatom
        coord%fx(iatom)=(coord%fx(iatom)-fxo)*1.d-10*correct
        coord%fy(iatom)=(coord%fy(iatom)-fyo)*1.d-10*correct
        coord%fz(iatom)=(coord%fz(iatom)-fzo)*1.d-10*correct
  2010 continue
  !
     do 3000 iatom=1,inp%inatom
        coord%ox(iatom)=coord%fx(iatom)
        coord%oy(iatom)=coord%fy(iatom)
        coord%oz(iatom)=coord%fz(iatom)
  3000 continue
  !
  !     determine structural asymmetry parameter
  !
     ang%theta=0.d0
     asymp=0.d0
     do igamma=0,360,2
        do iphi=0,180,2
           ang%gamma=dble(igamma)/cang
           ang%phi=dble(iphi)/cang
           call rotate
           xyzsum=0.d0
           yzsum=0.d0
           do iatom=1,inp%inatom
              xyz=dsqrt(coord%fx(iatom)**2+coord%fy(iatom)**2+coord%fz(iatom)**2)
              yz=dsqrt(coord%fy(iatom)**2+coord%fz(iatom)**2)
              xyzsum=xyzsum+xyz
              yzsum=yzsum+yz
           end do
           hold=((pi/4.d0)*xyzsum)/yzsum
           if(hold.gt.asymp) asymp=hold
         end do
      end do
  !
  624 format(1x,'masses do not add up')
  !
     return
  end
  !
  !     ***************************************************************
  !
