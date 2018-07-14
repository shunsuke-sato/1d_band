module global_variables

! parameter
  real(8),parameter :: pi=3.14159265358979323846d0
  complex(8),parameter :: zI=(0d0,1d0)
  real(8),parameter :: a_B=0.529177d0,Ry=13.6058d0

!Grid
  real(8) :: length_x,length_y,length_z,square_yz,H,dKx
  integer :: Nx,Nk,NB,NBocc
  real(8),allocatable :: Lx(:),Kx(:)
  real(8) :: dt
  integer :: Nt

!Common
  complex(8),allocatable :: zpsi(:,:,:),zpsi_GS(:,:,:),ztpsi(:),zhtpsi(:)
  real(8),allocatable :: Veff(:)
  real(8) :: V0,V1,V2,V3
  real(8),allocatable :: spe(:,:),mass_e(:,:)
  real(8) :: L_coef(-2:2),G_coef(-2:2)
  integer :: Nelec
  

! Electro-magnetic field
  real(8) Acx
  real(8),allocatable :: Ac(:),curr(:)
  real(8) :: Tp_fs1,omega_ev1,IWcm2_1,Tp_fs2,omega_ev2,IWcm2_2,T1_T2_fs
  real(8) :: f0_VAA_1
  real(8) :: phi_CEP_1,phi_CEP_2
  character(50) :: filename,filename_dm

  real(8) :: n_ex,e_ex

  complex(8),allocatable :: zWaf(:),ekx(:,:)
  real(8),allocatable :: xWaf(:),VWaf(:)
  
  complex(8),allocatable :: zdensity_matrix(:,:),zdensity_matrix0(:),zdensity_matrix0_GS(:)


end module global_variables
!--------------------------------------------------------------------------------------------------
program main
  use global_variables
  implicit none
  integer :: iter
  integer :: ik,ib,ix1,ix2,ib2
  integer :: N_ex_elec
  real(8) :: s
  character(50) :: cnum
  real(8),allocatable :: nex_k(:,:,:)

  NB=6
  NBocc=2
  Nelec=2
  dt=0.08d0
  Nt=12000
  Acx=0d0
  write(*,*)'ok1'
  read(*,*)filename
  read(*,*)Nx
  read(*,*)NK
  read(*,*)length_x,length_y,length_z
!  length_x=length_x/0.529d0
!  length_y=length_y/0.529d0
!  length_z=length_z/0.529d0
  read(*,*)V1,V2,V3
  read(*,*)dt
  read(*,*)Nt
  read(*,*)Tp_fs1
  read(*,*)omega_ev1
  read(*,*)f0_VAA_1
  read(*,*)phi_CEP_1

  write(*,*)'ok2'
  call preparation

  write(*,*)'ok3'
  call ground_state

!  stop
  write(*,*)'ok4'
  call input_Ac


!  stop

  open(10,file=trim(filename)//'_Ac_j.out')
  do iter=0,Nt

! density matrix
    
    Acx=Ac(iter)
    call current(iter)
    Acx=0.5d0*(Ac(iter) + Ac(iter+1))
    call dt_evolve    

    write(10,'(100e26.16)')dt*dble(iter),Ac(iter),curr(iter),-(Ac(iter+1)-Ac(iter-1))/dt*0.5d0
  end do
  close(10)


  allocate(nex_k(NB,NBocc,NK))
  do ik = 1,NK
    do ib = 1,NB
      do ib2 = 1,NBocc
        nex_k(ib,ib2,ik) = abs(sum(conjg(zpsi(:,ib2,ik))*zpsi_GS(:,ib,ik))*H)**2
      end do
    end do
  end do

  write(*,'(A)')'# of ecited electron (%)'
  write(*,'(999e26.16e3)')(dble(NBocc*NK)-sum(nex_k(1:2,:,:)))/dble(NBocc*NK)
  
!  open(10,file=trim(filename)//'_excited_elec.out')
!  do ik=1,NK
!    write(10,'(100e26.16)')kx(ik),(abs(sum(conjg(zpsi(:,1,ik))*zpsi_GS(:,ib,ik))*H)**2,ib=1,NB)
!  end do
!  close(10)

!  call  analyze
!  write(*,*)'Kx =',2d0*pi/(dble(Nk)*length_x)

end program main
 !--------------------------------------------------------------------------------------------------
subroutine preparation
  use global_variables
  implicit none
  integer :: ix,ik
  real(8) :: x
  real(8) :: lambda
  real(8) :: sigma_t,xp1,xp2,xp3,xp4,xm1,xm2,xm3,xm4


   H=length_x/dble(Nx)
   dKx=2d0*pi/(dble(Nk)*length_x)

   L_coef(0)=-30d0/12d0/(H**2)
   L_coef(1)=16d0/12d0/(H**2)
   L_coef(2)=-1d0/12d0/(H**2)
   L_coef(-1)=16d0/12d0/(H**2)
   L_coef(-2)=-1d0/12d0/(H**2)

   G_coef(0)=0d0
   G_coef(1)=8d0/12d0/H
   G_coef(2)=-1d0/12d0/H
   G_coef(-1)=-8d0/12d0/H
   G_coef(-2)=1d0/12d0/H

   allocate(Lx(Nx),Kx(NK))


   allocate(ztpsi(Nx),zhtpsi(Nx))
   allocate(Veff(Nx),spe(NB,NK),mass_e(NB,NK))
   allocate(zpsi(Nx,NBocc,NK),zpsi_GS(Nx,NB,NK))    
   allocate(curr(0:Nt))


   allocate(zWaf(Nx*NK),ekx(Nx*NK,NK),xWaf(Nx*NK))

   allocate(zdensity_matrix(Nx*NK,Nx*NK),zdensity_matrix0(Nx*NK))
   allocate(zdensity_matrix0_GS(Nx*NK))

   do ix=1,Nx
     Lx(ix)=dble(ix)*H -0.5d0*length_x
   end do


   do ik=1,NK
      kx(ik)=dKx*dble(ik)-pi/length_x-dkx/2d0
!      kx(ik)=dKx*dble(ik)-pi/length_x!-dkx/2d0      
   end do

   do ix=1,Nx*NK
     xWaf(ix)=dble(ix)*H
   end do
   
   do ik=1,NK
     do ix=1,Nx*NK
       ekx(ix,ik)=exp((0d0,1d0)*kx(ik)*xWaf(ix))
     end do
   end do


   sigma_t=1.5d0
   do ix=1,Nx

     Veff(ix) = -V1*(1d0+cos(2d0*pi*Lx(ix)/length_x))

!     Veff(ix)=V1*sin(pi/(length_x)*Lx(ix))**2 &
!       &+V2*sin(pi/(length_x)*Lx(ix))**4 &
!       &+V3*sin(pi/(length_x)*Lx(ix))**6 
!     Veff(ix)=V1*cos(pi/length_x*Lx(ix))**2
      !     Veff(ix)=-0.7d0*(1d0+tanh(Lx(ix)+0.8d0))*(1d0+tanh(-Lx(ix)+0.8d0))
!      xp1 = Lx(ix)+length_x
!      xm1 = Lx(ix)-length_x
!      xp2 = Lx(ix)+2d0*length_x
!      xm2 = Lx(ix)-2d0*length_x
!      xp3 = Lx(ix)+3d0*length_x
!      xm3 = Lx(ix)-3d0*length_x
!      xp4 = Lx(ix)+4d0*length_x
!      xm4 = Lx(ix)-4d0*length_x
! 0.1083
! 0.1085
!      Veff(ix)=0.1084d0*(&
!           Lx(ix)*exp(-0.5d0*(Lx(ix)/sigma_t)**2) &
!           +xp1*exp(-0.5d0*(xp1/sigma_t)**2) &
!           +xm1*exp(-0.5d0*(xm1/sigma_t)**2) &
!           +xp2*exp(-0.5d0*(xp2/sigma_t)**2) &
!           +xm2*exp(-0.5d0*(xm2/sigma_t)**2) &
!           +xp3*exp(-0.5d0*(xp3/sigma_t)**2) &
!           +xm3*exp(-0.5d0*(xm3/sigma_t)**2) &
!           +xp4*exp(-0.5d0*(xp4/sigma_t)**2) &
!           +xm4*exp(-0.5d0*(xm4/sigma_t)**2) &           
!      )


   end do

   return
 end subroutine preparation
 !--------------------------------------------------------------------------------------------------
 subroutine ground_state
   use global_variables
   implicit none
   complex(8) :: za(Nx,Nx),zv(Nx,Nx)
   real(8) :: ar(Nx,Nx),ai(Nx,Nx),vr(Nx,Nx),vi(Nx,Nx)
   real(8) :: e(Nx)
   integer :: ix,jx,jx2,ik,ib
   real(8) :: ce,ve,s
!LAPACK
  integer :: lwork
  complex(8),allocatable :: work_lp(:)
  real(8),allocatable :: rwork(:),w(:)
  integer :: info

  lwork=6*Nx
  allocate(work_lp(lwork),rwork(3*Nx-2),w(Nx))


   do ik=1,NK
     za=0d0
     do ix=1,Nx
       do jx=ix-2,ix+2
         jx2=jx
         if(jx<=0)jx2=Nx+jx
         if(jx>Nx)jx2=jx-Nx
         za(ix,jx2)=-0.5d0*(L_coef(jx-ix)+zI*2d0*(Kx(ik)+Acx)*G_coef(jx-ix))
       end do
       za(ix,ix)=za(ix,ix)+0.5d0*(kx(ik)+Acx)**2+Veff(ix)
     end do

!     ar=real(za)
!     ai=imag(za)
!!     call hermit_jacobi(Nx,ar,ai,vr,vi,e)
!     zv=vr+zI*vi
     Call zheev('V', 'U', Nx, za, Nx, w, work_lp, lwork, rwork, info)

     do ib=1,NB
       s=sum(abs(za(:,ib))**2)*H
       zpsi_GS(:,ib,ik)=za(:,ib)/sqrt(s)
     end do

     do ib=1,NB
       spe(ib,ik)=w(ib)
     end do

   end do

   do ib=1,NB
     mass_e(ib,1)=dKx**2/(spe(ib,NK)-2d0*spe(ib,1)+spe(ib,2))
     mass_e(ib,NK)=dKx**2/(spe(ib,1)-2d0*spe(ib,NK)+spe(ib,NK-1))
     do ik=2,NK-1
       mass_e(ib,ik)=dKx**2/(spe(ib,ik-1)-2d0*spe(ib,ik)+spe(ib,ik+1))
     end do
   end do

   open(10,file=trim(filename)//'_gswf.out')
   do ix=1,Nx
     write(10,'(100e26.16E3)')Lx(ix),abs(zpsi_GS(ix,1,1))**2,abs(zpsi_GS(ix,1,Nk/2))**2 &
       ,abs(zpsi_GS(ix,2,1))**2,abs(zpsi_GS(ix,2,Nk/2))**2
   end do
   close(10)
   

   open(10,file=trim(filename)//'_band_map.out')
     do ik=1,NK
       write(10,'(100e26.16E3)')kx(ik),(spe(ib,ik),ib=1,NB)
     end do
   close(10)

   open(10,file=trim(filename)//'_pot.out')
     do ix=1,Nx
       write(10,'(100e26.16E3)')Lx(ix),Veff(ix)
     end do
   close(10)

   open(10,file=trim(filename)//'_mass_map.out')
     do ik=1,NK
       write(10,'(100e26.16E3)')kx(ik),(mass_e(ib,ik),ib=1,NB)
     end do
   close(10)

   write(*,*)'Band gap       =',spe(NBocc+1,NK)-spe(NBocc,NK),(spe(NBocc+1,NK)-spe(NBocc,NK))*Ry*2d0
   write(*,*)'effective mass =',1d0/(1d0/mass_e(NBocc+1,NK)-1d0/mass_e(NBocc,NK)),mass_e(NBocc,NK),mass_e(NBocc+1,Nk)

   zpsi(:,1:NBocc,:)=zpsi_GS(:,1:NBocc,:)

   return
 end subroutine ground_state
 !--------------------------------------------------------------------------------------------------
 subroutine hpsi(ik)
   use global_variables
   implicit none
   integer :: ix,ik



   zhtpsi(1)=-0.5d0*( &
     &+L_coef(0)*ztpsi(1) &
     &+L_coef(1)*(ztpsi(2)+ztpsi(Nx)) &
     &+L_coef(2)*(ztpsi(3)+ztpsi(Nx-1))) &
     &-zI*(Kx(ik)+Acx)*( &
     &+G_coef(1)*(ztpsi(2)-ztpsi(Nx)) &
     &+G_coef(2)*(ztpsi(3)-ztpsi(Nx-1)))
   zhtpsi(2)=-0.5d0*( &
     &+L_coef(0)*ztpsi(2) &
     &+L_coef(1)*(ztpsi(3)+ztpsi(1)) &
     &+L_coef(2)*(ztpsi(4)+ztpsi(Nx))) &
     &-zI*(Kx(ik)+Acx)*( &
     &+G_coef(1)*(ztpsi(3)-ztpsi(1)) &
     &+G_coef(2)*(ztpsi(4)-ztpsi(Nx)))

   do ix=3,Nx-2
   zhtpsi(ix)=-0.5d0*( &
     &+L_coef(0)*ztpsi(ix) &
     &+L_coef(1)*(ztpsi(ix+1)+ztpsi(ix-1)) &
     &+L_coef(2)*(ztpsi(ix+2)+ztpsi(ix-2))) &
     &-zI*(Kx(ik)+Acx)*( &
     &+G_coef(1)*(ztpsi(ix+1)-ztpsi(ix-1)) &
     &+G_coef(2)*(ztpsi(ix+2)-ztpsi(ix-2)))
   end do

   zhtpsi(Nx-1)=-0.5d0*( &
     &+L_coef(0)*ztpsi(Nx-1) &
     &+L_coef(1)*(ztpsi(Nx)+ztpsi(Nx-2)) &
     &+L_coef(2)*(ztpsi(1)+ztpsi(Nx-3))) &
     &-zI*(Kx(ik)+Acx)*( &
     &+G_coef(1)*(ztpsi(Nx)-ztpsi(Nx-2)) &
     &+G_coef(2)*(ztpsi(1)-ztpsi(Nx-3)))
   zhtpsi(Nx)=-0.5d0*( &
     &+L_coef(0)*ztpsi(Nx) &
     &+L_coef(1)*(ztpsi(1)+ztpsi(Nx-1)) &
     &+L_coef(2)*(ztpsi(2)+ztpsi(Nx-2))) &
     &-zI*(Kx(ik)+Acx)*( &
     &+G_coef(1)*(ztpsi(1)-ztpsi(Nx-1)) &
     &+G_coef(2)*(ztpsi(2)-ztpsi(Nx-2)))


   zhtpsi(:)=zhtpsi(:)+(0.5d0*(kx(ik)+Acx)**2+Veff(:))*ztpsi(:)

   return
 end subroutine hpsi
 !--------------------------------------------------------------------------------------------------
 subroutine input_Ac
   use global_variables
   implicit none
   real(8) :: tpulse_1,omega_1,f0_1,tpulse_2,omega_2,f0_2,T1_T2
   integer :: iter
   real(8) :: tt

!   f0_1=5.338d-9*sqrt(IWcm2_1)      ! electric field in a.u.
   f0_1 = f0_VAA_1*a_B/(2d0*Ry)
   omega_1=omega_ev1/(2d0*Ry)  ! frequency in a.u.
   tpulse_1=Tp_fs1/0.02418d0 ! pulse duration in a.u.
!   tpulse_1=2d0*pi/omega_1 ! (single-cycle) pulse duration in a.u.

   allocate(Ac(-1:Nt+1)); Ac = 0d0
 ! pulse shape : A(t)=f0/omega*sin(Pi t/T)**2 *cos (omega t+phi_CEP) 
 ! pump laser
       do iter=0,Nt+1
         tt=iter*dt
!         if (tt<tpulse_1*0.5d0) then !CW
         if (tt<tpulse_1) then !Pulse
           Ac(iter)=-f0_1/omega_1*(cos(pi*(tt-0.5d0*tpulse_1)/tpulse_1))**4*cos(omega_1*(tt-0.5d0*tpulse_1)+phi_CEP_1*2d0*pi)
!           Ac(iter)=-f0_1/omega_1*(sin(pi*tt/tpulse_1))**2
!           Ac(iter)=-f0_1/omega_1*(sin(0.5d0*pi*tt/tpulse_1))**2
         else
!           Ac(iter)=-f0_1/omega_1*cos(omega_1*tt+phi_CEP_1) ! CW
           Ac(iter)=0d0 ! Pulse
         end if
       enddo

   return
 end subroutine input_Ac
 !--------------------------------------------------------------------------------------------------
 subroutine input_Ac_kick
   use global_variables
   implicit none
   real(8) :: tpulse_1,omega_1,f0_1,tpulse_2,omega_2,f0_2,T1_T2
   integer :: iter
   real(8) :: tt

   f0_1=5.338d-9*sqrt(IWcm2_1)      ! electric field in a.u.
   omega_1=omega_ev1/(2d0*Ry)  ! frequency in a.u.
   tpulse_1=Tp_fs1/0.02418d0 ! pulse duration in a.u.

   allocate(Ac(-1:Nt+1)); Ac = 0d0
   Ac = 0.005d0


   return
 end subroutine input_Ac_kick
 !--------------------------------------------------------------------------------------------------
 ! To Do*
 ! Lanchosz + exact diagonalization
 ! Direct Krylov + Taylor
 subroutine dt_evolve
   use global_variables
   implicit none
   integer,parameter :: n_Taylor   = 0
   integer,parameter :: n_Lanchosz = 1
   integer,parameter :: n_Krylov_direct = 2
   integer,parameter :: n_propagator = n_Taylor

   select case(n_propagator)
   case(n_Taylor)
      call dt_evolve_Taylor
   case default
      stop "Invalid propagator"
   end select
   

 end subroutine dt_evolve
  !--------------------------------------------------------------------------------------------------
 subroutine dt_evolve_Taylor
   complex(8) :: zfac
   integer :: ik,ib,iexp

   do ik=1,NK
     do ib=1,NBocc
       ztpsi(:)=zpsi(:,ib,ik)
       zfac=1d0
       do iexp=1,4
         zfac=zfac*(-zI*dt)/dble(iexp)
         call hpsi(ik)
         zpsi(:,ib,ik)=zpsi(:,ib,ik)+zfac*zhtpsi(:)
         ztpsi(:)=zhtpsi(:)
       end do
     end do
   end do

   return
 end subroutine dt_evolve_Taylor
 !--------------------------------------------------------------------------------------------------
 subroutine current(iter)
   use global_variables
   implicit none
   integer :: ik,ib
   integer :: ix,iter
   real(8) :: jxt,jxt2

   jxt=0d0

   do ik=1,NK
     do ib=1,NBocc
       jxt2=0d0
       jxt2=jxt2+aimag(conjg(zpsi(1,ib,ik))*( &
         &+G_coef(1)*(zpsi(2,ib,ik)-zpsi(Nx,ib,ik)) &
         &+G_coef(2)*(zpsi(3,ib,ik)-zpsi(Nx-1,ib,ik)))) &
         +aimag(conjg(zpsi(2,ib,ik))*( & 
         &+G_coef(1)*(zpsi(3,ib,ik)-zpsi(1,ib,ik)) &
         &+G_coef(2)*(zpsi(4,ib,ik)-zpsi(Nx,ib,ik))))
       do ix=3,Nx-2
         jxt2=jxt2+aimag(conjg(zpsi(ix,ib,ik))*( &
           &+G_coef(1)*(zpsi(ix+1,ib,ik)-zpsi(ix-1,ib,ik)) &
           &+G_coef(2)*(zpsi(ix+2,ib,ik)-zpsi(ix-2,ib,ik))))
       end do
       jxt2=jxt2+aimag( conjg(zpsi(Nx-1,ib,ik))*( &
         &+G_coef(1)*(zpsi(Nx,ib,ik)-zpsi(Nx-2,ib,ik)) &
         &+G_coef(2)*(zpsi(1,ib,ik)-zpsi(Nx-3,ib,ik)))) &
         +aimag(conjg(zpsi(Nx,ib,ik))*( &
         &+G_coef(1)*(zpsi(1,ib,ik)-zpsi(Nx-1,ib,ik)) &
         &+G_coef(2)*(zpsi(2,ib,ik)-zpsi(Nx-2,ib,ik))))

       jxt=jxt+jxt2*H+kx(ik)+Acx
     end do
   end do

   curr(iter)=jxt*dble(Nelec)/dble(NK*NBocc)/(length_x*length_y*length_z)
 !  curr(iter)=sum(abs(zpsi)**2)*H*dble(Nelec)/dble(NBocc*NK)

   return
 end subroutine current
 !--------------------------------------------------------------------------------------------------
 subroutine analyze
   use global_variables
   implicit none
   complex(8) :: Aw,jw,Ew
   complex(8) :: sigma,eps
   real(8) :: omega,domega,omega_s,omega_e
   integer :: it,iw
   real(8) :: f1,f2,f3,f4,tt

   domega=0.005d0/27.2d0
  open(10,file=trim(filename)//'eps.out')
  do iw=1,20000
    omega=iw*domega

    Aw=0d0
    jw=0d0
    do it=1,Nt
      tt=dt*dble(it)
      Aw=Aw+exp(zI*omega*tt)*Ac(it)*(1.d0-3.d0*(tt/(Nt*dt))**2+2.d0*(tt/(Nt*dt))**3)!*exp(-0.05d0/27.2d0*tt)
      jw=jw+exp(zI*omega*tt)*curr(it)*(1.d0-3.d0*(tt/(Nt*dt))**2+2.d0*(tt/(Nt*dt))**3)!*exp(-0.05d0/27.2d0*tt)
    end do
    Ew=-zI*omega*Aw
    sigma=jw/Ew
    eps=1+4d0*pi*zI*sigma/omega

    write(10,'(100e26.16E3)')omega*27.2d0,real(eps),aimag(eps),abs(Aw),abs(jw),abs(Ew)

  end do
  close(10)

  return
end subroutine analyze
!--------------------------------------------------------------------------------------------------
subroutine excited_elec
  use global_variables
  implicit none
   complex(8) :: za(Nx,Nx),zv(Nx,Nx)
   real(8) :: ar(Nx,Nx),ai(Nx,Nx),vr(Nx,Nx),vi(Nx,Nx)
   real(8) :: e(Nx)
   integer :: ix,jx,jx2,ik,ib,ib1,ib2
   real(8) :: ce,ve,s


   do ik=1,NK
     za=0d0
     do ix=1,Nx
       do jx=ix-2,ix+2
         jx2=jx
         if(jx<=0)jx2=Nx+jx
         if(jx>Nx)jx2=jx-Nx
         za(ix,jx2)=-0.5d0*(L_coef(jx-ix)+zI*2d0*(Kx(ik)+Acx)*G_coef(jx-ix))
       end do
       za(ix,ix)=za(ix,ix)+0.5d0*(kx(ik)+Acx)**2+Veff(ix)
     end do

     ar=real(za)
     ai=imag(za)
!     call hermit_jacobi(Nx,ar,ai,vr,vi,e)
     zv=vr+zI*vi

     do ib=1,NBocc
       s=sum(abs(zv(:,ib))**2)*H
       zpsi_GS(:,ib,ik)=zv(:,ib)/sqrt(s)
     end do

     do ib=1,NB
       spe(ib,ik)=e(ib)
     end do

   end do

   n_ex=0d0
   do ik=1,NK
     do ib1=1,NBocc
     do ib2=1,NBocc
       n_ex=n_ex+abs(sum(conjg(zpsi(:,ib1,ik))*zpsi_GS(:,ib2,ik))*H)**2
     end do
     end do
   end do
   n_ex=(dble(NK*NBocc)-n_ex)/dble(NK*NBocc)

  return
end subroutine excited_elec
!--------------------------------------------------------------------------------------------------
subroutine Wannier_state
  use global_variables
  implicit none
  integer :: ik1,ik2,ix1,ix2

  zWaf=0d0
  do ik1=1,NK
    ix2=0
    do ik2=1,NK
      do ix1=1,NX
        ix2=ix2+1
        zWaf(ix2)=zWaf(ix2)+zpsi(ix1,1,ik1)*ekx(ix2,ik1)
      end do
    end do
  end do


  return
end subroutine Wannier_state
!--------------------------------------------------------------------------------------------------
subroutine density_matrix
  use global_variables
  implicit none
  integer :: ik1,ix1,ix2

  zdensity_matrix(:,:)=0d0

  do ik1=1,NK
    do ix1=1,Nx*NK
      do ix2=1,Nx*NK
        zdensity_matrix(ix1,ix2)=zdensity_matrix(ix1,ix2) &
          &+zpsi(mod(ix1-1,Nx)+1,1,ik1)*ekx(ix1,ik1)*conjg(zpsi(mod(ix2-1,Nx)+1,1,ik1)*ekx(ix2,ik1))
      end do
    end do
  end do

  zdensity_matrix(:,:)=zdensity_matrix(:,:)/dble(NK)

  return
end subroutine density_matrix
!--------------------------------------------------------------------------------------------------
subroutine density_matrix0
  use global_variables
  implicit none
  integer :: ik1,ix1,ix2

  zdensity_matrix0(:)=0d0

  do ik1=1,NK
    do ix1=1,Nx*NK
      zdensity_matrix0(ix1)=zdensity_matrix0(ix1) &
        &+zpsi(mod(ix1-1,Nx)+1,1,ik1)*ekx(ix1,ik1)*conjg(zpsi(Nx,1,ik1))
    end do
  end do

  zdensity_matrix0(:)=zdensity_matrix0(:)/dble(NK)

  return
end subroutine density_matrix0
