      ! This is a 200 line 2-D Lid-Driven Cavity incompressible flow solver in Fortran. Details:
      !      - Staggered grid with uniform spacing: dx=dy=h
      !      - The PPE is solved, to enforce div(u)=0, with red-black Gauss-Seidel method
      !      - Time marching: Explicit Euler (1st order accurate)
      !      - Spatial discretization: 2nd order accurate, 1st order accurate boundary treatment
      !      - Pressure is treated purely implicitly
      !      - One ghost cell surrounds all interior cells
      !      - OpenMP parallelization is implemented
      !      - Compile command: "gfortran -fopenmp -std=gnu -O3 -g main.f90 -o main"
      ! Developed by Charlie Kawczynski. Email: charliekawczynski@gmail.com
      module PPE_solver
      implicit none
      integer,parameter :: cp = selected_real_kind(14)
      contains
      subroutine red_black(p,divU,fact,N,h,odd)
      implicit none
      integer,intent(in) :: N
      integer,dimension(2),intent(in) :: odd
      real(cp),intent(in) :: fact,h
      real(cp),dimension(N+2,N+2),intent(in) :: divU
      real(cp),dimension(N+2,N+2),intent(inout) :: p
      integer :: i,j
      !$OMP PARALLEL DO
      do j=2+odd(2),N+1,2; do i=2+odd(1),N+1,2
      p(i,j)=fact*(p(i+1,j)+p(i,j+1)+p(i-1,j)+p(i,j-1)-divU(i,j)*h)
      enddo; enddo
      !$OMP END PARALLEL DO
      end subroutine
      end module

      program main
      use PPE_solver
      implicit none
      real(cp),dimension(:,:),allocatable :: p,divU,u,ustar,v,vstar,E1_z,E2_z
      real(cp) :: h,dt,Re,hdt_inv,h_inv,dt_inv,Re_inv,h2,h2_inv,dt4,dt_h_inv,dV8
      real(cp) :: fact,KE,KE_old,KE_temp,dV,tol,delta_BL,max_divU,dt_h2Re_inv
      integer :: i,j,N,N_output,N_PPE,iter_PPE,iter,N_iter

      Re = 1000.0_cp                                ! Reynolds number
      dt = 2.0_cp*10.0_cp**(-3.0_cp)               ! time step
      N_iter = 10**5                               ! number of time steps
      N_PPE = 5                                    ! number of PPE iterations

      delta_BL = 1.0_cp/sqrt(Re)                   ! approximate boundary layer thickness
      N = maxval((/floor(4.0_cp/delta_BL),40/))    ! number of cells in each direction, N ~ Re**(3/4) for DNS
      h = 1.0_cp/real(N,cp)                        ! spatial step size (hard coded and uniform)
      N_output = 100                               ! output transient data every N_output time steps
      tol = 1.0_cp*10.0_cp**(-5.0_cp)              ! stops simulation when |KE-KE_old|/dt < tol
      if (N.gt.300) then; write(*,*) 'are you sure you want this large of a mesh? N=',N; stop 'done'; endif

      ! Initialize data
      KE_old=0.0_cp; KE=0.0_cp;fact=1.0_cp/4.0_cp;h_inv=1.0_cp/h;h2=h**2.0_cp;dt_h2Re_inv=dt/(h2*Re);dt_h_inv=dt/h
      Re_inv=1.0_cp/Re;dt_inv=1.0_cp/dt;dV=h**2.0_cp;hdt_inv=h_inv*dt_inv;h2_inv=1.0_cp/h2; dt4=0.25_cp*dt;dV8=dV/8.0_cp
      write(*,*) '3-D Lid-driven cavity flow. Re,N,h,dt,N_iter,N_PPE:'; write(*,*) Re,N,h,dt,N_iter,N_PPE
      write(*,*) ''
      allocate(u(N+3,N+2),ustar(N+3,N+2),E1_z(N+3,N+3),divU(N+2,N+2))
      allocate(v(N+2,N+3),vstar(N+2,N+3),E2_z(N+3,N+3),p(N+2,N+2))
      u = 0.0_cp; ustar = u; E1_z = 0.0_cp; p = 0.0_cp
      v = 0.0_cp; vstar = v; E2_z = 0.0_cp; divU = p

      call system('mkdir output')
      open(1,file='output/KE.dat')
      write(1,*) 'TITLE="KINETIC ENERGY VS TIME"'
      write(1,*) 'VARIABLES = "Time","Kinetic Energy","max(divU)"'
      write(1,*) 'ZONE DATAPACKING = POINT'

      do iter=1,N_iter ! Momentum equation
        !$OMP PARALLEL DO
        do j=2,N+2; do i=2,N+2 ! Advection term in form d/dxj (uj ui)
        ! d/dxj (uj ui) for i=j
        ustar(i,j) = -( (u(i,j)+u(i+1,j))**2 - (u(i-1,j)+u(i,j))**2 )*h_inv
        vstar(i,j) = -( (v(i,j)+v(i,j+1))**2 - (v(i,j-1)+v(i,j))**2 )*h_inv
        ! d/dxj (uj ui) for i≠j
        E1_z(i,j) = (u(i,j-1)+u(i,j))*(v(i-1,j)+v(i,j)) ! x (z edge)
        E2_z(i,j) = (v(i-1,j)+v(i,j))*(u(i,j-1)+u(i,j)) ! y (z edge)
        enddo; enddo
        !$OMP END PARALLEL DO
        !$OMP PARALLEL DO
        do j=2,N+1; do i=2,N+1 ! Advection term d/dxj (uj ui) (continued)
        ustar(i,j)=ustar(i,j)-(E1_z(i,j+1)-E1_z(i,j))*h_inv
        vstar(i,j)=vstar(i,j)-(E2_z(i+1,j)-E2_z(i,j))*h_inv
        enddo; enddo
        !$OMP END PARALLEL DO

        !$OMP PARALLEL DO
        do j=2,N+1; do i=2,N+1 ! Diffusion term
        ustar(i,j)=u(i,j)+(dt4*ustar(i,j)+dt_h2Re_inv*(u(i+1,j)-2.0_cp*u(i,j)+u(i-1,j) + &
                                                       u(i,j+1)-2.0_cp*u(i,j)+u(i,j-1)))
        vstar(i,j)=v(i,j)+(dt4*vstar(i,j)+dt_h2Re_inv*(v(i+1,j)-2.0_cp*v(i,j)+v(i-1,j) + &
                                                       v(i,j+1)-2.0_cp*v(i,j)+v(i,j-1)))
        enddo; enddo
        !$OMP END PARALLEL DO

        ustar( 2 ,:) = 0.0_cp; vstar(:, 2 ) = 0.0_cp ! Remove forcing on wall (kinematic)
        ustar(N+2,:) = 0.0_cp; vstar(:,N+2) = 0.0_cp ! Remove forcing on wall (kinematic)

        !$OMP PARALLEL DO
        do j=2,N+1; do i=2,N+1 ! Compute PPE source
        divU(i,j) = dt_inv*(ustar(i+1,j)-ustar(i,j)+vstar(i,j+1)-vstar(i,j))
        enddo; enddo
        !$OMP END PARALLEL DO

        do iter_PPE=1,N_PPE ! solve PPE: red-black Gauss-Seidel
          call red_black(p,divU,fact,N,h,(/0,0/)); call red_black(p,divU,fact,N,h,(/1,0/))
          call red_black(p,divU,fact,N,h,(/1,1/)); call red_black(p,divU,fact,N,h,(/0,1/))
          p( 1 ,:) = p( 2 ,:); p(:, 1 ) = p(:, 2 ) ! Apply p BCs
          p(N+2,:) = p(N+1,:); p(:,N+2) = p(:,N+1) ! Apply p BCs
        enddo

        !$OMP PARALLEL DO
        do j=2,N+1; do i=2,N+1 ! Pressure correction
          u(i,j) = ustar(i,j) - dt_h_inv*(p(i,j) - p(i-1,j))
          v(i,j) = vstar(i,j) - dt_h_inv*(p(i,j) - p(i,j-1))
        enddo; enddo
        !$OMP END PARALLEL DO

        u( 2 ,:) = 0.0_cp; v(:, 2 ) = 0.0_cp ! Apply u BCs (wall coincident)
        u(N+2,:) = 0.0_cp; v(:,N+2) = 0.0_cp ! Apply u BCs (wall coincident)

        ! Apply u BCs (ghost cells, including sliding lid)
        u(:,1)=-u(:,2);u(:,N+2)=2.0_cp-u(:,N+1)
        v(1,:)=-v(2,:);v(N+2,:)=-v(N+1,:)

        if (mod(iter+1,N_output).eq.1) then
          KE_temp = 0.0_cp
          !$OMP PARALLEL DO REDUCTION(+:KE_temp)
          do j=2,N+1; do i=2,N+1 ! Pressure correction
          KE_temp = KE_temp + (u(i,j)+u(i+1,j))**2+(v(i,j)+v(i,j+1))**2
          enddo; enddo
          !$OMP END PARALLEL DO
          !$OMP PARALLEL DO
          do j=2,N+1; do i=2,N+1 ! Compute divU
          divU(i,j) = h_inv*(u(i+1,j)-u(i,j)+v(i,j+1)-v(i,j))
          enddo; enddo
          !$OMP END PARALLEL DO
          KE_old = KE_temp; KE_old = KE_old*dV8
          if ((iter.gt.1).and.(abs(KE-KE_old)*dt_inv.lt.tol)) then; write(*,*) 'Exited early'; exit
          endif; max_divU = maxval(divU)
          write(1,*) iter*dt,KE,max_divU; flush(1)
          write(*,'(A43,F8.4,I8,1F9.3,3E15.4E2)') 't,iter,\% complete,max(divU),KE,|dKE/dt| = ',&
          iter*dt,iter,real(iter,cp)/real(N_iter,cp)*100.0_cp,max_divU,KE,abs(KE-KE_old)*dt_inv
          KE = KE_old
        endif
      enddo
      close(1) ! Close KE unit

      !$OMP PARALLEL DO
      do j=2,N+1; do i=2,N+1 ! Compute divU
      divU(i,j) = h_inv*(u(i+1,j)-u(i,j)+v(i,j+1)-v(i,j))
      enddo; enddo
      !$OMP END PARALLEL DO

      open(2,file='output/solution_c.dat') ! Export solution at center plane
      write(2,*) 'TITLE="3-D CUBIC LDC AT CELL CENTERS, Re=',Re,'dt=',dt,'N_PPE=',N_PPE,'"'
      write(2,*) 'VARIABLES = "x","y","u","v","p","divU"'
      write(2,*) 'ZONE, I=',N,',J=',N,' DATAPACKING = POINT'
      do j=2,N+1; do i=2,N+1
      write(2,*) (i-1)*h-0.5*h,(j-1)*h-0.5*h,0.5_cp*(u(i+1,j)+u(i,j)),0.5_cp*(v(i,j+1)+v(i,j)),p(i,j),divU(i,j)
      enddo; enddo; close(2)

      open(4,file='output/solution_n.dat') ! Export solution at nodes
      write(4,*) 'TITLE="3-D CUBIC LDC AT CELL CORNERS, Re=',Re,'dt=',dt,'N_PPE=',N_PPE,'"'
      write(4,*) 'VARIABLES = "x","y","u","v"'
      write(4,*) 'ZONE, I=',N+1,',J=',N+1,' DATAPACKING = POINT'
      do j=2,N+2; do i=2,N+2
      write(4,*) (i-2)*h,(j-2)*h,0.5_cp*(u(i,j-1)+u(i,j)),0.5_cp*(v(i-1,j)+v(i,j))
      enddo; enddo; close(4)

      deallocate(p,divU,u,v,ustar,vstar,E1_z,E2_z)
      end program