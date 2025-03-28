program isingchimera

  integer narg ! numero de argumentos por linea de comandos
  integer stat, cont, isor
  real*8 dranu, sumap, suma, sumah, sumath,ss,sss
  real*8, dimension (:), allocatable :: argumento
  character(50) arg
  character(50), dimension (:), allocatable :: name ! esto es el vector de nombres

  integer, dimension (:), allocatable :: ivec, ivec2,iindex, itime
  integer AllocateStatus, allo1, allo2, allo3, allo4,allo5, allo6, allo7, allo8, allo9, allo10, allo11, allo12, allo13
  integer allo15, allo16,allo2hk
  real*8, dimension (:), allocatable :: hg,th, over,s, phi,mean, meanlocal
  real*8, dimension (:,:), allocatable :: ee, ec,hk
  real*8, dimension (:,:), allocatable :: xi
  real*8, dimension (:,:), allocatable :: pesos
  real*8 meanfiring, ap, amean, amean2,pphi,p0
  open (unit = 10, file = "file.dat", status="unknown")
  open (unit = 11, file="meanfiring.dat", status="unknown")
  open (unit = 12, file="red.dat", status="unknown")
  open(unit=13, file="velocidad.dat", status="unknown")
  narg=9
  allocate (name(narg))
  allocate (argumento(narg))
  name(1)= 'numero de neuronas'
  name(2)= 'numero de vecinos por arriba electricas'
  name(3)= 'numero de vecinos por arriba quimicas'
  name(4)= 'semilla del numero aleatorio'
  name(5)= 'Temperatura'
  name (6)= 'numero de patrones'
  name (7)= 'acoplamiento electrico'
  name (8)= 'acoplamiento quimico'
  name (9)= 'Numero de pasos MC'

  ! entrada de parametros por linea de comandos
  !===========================================
  cont=0
  do i=1,narg
     call get_command_argument(number=i, value=arg, status=stat)
     IF (LEN_TRIM(arg) == 0) then
        write(*,*)'Faltan comandos:'
        do j=1,narg
           write(*,*) name(j)
        enddo
        stop
     end IF
     read(arg,*) argumento(i)
     ! write(*,*) argumento(i)
  end do
  !=============================================

  nneu=argumento(1) ! numero de neuronas
  npat=argumento(6) ! numero de patrones

  allocate(ee(nneu+1,nneu+1), stat=allo1)
  IF (allo1 /= 0) STOP "*** Not enough memory ***"

  allocate(ec(nneu+1,nneu+1), stat=allo12)
  IF (allo12 /= 0) STOP "*** Not enough memory ***"

  allocate(pesos(nneu+1,nneu+1), stat= allo2)
  IF (allo2 /= 0) STOP "*** Not enough memory ***"

  allocate(hk(nneu+1,nneu+1), stat= allo2hk)
  IF (allo2hk /= 0) STOP "*** Not enough memory ***"

  allocate(iindex(nneu+1), stat=allo3)
  IF (allo3 /= 0) STOP "*** Not enough memory ***"

  allocate(itime(nneu+1), stat=allo4)
  IF (allo4 /= 0) STOP "*** Not enough memory ***"

  allocate(hg(nneu+1), stat=allo5)
  IF (allo5 /= 0) STOP "*** Not enough memory ***"

  allocate(th(nneu+1), stat=allo6)
  IF (allo6 /= 0) STOP "*** Not enough memory ***"

  allocate(s(nneu+1), stat=allo7)
  IF (allo7 /= 0) STOP "*** Not enough memory ***"

  allocate(xi(npat+1,nneu+1), stat=allo8)
  IF (allo8 /= 0) STOP "*** Not enough memory ***"

  allocate(over(npat+1), stat=allo9)
  IF (allo9 /= 0) STOP "*** Not enough memory ***"


  nvec=argumento(2) ! numero de vecinos por arriba
  write(*,*)'numero de vecinos electricos', 2*nvec

  nnvec=2*nvec ! numero de vecinos de cada neurona


  nvec2=argumento(3) ! numero de vecinos por arriba
  write(*,*)'numero de vecinos quimicos', 2*nvec2

  nnvec2=2*nvec2 ! numero de vecinos de cada neurona

  allocate (ivec(nnvec+1),stat=allo10)
  IF (allo10 /= 0) STOP "*** Not enough memory ***"

  allocate (ivec2(nnvec2+1),stat=allo13)

  IF (allo13 /= 0) STOP "*** Not enough memory ***"

  allocate(phi(nneu+1), stat=allo11)
  IF (allo11 /= 0) STOP "*** Not enough memory ***"

  allocate(mean(nneu+1), stat=allo15)
  IF (allo15 /= 0) STOP "*** Not enough memory ***"

  allocate(meanlocal(nneu+1), stat=allo16)
  IF (allo16 /= 0) STOP "*** Not enough memory ***"

  !write(*,*) nneu, nvec

  ! elegimos los vecinos de las neuronas

  do l=1, nvec
     ivec(l)=l
     ivec(nvec+l)=-l
     !write(*,*)  ivec(l), ivec(nvec+l)
  end do

  do l=1, nvec2
     ivec2(l)=nvec+l
     ivec2(nvec2+l)=-nvec-l
     !write(*,*)  ivec2(l), ivec2(nvec2+l)
  end do



  ! semilla del numero aleatorio
  iseed=argumento(4)

  call dranini(iseed)

  !write(*,*) iseed

  !calentamiento del numero aleatorio

  do i=1, 10000
     xr=dranu()
     !write(*,*) xr
  end do



  temp=argumento(5)
  tempdifusion=1
  beta=1.d0/temp
  betadifusion=1./tempdifusion

  write(*,*) temp, beta


  ! inicializacion de pesos y matriz de connectividad

  do i=1,nneu
     do j=1, nneu
        ee(i,j)=0.d0
        ec(i,j)=0.d0
        !write(*,*) ee(i,j), i ,j
        pesos(i,j)=0.d0
     end do
     ee(i,i)=0.d0
     ec(i,i)=0.d0
     pesos(i,i)=0.d0
  end do

  gel=argumento(7)
  gch=argumento(8)




  ! umbrales
  do i=1, nneu

     th(i)=0d0

  end do


  !write(*,*) nnvec


  ! construccion de la red nolocal para las quimicas

  write(*,*) 'quimicas'

  do i=1,nneu
     do ii=1, nnvec2
        l=i+ivec2(ii)
        !write(*,*) 'esto',i, l
        l=mod(l,nneu)
        if (l.le.0) then
           l=l+nneu
        endif
        !write(*,*) 'neurona', i, 'vecino', l
        ec(l,i)=1.d0
        ec(i,l)=1.d0
     end do
     ec(i,i)=0.d0
  end do

  ! red regular nolocal generada

  write(*,*) 'Red regular no-local generada para las quimicas'




  ! construccion de la red nolocal para las electricas

  write(*,*) 'Electricas'

  do i=1,nneu
     do ii=1, nnvec
        l=i+ivec(ii)
        !write(*,*) 'esto',i, l
        l=mod(l,nneu)
        if (l.le.0) then
           l=l+nneu
        endif
        !write(*,*) 'neurona', i, 'vecino', l
        ee(l,i)=1.d0
        ee(i,l)=1.d0
     end do
     ee(i,i)=0.d0
  end do

  ! red regular nolocal generada

  write(*,*) 'Red regular no-local generada para las electricas'

  ! escribimos la red en un fichero
  do i=1,nneu
     do j=1,nneu
        !        if (ee(i,j).gt.0) then
        write(12,*) i , j, ee(i,j), ec(i,j)
        !       end if
     end do
  end do


  !stop

  ! condiciones iniciales de las neuronas

  do i=1, nneu

     xr=dranu()
     if (xr.gt.0.5) then
        s(i)=1.d0
     else
        s(i)=-1.d0
     end if

     !write(*,*) s(i), i

  end do


  ! otras condiciones iniciales


!  do i=1,nneu
!     aimenos=0.5*nneu-20
!     aimas=0.5*nneu+20

!     p0=0
!     s(i)=1.
!     xr=dranu()
!     if ((i.gt.aimenos).and.(i.lt.aimas)) then
!        p0=1
!     end if

!     if (p0>xr) then
!        s(i)=-1.
!     end if
!  end do

  ! campos locales iniciales
  ! campo para la dinamica de glauber
  do i=1, nneu
     hg(i)=0.d0

     sumah=0.d0
     do j=1, nneu
        sumah=sumah+2.*gch*ec(i,j)*s(j)
     end do
     hg(i)=sumah*s(i)

     ! write(*,*) h(i),th(i), i
  enddo

  ! campo para la dinamica de kawasaki

  do i=1,nneu
     do j=1, nneu
        hk(i,j)=0.d0
        sumah=0
        do ll=1,nneu
           sumah=sumah+gch*ec(i,ll)*s(ll)-gch*ec(j,ll)*s(ll)
        end do
        sumah=sumah-gch*ec(i,j)*s(j)+gch*ec(j,i)*s(i)
        sumah=sumah*(s(i)-s(j))
        hk(i,j)=sumah
     end do
  end do



  ! empezamos la simulacion

  niter=argumento(9)



  do i=1, nneu
     mean(i)=0.d0
  end do

  npinta=nneu



  do l=1,niter



     ir=1+dranu()*nneu
     jr=1+dranu()*nneu

     isigno=s(ir)*s(jr)

     ! decidimos que tipo de movimiento hacemos o spin-flip o kawasaki

     xr=dranu()


     pp=0.999

     if (xr.gt.pp) then
        delta=hg(ir)
        ss=1

     else if ((xr.lt.pp).and.(isigno.lt.0)) then
        delta=ec(ir,jr)*hk(ir,jr)
        ss=0
     end if




     xr=dranu()



     pphi=min(exp(-beta*delta),1.) ! rate de metropolis
     if (delta.eq.0) then
        pphi=0
     end if

     !pphi=0.5*(1.-tanh(beta*delta))

     if (pphi.gt.xr) then

        sss=s(ir)
        s(ir)=-ss*s(ir)+(1-ss)*s(jr)
        s(jr)=ss*s(jr)+(1-ss)*sss
     end if






     ! nuevos campos locales



     ! campo para la dinamica de glauber
     do i=1, nneu
        hg(i)=0.d0

        sumah=0.d0
        do j=1, nneu
           sumah=sumah+2.*gch*ec(i,j)*s(j)
        end do
        hg(i)=sumah*s(i)

        ! write(*,*) h(i),th(i), i
     enddo

     ! campo para la dinamica de kawasaki

     do i=1,nneu
        do j=1, nneu
           hk(i,j)=0.d0
           sumah=0.
           do ll=1,nneu
              sumah=sumah+gch*ec(i,ll)*s(ll)-gch*ec(j,ll)*s(ll)
           end do
           sumah=sumah-gch*ec(i,j)*s(j)+gch*ec(j,i)*s(i)
           sumah=sumah*(s(i)-s(j))
           hk(i,j)=sumah
        end do
     end do




     if (l>npinta) then

        do j=1,nneu
           if (s(j).gt.0) then
              write(10,*) j,l/nneu
           end if
        end do
        npinta=npinta+nneu
        flush(10)

        sumafiring=0
        do j=1,nneu
           sumafiring=sumafiring+s(j)/(0.+nneu)
        end do

        write(11,*) sumafiring, l/nneu
        flush(11)

     end if

  end do



end program isingchimera
