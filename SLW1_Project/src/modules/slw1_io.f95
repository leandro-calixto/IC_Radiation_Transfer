module slw1_io
    implicit none
    
    contains
    ! Rotina para ler arquivo de entrada
    subroutine read_input(filename, T, p, x_h2o, x_co2, nx, ny, nz, success)
      character(len=*), intent(in) :: filename
      real, dimension(:,:,:), allocatable, intent(out) :: T, p, x_h2o, x_co2
      integer, intent(out) :: nx, ny, nz
      logical, intent(out) :: success
      integer :: i, j, k, ierr
      
      ! Inicialização
      success = .false.
      
      ! Abre o arquivo
      open(unit=10, file=filename, status='old', action='read', iostat=ierr)
      if (ierr /= 0) then
        print *, "Erro ao abrir arquivo ", filename
        return
      end if
      
      ! Lê dimensões do domínio
      read(10, *, iostat=ierr) nx, ny, nz
      if (ierr /= 0) then
        print *, "Erro ao ler dimensões"
        close(10)
        return
      end if
      
      ! Aloca arrays
      allocate(T(nx,ny,nz), p(nx,ny,nz), x_h2o(nx,ny,nz), x_co2(nx,ny,nz), stat=ierr)
      if (ierr /= 0) then
        print *, "Erro ao alocar memória"
        close(10)
        return
      end if
      
      ! Lê dados (formato simples: valores em ordem i,j,k)
      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            read(10, *, iostat=ierr) T(i,j,k), p(i,j,k), x_h2o(i,j,k), x_co2(i,j,k)
            if (ierr /= 0) then
              print *, "Erro ao ler dados na posição ", i, j, k
              close(10)
              return
            end if
          end do
        end do
      end do
      
      close(10)
      success = .true.
      print *, "Arquivo lido com sucesso. Dimensões:", nx, "x", ny, "x", nz
    end subroutine read_input
    
    ! Rotina para escrever resultados
    subroutine write_results(filename, flux, source)
      character(len=*), intent(in) :: filename
      real, dimension(:,:,:,:), intent(in) :: flux  ! flux(3,nx,ny,nz)
      real, dimension(:,:,:), intent(in) :: source
      integer :: i, j, k
      
      open(unit=20, file=filename, action='write')
      
      write(20, '(A)') "Resultados SLW1 - Formato: i j k flux_x flux_y flux_z source"
      
      do k = 1, size(source,3)
        do j = 1, size(source,2)
          do i = 1, size(source,1)
            write(20, '(3I5,4E15.6)') i, j, k, flux(1,i,j,k), flux(2,i,j,k), flux(3,i,j,k), source(i,j,k)
          end do
        end do
      end do
      
      close(20)
      print *, "Resultados escritos em ", filename
    end subroutine write_results
  end module slw1_io