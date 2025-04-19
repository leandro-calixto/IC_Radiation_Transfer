program radiative_transfer
    use slw1_io
    use slw1_physics
    use radiation_solver
    implicit none

    ! Variáveis do domínio 3D
    integer, parameter :: nx = 50, ny = 50, nz = 50
    real, dimension(nx,ny,nz) :: T, p, x_h2o, x_co2
    real, dimension(3,nx,ny,nz) :: flux
    real, dimension(nx,ny,nz) :: source_term

    ! 1. Leitura de dados
    call read_input("input.dat", T, p, x_h2o, x_co2)

    ! 2. Cálculo das propriedades radiativas
    call calculate_radiative_properties(T, p, x_h2o, x_co2)

    ! 3. Solução da equação de transferência radiativa
    call solve_radiative_transfer(flux, source_term)

    ! 4. Saída dos resultados
    call write_output("flux.dat", flux, "source.dat", source_term)

    print *, "Simulação concluída com sucesso!"
end program radiative_transfer