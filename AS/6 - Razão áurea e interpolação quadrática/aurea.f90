program aurea

implicit none

real(8) ra_opt, p_max, ra1, ra2
real(8) r1, r2, r3
real(8) v
real(8) razao
real(8) ra_lower, ra_upper, length_1
real(8) erro, tol
integer(4) iter

! Para melhor compreensão, p(ra) = f(x)

iter = 0

tol = 1.0d-4
erro = 1.0d0

r1 = 8.0
r2 = 12.0
r3 = 10.0
v = 80.0

razao = (sqrt(5.0)-1.0)/2.0

ra_upper = 100
ra_lower = 0

write(*,*) "--------------------------------------------------------------------------------"
write(*,*) "                      OTIMIZACAO: METODO DA RAZAO AUREA                         "
write(*,*) "--------------------------------------------------------------------------------"
write(*,"(A)") " Iter |    ra_lower    |    ra_upper    |   ra_otimizado |      Erro      "
write(*,*) "--------------------------------------------------------------------------------"

do while (abs(erro).gt.tol)

iter = iter+1

length_1 = (ra_upper-ra_lower)*razao
ra1 = ra_lower+length_1
ra2 = ra_upper-length_1

    if (p(ra2).gt.p(ra1)) then
        ra_upper=ra1
        else
        ra_lower=ra2
    end if    

erro = ra1-ra2

ra_opt=(ra_upper+ra_lower)/2.0

write(*,*) iter, ra1, ra2, ra_opt, erro 

end do

ra_opt=(ra_upper+ra_lower)/2.0

p_max = p(ra_opt)
write(*,*) "-----------------------------------------------------------------------------------"
write(*,*) "Ra otimizado:", ra_opt
write(*,*) "Potência máxima:", p_max
write(*,*) "-----------------------------------------------------------------------------------"

contains
    
    function p(ra)
    real(8) p, ra, num, den
    num = v*r3*ra
    den = r1*(ra+r2+r3)+r3*ra+r3*r2
    p = ((num/den)**2)/ra
    end function
    
end program aurea
