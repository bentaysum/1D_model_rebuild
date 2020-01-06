!***********************************************

!	param.par

!	Parameters for paramhr.f
!***********************************************

	integer    ninter
	parameter  (ninter=36)

	integer    nabs
	parameter  (nabs=13)

	integer    nz2
	parameter  (nz2=253)

	integer    ninter2
        parameter  (ninter2=16)

	real       kboltzman                  !cte Boltzman
	parameter  (kboltzman = 1.381e-16)
	
	real       n_avog                    !# de Avogadro
	parameter  (n_avog = 6.023e23)

	real       gg                        !cte gravitacion
	parameter  (gg = 6.67259e-8)

	real       masa                      !masa de Marte(g)
	parameter  (masa = 6.4163e26)

	real       radio
	parameter  (radio = 3390.)           !radio de Marte(km)

	integer	   nreact
	parameter  (nreact=93)




