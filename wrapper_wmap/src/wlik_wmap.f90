MODULE WMAP_EXTRA

	IMPLICIT NONE

	INTEGER:: BOK = 0
	INTEGER:: CLIK_LMAX,CLIK_LMIN
	real(8), dimension(:), allocatable :: cltt,clte,clee,clbb
	
END MODULE WMAP_EXTRA



SUBROUTINE WMAP_EXTRA_ONLY_ONE(MOK)
	USE WMAP_EXTRA
	INTEGER,INTENT(OUT)::MOK
	MOK = BOK
	BOK = 1
END SUBROUTINE 	WMAP_EXTRA_ONLY_ONE

SUBROUTINE WMAP_EXTRA_FREE()
	USE WMAP_EXTRA
	BOK =0
END SUBROUTINE 	WMAP_EXTRA_FREE

SUBROUTINE WMAP_EXTRA_LKL(LKL,CL)
	USE WMAP_EXTRA
	use wmap_likelihood_7yr
	use wmap_options
	use wmap_util
	use healpix_types
  
	REAL(8),INTENT(OUT)::LKL
	REAL(8) :: like(num_WMAP)
	REAL(8),INTENT(IN),DIMENSION(0:4*CLIK_LMAX+3)::CL
	INTEGER::i,cur

	!TT
	cur = 0
	cltt = 0
	clee = 0
	clbb = 0
	clte = 0
	
	DO i = ttmin,ttmax
		cltt(i)=CL(cur+i)*(i*(i+1.))/TWOPI
	END DO	
	cur = cur+clik_lmax+1


	!EE	
	DO i = temin,temax
		clee(i)=CL(cur+i)*(i*(i+1.))/TWOPI
	END DO	
	cur = cur+clik_lmax+1

	!BB
	DO i = temin,temax
		clbb(i)=CL(cur+i)*(i*(i+1.))/TWOPI
	END DO	
	cur = cur+clik_lmax+1
	
	!TE
	DO i = temin,temax
		clte(i)=CL(cur+i)*(i*(i+1.))/TWOPI
	END DO	
	cur = cur+clik_lmax+1
	
	CALL wmap_likelihood_compute(cltt,clte,clee,clbb,like)	
	
	LKL = -sum(like(1:num_WMAP))
END SUBROUTINE 	WMAP_EXTRA_LKL

SUBROUTINE WMAP_EXTRA_PARAMETER_INIT(tt_min,tt_max,te_min,te_max,m_use_gibbs,m_use_lowl_pol)
	USE WMAP_EXTRA
	use wmap_likelihood_7yr
	use wmap_options
	use wmap_util
	INTEGER,INTENT(IN)::tt_min,tt_max,te_min,te_max,m_use_gibbs,m_use_lowl_pol
	
	
	ttmin = tt_min
	ttmax = tt_max
	temin = te_min
	temax = te_max
	lowl_max = 30
	use_gibbs = .false.
	
	if (temax>800) then
		temax = 800
	endif
	if (ttmax>1200) then
		ttmax = 1200
	endif

	if (m_use_gibbs.EQ.1) then
		lowl_max = 32
		use_gibbs = .true.
	endif
	
	if (ttmin>lowl_max) then
		use_lowl_TT = .false.
	endif
	
	use_lowl_pol = .false.
	if (m_use_lowl_pol==1) then
		use_lowl_pol = .true.
	endif
	
	if (temin.ge.temax) then
		use_TE = .false.
	endif

	if (ttmin.ge.ttmax) then
		use_TT = .false.
	endif
	
	
	clik_lmax = tt_max
	if (te_max>clik_lmax) then
		clik_lmax = te_max
	endif
	clik_lmin = tt_min
	if (te_max<clik_lmin) then
		clik_lmin = te_min
	endif
	
	allocate( cltt(2:clik_lmax) )
	allocate( clte(2:clik_lmax) )
	allocate( clee(2:clik_lmax) )
	allocate( clbb(2:clik_lmax) )

	
	CALL wmap_likelihood_init()

END SUBROUTINE 	WMAP_EXTRA_PARAMETER_INIT
