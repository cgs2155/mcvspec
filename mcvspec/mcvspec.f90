!
!  MCVSPEC
!  ----
!
!  Nov 2021: Accretion column temperature/density profile calculator for magnetic CVs.
!
!  INPUT PARAMETERS:
!  PARAM(1) - B-field
!  PARAM(2) - Accretion Rate
!  PARAM(3) - WD mass in solar masses
!  PARAM(4) - Abundance
!  PARAM(5) - SIGMA_S -let this fit freely if
!  PARAM(6) - Alfven Radius -set to 1 if the system is an IP, it includes Alfven Radius in the calculation of x_s


SUBROUTINE MCVSPEC(EAR,NE,PARAM,IFL,PHOTAR,PHOTER)

  IMPLICIT NONE

  INTEGER    NE,IFL

  INTEGER    VGRIDMAX,RGRIDMAX

  PARAMETER  (VGRIDMAX=15000)
  PARAMETER  (RGRIDMAX=5)

  REAL     EAR(0:NE),PARAM(*),PHOTAR(NE),PHOTER(NE)
  REAL     METABUN

  REAL     M, B, B_approx,B_derived, SIGMA
  INTEGER  MSWITCH
  REAL     DISTNORM
  REAL     RHO(VGRIDMAX),P(VGRIDMAX),TK(VGRIDMAX),                &
       &               X(VGRIDMAX),NELEC(VGRIDMAX), SOLN(VGRIDMAX), &
       &   TAU(VGRIDMAX), PI_E(VGRIDMAX)

  INTEGER    J,VGRID,RGRID
  REAL, EXTERNAL :: func, func_prime


  REAL      Msun, Rsun, G, mue, mmw, A, mu, k, mH, alpha, beta, gamma
  REAL M_3 , M_wd, R_wd, shock_height, vff, n_elec_shock, Tshock, ES0, MDOT0, XS0, coeff, R_a, Bm, L
  REAL Lpart1, Lpart2

  REAL M_eff, C_log, m_e, pi, h, e, E0, twostreamcheck, X__, chi, p_a
  REAL mbar, zbar, zbarsqr, g_b
  common  /wdparams/ M_3, M_wd, R_wd, shock_height, vff, n_elec_shock, Tshock, ES0, MDOT0, XS0, coeff, R_a, Bm, L, p_a, SIGMA
  common /constants/ Msun, Rsun, G, mue, mmw, A, mu, k, mH, alpha, beta, gamma, X__, chi


  Msun     = 1.989100e+33 ! solar mass [gram]
  Rsun     = 6.9599e10 ! solar radius [cm]
  G        = 6.672590e-8 ! Gravitational constant
  mue      = 1.173 ! - mean number of electrons per baryons
  mmw      = 2.0    !- Mean molecular mass in white dwarf-
  mu       = 0.615  ! - mean molecular weight in accretion column (some papers use mu = 0.5 and Cropper 1999 did comparison)
  A        = 6.99e16 !used by Saxton
  k        = 1.380658e-16 ! -Boltzmann constant-
  mH       = 1.672623e-24 ! -Mass of hydrogen-
  alpha    = 2.0 ! -Thermodynamic Constant-
  beta     = 3.85 ! -Thermodynamic Constant-
  gamma    = 5.0/3.0  ! -Adidabatic Constant-
  !TSWITCH  = 0 ! parameter to choose whether or not to account for a tall accretion column and lower vff
  e        = 1.60218e-19 !charge of an electron
  E0       = 8.75419e-12 !permissivity of free space
  h        = 6.62607e-34 !planck's constant
  m_e      = 9.11e-28 !electron mass
  pi       = 3.141592654
  X__      = 2.7e34
  chi      = 1.9091 !constant that depends on the abundance weighted mean charge of ions
  mbar = 1.2886*mH
  zbar = 1.09987
  zbarsqr = 1.391
  g_b = 1.25 !gaunt factor

!
! RGRID is the number of radial zones in the accretion column.
!
  RGRID = 1
!
! These are the input parameters:
!  Mdot(0) is the specific accretion rate on the central axis of the col
!  Es(o) is the ratio of cyclotron to brems cooling at the shock front
!  Rc is the radius of the column in 10^9 cm
!  M is the mass of the primary in solar masses
!  METABUN is the metal abundance: solar = 1.0
!  VGRID is the number of vertical layers in the accretion column
!  ANG is the viewing angle used to modify the albedo in degrees
!
  B       = PARAM(1) ! Surface B-field [10^6 G]
  MDOT0   = PARAM(2) ! mass accretion rate [g/cm2/s]
  M       = PARAM(3) ! WD mass in solar mass
  METABUN = PARAM(4) ! Abundance
  SIGMA   = PARAM(5) ! calculate vff with Alfven radius
  MSWITCH = PARAM(6)
  VGRID   = 250 ! Number of vertical grids fixed to 50
!
! Flux normalization factor for APEC model as we define our flux norm as (R [km] / d [kpc])^2
! Note: PI*R_C^2/(4*PI*Distance^2) where R_C = accretion column radius [cm] and distance [cm]
!
  DISTNORM = 2.62511E-34

  DO J=1,NE
     PHOTER(J)=0.0
  ENDDO

! **This checks to make sure that there aren't too many vertical layers
! Checks that VGRID is not too large

  IF(VGRID.GT.VGRIDMAX)THEN
     PRINT *,'Too many vertical elements, max = ',VGRIDMAX
     VGRID = VGRIDMAX
  ENDIF

!           Calculating Neccesary WD Characteristics
  M_3  = (5.816*Msun)/(mmw**2.0) ! -Chandresekar Mass [gram]
  M_wd = M*Msun ! WD mass [grams]

  R_wd = Rsun*(0.0225/mmw)*SQRT(1.0-(M_wd/M_3)**(4.0/3.0))/((M_wd/M_3)**(1.0/3.0))  ! WD radius [cm]
  Bm   = (B*1e5)*((R_wd)**3) !Magnetic Moment

  Lpart1 = G*M_wd/R_wd
  Lpart2 = (8.6e17*(MDOT0**(7./5.))*(M**(1./7.))*(((R_wd/(1.e9))**(9./5.)))*(B**(-4./5.)))

  L    = Lpart1*Lpart2 !Accretion Luminosity

  R_a  = (2.75e8)*((M)**(1.0/7.0))*(((R_wd/(1.e9))**(-2.0/7.0))) *((L/1.0e33)**(-2.0/7.0))*((Bm/1.0e30)**(4.0/7.0)) !Alfven Radius
  PRINT *,"R_a = ", R_a

  IF (MSWITCH.EQ.1) THEN
    vff = ((2.*G*M_wd)*((1./(shock_height+R_wd))-(1./R_a)))**(1./2.) ! Chuck's calculation for vff
  ELSE
    vff = ((2.*G*M_wd)/(R_wd))**(1./2.) ! -Free Fall Velocity- at the shock height, Mukai 2017
  ENDIF

  p_a = MDOT0/vff

  n_elec_shock = 4.0*7.013e23*MDOT0/vff !electron number density at shock

  TShock = (3./16.)* (m_e*(vff**2)/k)*(1/(1/SIGMA+1))*((zbar+(mbar/m_e))/zbar)

  XS0 = vff**3*0.049/(2.*A*MDOT0) ! shock height [cm] when B = 0

  coeff = 9.1e-3*(B/10.)**2.85*(TShock/1.e8)**2.0*(n_elec_shock/1.e16)**(-1.85)*(XS0/1.e7)**(-0.85) !coefficient for ratio of power law cooling to bremmstraulung cooling

  ES0 =  ((vff/1e8)**4)*((p_a/1e-8)**(-1.85))* &
        & ((B/10)**2.85)*((1)**(-.425))*( 2.13e-16)* &
        & ((zbar+ (mbar/m_e))**3.85)/(g_b*(zbar**2.85)*zbarsqr*((1+1/SIGMA)**2))


  CALL SHOCK_APPROX( (3./(4. *( 1.+1./SIGMA))), VGRIDMAX, TAU, PI_E)
!           Checking WD characteristics calculation
  write(*,*) M_3, M_wd, R_wd, shock_height, vff, MDOT0
  write(*,*) "B [MG] = ", B, "  ES0 = ", ES0

  B_approx = 52.2*(ES0**0.35)*(TShock/1.0E8)**(-0.7)*(n_elec_shock/1.0E16)**0.65*(shock_height/1.E7)**0.3

  Call TWO_TEMPERATURE(VGRID,VGRIDMAX,RHO,P,TK,X,NELEC,SOLN, TAU, PI_E)


! Calculate the magnetic field. This is done using eqn(10) of Wu, Chanmu
! and Shaviv (1994) ApJ 426, 664, inverted to determine B from the other
! parameters. KM note: I think B is in unit of MegaGauss.
!


  B_derived = 52.2*ES0**0.35*                                          &
       &             (1.0E8/TK(VGRID))**0.7*                              &
       &             (NELEC(VGRID)/1.0E16)**0.65*                         &
       &             (X(VGRID)/1.0E7)**0.3

  PRINT '(A, F10.2)', 'White dwarf radius [10^7 cm] = ', R_wd*1e-7
  PRINT '(A, F10.2)', 'Shock height [10^7 cm] (numerical) = ', X(VGRID)*1e-7
  PRINT '(A, F10.2)', 'Shock height [10^7 cm] (approximate) = ', shock_height*1e-7
  PRINT '(A, F10.2)', 'Tshock [keV] = ', TShock*8.618e-8
  PRINT '(A, F10.2)', 'B [MG] (numerical) = ', B_derived
  PRINT '(A, F10.2)', 'B [MG] (approximate) = ', B_approx
  PRINT '(A, F10.2)', 'V_ff [10^8 cm/s] = ', vff*1e-8
  write(*,*) ES0, TK(VGRID), TSHOCK, NELEC(VGRID), X(VGRID), shock_height
  !
  IF(X(VGRID).GT.(R_wd*.8)) THEN
    PRINT *,'WARNING WARNING WARNING'
    PRINT *,'SHOCK HEIGHT IS COMPARABLE TO WHITE DWARF RADIUS'
    PRINT *,'1D ACCRETION COLUMN IS NO LONGER VALID'
  ENDIF

  C_log     = log(A*(RHO(VGRID)**2.) *((P(VGRID)/RHO(VGRID))**(1./2.))*(1.+ES0))
  M_eff     = (1-(MDOT0/(1.12e21*R_wd/1e9)))*M_wd !from Wu 1994

  twostreamcheck = (4.06e-2)*ES0*(1e9/R_wd)*(M_eff/Msun)*(15/C_log)

  PRINT *,'Ratio of electron-ion energy exchange timescale to cyclotron timescale = ', twostreamcheck
  IF(twostreamcheck.GT.1) THEN
    PRINT *,'WARNING WARNING WARNING'
    PRINT *,'COOLING IS NO LONGER CYCLOTRON DOMINATED'
    PRINT *,'ELECTRON THERMAL CONDUCTION AND COMPTON COOLING SHOULD BE CONSIDERED'
  ENDIF

  PRINT *,"Alfven Radius (10^7 cm) = ", R_a/(1e7)
  PRINT *,"Alfven Radius/White Dwarf Radius = ", R_a/R_wd

  CALL MCVSPEC_APEC(VGRID,X,TK,NELEC,METABUN,    &
       &              DISTNORM,EAR,NE,IFL,PHOTAR)

  RETURN
END SUBROUTINE MCVSPEC

SUBROUTINE MCVSPEC_APEC(VGRID,X,TK,NELEC,METABUN,    &
     &                        DISTNORM,EAR,NE,IFL,PHOTAR)

  IMPLICIT NONE

  INTEGER    IFL,NE
  INTEGER    VGRID,J,L
  REAL     EAR(0:NE),PHOTAR(NE)!,ANG!,FACT
  REAL     PARAM1(3)
  REAL     METABUN,FLX(NE),FLXERR(NE)
  REAL     TK(VGRID),X(VGRID),NELEC(VGRID)
  !Two Temperature Additions
  !REAL

  REAL     KK,PI,NENH
  REAL     DISTNORM!,DIST,NORM_FACTOR

!  REAL     GAMMA,THETA,ALPHA,ALB

  PI = 3.141592654
  !
  ! Constant to convert T(Kelvin) to kT
  !
  KK  = 11.6048E6
  !
  ! Constant to convert Ne to Nh
  !
  NENH = 1.21
  !
  ! First zero the flux array
  !
  DO L=1,NE
     PHOTAR(L) = 0.0
  ENDDO
  !
  ! Main loop to label 100.
  ! Loop over each vertical element to calculate and add the flux into the
  ! array
  !
  DO 100 J=1,VGRID
     !
     ! Calculates the Mewe spectrum for each vertical element on the energy
     ! grid passed into it, using
     !  a) PARAM(1) the temperature in keV of the vertical element
     !  b) PARAM(2) the density of electrons nH (cm^-3) of the vertical eleme
     !  c) PARAM(3) the heavy metal abundance
     !  d) PARAM(4) is the redshift, here set to zero
     !

     ! MEKAL parameter set
     !       PARAM1(1) = TK(J)/KK
     !       PARAM1(2) = NELEC(J)
     !       PARAM1(3) = METABUN
     !       PARAM1(4) = 0.0
     !       PARAM1(5) = 0.0
     !       PARAM1(6) = 0.0

     ! APEC parameter set
     PARAM1(1) = TK(J)/KK
     PARAM1(2) = METABUN
     PARAM1(3) = 0.0

     IF(TK(J)/KK.LT.86.)THEN
        CALL APEC(EAR,NE,PARAM1,IFL,FLX,FLXERR)
     ELSE
        CALL BREMSS(EAR,NE,PARAM1,IFL,FLX,FLXERR) ! added by KM if kT > 86 keV as it exceeds APEC's temperature upper limit.
     ENDIF

     ! Now multiplies the calculated spectrum by the volume and density^2
     ! of the particular element for each energy bin. This is required becaus
     ! the XSMEKL code has the density stripped out of the flux it returns,
     ! expecting this to be applied in the XSPEC normalisation. The 10^14 ari
     ! because of the units of 10^50 cm^3 at a distance of 1 pc = 3.086x10^18
     ! (ie 10^50/(10^(18^2))).
     !
     ! Explicitly we have:
     !  lumin(fmekal) = flux(fmekal)*(4pi*(3.086^18)^2)
     ! so for 1 cm^3 of gas
     !  lumin(fmekal) = flux(fmekal)*(4pi*(3.086^18)^2)/10^50
     ! as this passes to XSMEKL through XSVMKL, flux(fmekal) gets multiplied
     ! by 4pi*3.086^2/Ne^2 for XSPEC normalisations, so this needs to be reve
     ! so that to get the correct luminosity from XSMEKL we must have
     !  lumin(xsmekl) = flux(fmekal)*[(4pi*(3.086^18)^2)/10^50]*[Ne^2/(4pi*3.
     !                = flux(fmekal)*Ne^2*10^-14
     ! and then to get the fluxes from a volume element dV we then have
     !  flux(xsmekl)  = lumin(xsmekl)*dV/(4pi*dist^2)   (dist in cm)
     !                = flux(fmekal)*Ne^2*(10^-14)*dV/(4pi*dist^2)
     !
     ! On 1 June 1998 this routine was changed to reflect the fact that the
     ! XSVMKL normalisation is acutally 4pi*3.086^2/(NeNh)
     !                                = 4pi*3.086^2/(Ne^2/1.21)
     ! rather than 4pi*3.086^2/(Ne^2)
     !
     !DIST = (DISTNORM*3.085678E18)**2
     DO L=1,NE
        IF(J.EQ.1)THEN
            FLX(L) = DISTNORM*X(J)*((NELEC(J)**2)/NENH)*1.0E-14*FLX(L)
        ELSE
           FLX(L) = DISTNORM*(X(VGRID)-X(J-1))*((NELEC(J)**2)/NENH)*1.0E-14*FLX(L)
        ENDIF

     ENDDO
     !
     ! adds the flux for this vertical element to the existing flux in each
     ! energy bin
     !
     DO L=1,NE
        PHOTAR(L) = PHOTAR(L) + FLX(L)
     ENDDO

100 END DO

  RETURN
END SUBROUTINE MCVSPEC_APEC

!! TWO Temperature EQUATIONS
SUBROUTINE TWO_TEMPERATURE(VGRID,VGRIDMAX, RHO,P,TK,X,NELEC,SOLN, TAU_, PI_E_)
  !shock temp is temperature at shock which is the same from ions and electrons
  common  /wdparams/ M_3, M_wd, R_wd, shock_height, vff, n_elec_shock, Tshock, ES0, MDOT0, XS0, coeff, R_a, Bm, L, p_a, SIGMA
  common /constants/ Msun, Rsun, G, mue, mmw, A, mu, k, mH, alpha, beta, gamma, X__, chi

  INTEGER steps, r, VGRID, VGRIDMAX, JUMP

  REAL    soln(VGRIDMAX), taufinal(VGRIDMAX), piefinal(VGRIDMAX), MDOT0
  REAL    RHO(VGRID), P_i(VGRID),T_i(VGRID), X(VGRID), NELEC(VGRID), NION(VGRID), TikeV(VGRID)
  REAL    P(VGRID),TK(VGRID), TkeV(VGRID)
  REAL    TAU_(VGRIDMAX), PI_E_(VGRIDMAX)


  REAL xi, xi_init, n
  REAL tau_init, tau_f

  xi_init = 0
  xi = xi_init

  tau_init = TAU_(VGRIDMAX)
  tau_f    = .25

  n = (tau_f - tau_init)/float(VGRIDMAX)

  taufinal(1)  = TAU_(VGRIDMAX)
  piefinal(1)  = PI_E_(VGRIDMAX)
  soln(1)      = 0

  DO steps = 2,VGRIDMAX

     taufinal(steps) = TAU_(VGRIDMAX - steps+1)
     piefinal(steps) = PI_E_(VGRIDMAX - steps+1)
     xi   = xi+ n*(dxidtau(taufinal(steps), piefinal(steps)))
     soln(steps)     = xi ! -This is the 'non-normalized position grid-

     !print *,"STEP NUMBER ", steps
     !print *,"xi   = ", xi
     !print *,"tau  =", taufinal(steps)
     !print *,"pi_e =", piefinal(steps)
  ENDDO

  JUMP = VGRIDMAX/VGRID

  DO r = 1, VGRID
    X(r)         = soln(r*JUMP)*shock_height ! -This is the position grid-
    RHO(r)       = p_a/taufinal(r*JUMP) ! -This is the density grid- in terms of mass per cubic centimeter
    P(r)         = piefinal(r*JUMP) * vff !electron pressure
    P_i(r)       = (1-taufinal(r*JUMP) - piefinal(r*JUMP))*(vff)!ion pressure
    NELEC(r)     = RHO(r) * 7.01e23!electron number density
    NION(r)      = NELEC(r) * 1.173!ion number density
    TK(r)        = (16.*(taufinal(r*JUMP)*(1. - taufinal(r*JUMP))*(1.0/3.0))) * TShock ! -This is the actual temperature grid-
    T_i(r)       = P_i(r)/(NION(r)*1.380658e-16 )!ion temperature
    TkeV(r)      = (8.6173e-8)*TK(r)
    TikeV(r)     = (8.6173e-8)*T_i(r)
  ENDDO
END SUBROUTINE TWO_TEMPERATURE

SUBROUTINE SHOCK_APPROX(p, tnumber, t, pi)
  common  /wdparams/ M_3, M_wd, R_wd, shock_height, vff, n_elec_shock, Tshock, ES0, MDOT0, XS0, coeff, R_a, Bm, L, p_a, SIGMA
  common /constants/ Msun, Rsun, G, mue, mmw, A, mu, k, mH, alpha, beta, gamma, X__, chi

  REAL pi_e_integral, step_size, pi_e, int_coeff1, int_coeff2, p
  REAL rt_coeff1, rt_coeff2
  REAL limit, m_s
  INTEGER tnumber, steps
  REAL pi(tnumber), t(tnumber)
  REAL M_wd, Msun

  m_s = M_wd/Msun
  !this is a quadratic regression fit for what the lower limit of tau
  !should be such that the ODE does not encounter stiffness
  limit = 0.005696503 + (0.06507075*m_s) - 0.1599213*(m_s**2) + 0.1441434*(m_s**3) - 0.05448718*(m_s**4) + 0.006410256*(m_s**5)
  !"safety padding"
  limit = limit + .002

  step_size = (.25-limit)/tnumber
  t_init = .25
  pi_e = p
  pi(1) = pi_e
  t(1) = t_init
  pi_e_integral = 0.;


  DO steps = 1, tnumber
      t(steps+1) = t(steps)-step_size
      pi(steps+1) = pi(steps) - dpiedtau(t(steps),pi(steps))*(step_size)
      int_coeff1 = (gamma*(1.-t(steps)) - t(steps)) /(1. + ES0*f_(t(steps),pi(steps)))
      rt_coeff1 = SQRT((t(steps)**3.)/pi(steps))
      int_coeff2 = (gamma*(1.-t(steps+1)) - t(steps+1))/(1. + ES0*f_(t(steps+1),pi(steps+1)))
      rt_coeff2 = SQRT((t(steps+1)**3.)/pi(steps+1))

      !Trapezoidal Riemmann Sum for Calculating the Integral to Find Shock Height
      pi_e_integral = pi_e_integral + (step_size*(((rt_coeff1*int_coeff1) + (rt_coeff2*int_coeff2))/2.))
  ENDDO


  rhs_coeff = (gamma-1.)*p_a*A/(vff**2.)

  shock_height = (pi_e_integral/rhs_coeff)
END SUBROUTINE SHOCK_APPROX

!Differential Equation for Normalized Height with respect to Normalized Velocity
REAL FUNCTION dxidtau(t,p)
  common  /wdparams/ M_3, M_wd, R_wd, shock_height, vff, n_elec_shock, Tshock, ES0, MDOT0, XS0, coeff, R_a, Bm, L, p_a, SIGMA
  common /constants/ Msun, Rsun, G, mue, mmw, A, mu, k, mH, alpha, beta, gamma, X__, chi

  REAL lambda, t, p

  lambda    = (gamma - 1.)*shock_height*(p_a)*(vff**(-2.))*A*SQRT(p/(t**3.)) * (1. + (ES0*f_(t,p)))

  dxidtau = (gamma*(1-t) - t)/lambda

  RETURN
END FUNCTION dxidtau

!Differential Equation for Normalized Electron Pressure with respect to Normalized Velocity
!Used for Shock Height Approximation
REAL FUNCTION dpiedtau(t,p)
  common  /wdparams/ M_3, M_wd, R_wd, shock_height, vff, n_elec_shock, Tshock, ES0, MDOT0, XS0, coeff, R_a, Bm, L, p_a, SIGMA
  common /constants/ Msun, Rsun, G, mue, mmw, A, mu, k, mH, alpha, beta, gamma, X__, chi

  REAL outterm, midterm, gamlam1, gamlam2

  outterm = (gamma*(1.-t)-t)/t
  midterm = (gamma*p)/(gamma*(1.-t)-t)
  gamlam1 =  (X__/(A*(vff**(2.))))
  gamlam2 = (1.-t-(chi*p))/(t*(p**2)*(1.+ (ES0*f_(t,p))))

  dpiedtau = outterm*(1. - midterm - (gamlam1*gamlam2))

  RETURN
END FUNCTION dpiedtau

!Uses ES0 as a coefficient for describing efficiency of a secondary cooling process
!relative to thermal bremsstrahlung cooling as a function of pressure and velocity
REAL FUNCTION f_(t,p)
  common  /wdparams/ M_3, M_wd, R_wd, shock_height, vff, n_elec_shock, Tshock, ES0, MDOT0, XS0, coeff, R_a, Bm, L, p_a, SIGMA
  common /constants/ Msun, Rsun, G, mue, mmw, A, mu, k, mH, alpha, beta, gamma, X__, chi
  f_ = (4.**(alpha+beta))/(3.**alpha)*(((1.+sigma)/(sigma))**alpha)*(p**alpha)*(t**beta)
  RETURN
END FUNCTION f_
