  SUBROUTINE STOMATA (VEGTYP  ,MPE     ,APAR    ,FOLN    ,ILOC    , JLOC, & !in
       &              TV      ,EI      ,EA      ,SFCTMP  ,SFCPRS  , & !in
       &              O2      ,CO2     ,IGS     ,BTRAN   ,RB      , & !in
       &              RS      ,PSN)                                   !out
    USE NOAHMP_VEG_PARAMETERS
    IMPLICIT NONE
    ! --------------------------------------------------------------------------------------------------
    ! input
    INTEGER,INTENT(IN)  :: ILOC   !grid index
    INTEGER,INTENT(IN)  :: JLOC   !grid index
    INTEGER,INTENT(IN)  :: VEGTYP !vegetation physiology type
  
    REAL, INTENT(IN)    :: IGS    !growing season index (0=off, 1=on)
    REAL, INTENT(IN)    :: MPE    !prevents division by zero errors
  
    REAL, INTENT(IN)    :: TV     !foliage temperature (k)
    REAL, INTENT(IN)    :: EI     !vapor pressure inside leaf (sat vapor press at tv) (pa)
    REAL, INTENT(IN)    :: EA     !vapor pressure of canopy air (pa)
    REAL, INTENT(IN)    :: APAR   !par absorbed per unit lai (w/m2)
    REAL, INTENT(IN)    :: O2     !atmospheric o2 concentration (pa)
    REAL, INTENT(IN)    :: CO2    !atmospheric co2 concentration (pa)
    REAL, INTENT(IN)    :: SFCPRS !air pressure at reference height (pa)
    REAL, INTENT(IN)    :: SFCTMP !air temperature at reference height (k)
    REAL, INTENT(IN)    :: BTRAN  !soil water transpiration factor (0 to 1)
    REAL, INTENT(IN)    :: FOLN   !foliage nitrogen concentration (%)
    REAL, INTENT(IN)    :: RB     !boundary layer resistance (s/m)
  
    ! output
    REAL, INTENT(OUT)   :: RS     !leaf stomatal resistance (s/m)
    REAL, INTENT(OUT)   :: PSN    !foliage photosynthesis (umol co2 /m2/ s) [always +]
  
    ! in&out
    REAL                :: RLB    !boundary layer resistance (s m2 / umol)
    ! ---------------------------------------------------------------------------------------------
  
    ! ------------------------ local variables ----------------------------------------------------
    REAL :: CI     !internal co2 (pa)
    REAL, PARAMETER :: CIERR = 5e-2  !threshold of terminating the bisection (Pa)
    REAL :: CIHI, CILOW              !intermediate inner leaf CO2 pressure (Pa)
    REAL :: FCIHI, FCILOW, FCI       !intermediate inner leaf CO2 pressure (Pa)
    INTEGER :: ITER                  !iteration index
    INTEGER, PARAMETER :: NITER = 20 !number of iterations
  
    REAL :: TC          !foliage temperature (degree Celsius)
    REAL :: KC          !co2 Michaelis-Menten constant (pa)
    REAL :: KO          !o2 Michaelis-Menten constant (pa)
    REAL :: FNF         !foliage nitrogen adjustment factor (0 to 1)
    REAL :: PPF         !absorb photosynthetic photon flux (umol photons/m2/s)
    REAL :: CP          !co2 compensation point (pa)
    REAL :: AWC         !intermediate calculation for wc
    REAL :: VCMX        !maximum rate of carbonylation (umol co2/m2/s)
    REAL :: J           !electron transport (umol co2/m2/s)
    REAL :: CEA         !constrain ea or else model blows up
    REAL :: CF          !s m2/umol -> s/m
  
    REAL :: T
    ! ---------------------------------------------------------------------------------------------
  
    ! initialize RS=RSMAX and PSN=0 because will only do calculations
    ! for APAR > 0, in which case RS <= RSMAX and PSN >= 0
  
    CF = SFCPRS / (8.314 * SFCTMP) * 1.0e06
    RS = 1.0 / BP(VEGTYP) * CF
    PSN = 0.0
    CI = CO2
  
    IF (APAR <= 0.0) RETURN
  
    FNF = MIN( FOLN/MAX(MPE,FOLNMX(VEGTYP)), 1.0 )
    TC  = TV - TFRZ
    PPF = 4.6 * APAR
    J   = PPF * QE25(VEGTYP)
    ! F1 = AB ** ((BC - 25.0) / 10.0)
    ! KC  = KC25(VEGTYP) * F1(AKC(VEGTYP),TC)
    ! KO  = KO25(VEGTYP) * F1(AKO(VEGTYP),TC)
    KC  = KC25(VEGTYP) * AKC(VEGTYP) ** ((TC - 25.0) / 10.0)
    KO  = KO25(VEGTYP) * AKO(VEGTYP) ** ((TC - 25.0) / 10.0)
    AWC = KC * (1.0 + O2 / KO)
    CP  = 0.5 * KC / KO * O2 * 0.21
    ! VCMX = VCMX25(VEGTYP) / F2(TC) * FNF * BTRAN * F1(AVCMX(VEGTYP),TC)
    ! F2 = 1.0 + EXP((-2.2E05 + 710.0 * (AB + 273.16)) / (8.314 * (AB + 273.16)))
    VCMX = VCMX25(VEGTYP) &
       & / (1.0 + EXP((-2.2E05 + 710.0 * (TC + TFRZ)) / (8.314 * (TC + TFRZ)))) &
       & * FNF * BTRAN &
       & * (AVCMX(VEGTYP) ** ((TC - 25.0) / 10.0))
  
    ! rb: s/m -> s m**2 / umol
    RLB = RB/CF
  
    ! endpoints of the search intervals
    CIHI = 1.5 * CO2
    CILOW = 0.0
    DO ITER = 1, NITER
       CI = 0.5 * (CIHI + CILOW)
       CALL CI2CI(CI, FCI, RS, PSN)
       IF (((CIHI - CILOW) <= CIERR) .OR. ABS(FCI - CI) <= MPE) THEN
          EXIT
       ELSEIF (FCI > CI) THEN
          CILOW = CI
       ELSE
          CIHI = CI
       END IF
    END DO
  
    RS = RS*CF
  
  CONTAINS
    SUBROUTINE CI2CI(CI, FCI, RS, PSN)
      !function for serching the fixed point of CI, that is CI = FCI(CI)
      IMPLICIT NONE
      REAL, INTENT(IN) :: CI
      REAL, INTENT(OUT) :: FCI
      REAL, INTENT(OUT) :: RS
      REAL, INTENT(OUT) :: PSN
      REAL :: WC          !Rubisco limited photosynthesis (umol co2/m2/s)
      REAL :: WJ          !light limited photosynthesis (umol co2/m2/s)
      REAL :: WE          !export limited photosynthesis (umol co2/m2/s)
      REAL :: CS          !co2 concentration at leaf surface (pa)
      REAL :: A, B, C, Q  !intermediate calculations for RS
      REAL :: R1, R2      !roots for RS
  
      WJ = MAX(CI-CP,0.0) * J / (CI + 2.0 * CP)*C3PSN(VEGTYP) + J*(1.-C3PSN(VEGTYP))
      WC = MAX(CI-CP,0.0) * VCMX / (CI + AWC)*C3PSN(VEGTYP) + VCMX*(1.-C3PSN(VEGTYP))
      WE = 0.5 * VCMX * C3PSN(VEGTYP) + 4000.0 * VCMX * CI / SFCPRS*(1.-C3PSN(VEGTYP))
      PSN = MIN(WJ, WC, WE) * IGS
  
      CS = MAX(CO2 - 1.37 * RLB * SFCPRS * PSN, MPE )
      A = MP(VEGTYP) * PSN * SFCPRS * EA / (CS * EI) + BP(VEGTYP)
      B = (MP(VEGTYP) * PSN * SFCPRS / CS + BP(VEGTYP) ) * RLB - 1.
      C = -RLB
      IF (B >= 0.0) THEN
         Q = -0.5 * (B + SQRT(B * B - 4.0 * A * C))
      ELSE
         Q = -0.5 * (B - SQRT(B * B - 4.0 * A * C) )
      END IF
      R1 = Q / A
      R2 = C / Q
      RS = MAX(R1, R2)
  
      FCI = MAX(CS - PSN * SFCPRS * 1.65 * RS, 0.0)
    END SUBROUTINE CI2CI
  END SUBROUTINE STOMATA
