$PROBLEM    

$INPUT      ID    ; Patient number
            TIME  ; Time of event
            EVID  ; Even tidentifier
            DV    ; Progression or death/Censored event
            SIM   ; Flag for simulation records
            C3D1  ; CRP concentration at cycle 3
	    ABC3C2 ; Difference in CRP concnetration between cycle 3 and 2	
            LM   ; Flag for landmark time (day 42)

$DATA      5-CRP_PFS_COV_LM.csv IGNORE=@ 
IGNORE=(ID.EQ.13028)         ;Patient with no CRP (all are BLQ)
IGNORE=(LM.EQ.1)             ;Patients with events before cycle 3

;Sim_start : Add/remove for simulation
      IGNORE=(SIM.EQ.1)
;Sim_end

$ABBR  PROTECT
$SUBR  ADVAN =13 TOL = 9
$MODEL COMP  = (HAZARD)

$PK
;**************Covariates*************

;[1] CRP cycle 3
HZC3D1= ( 1 + THETA(3)*(C3D1))
;--------------------------------

;[2] Difference in CRP concentration between cycle 3 and cycle 2
HZABC3C2= ((ABC3C2)**THETA(4))
;--------------------------------

HZCOV=HZC3D1*HZABC3C2

;**************Hazard function: lognormal************* 

SDTTP=THETA(1)*EXP(ETA(1))  ; standard deviation of lognormal hazard model
MUTTP=THETA(2)              ; mean of lognormal hazard fucntion
PI   = 3.14159265

$DES
DEL   = 1E-16              ; Small number to avoid LOG(0)
TIM   = T+DEL
LTIM  = LOG(TIM)

X1X    = (LTIM-MUTTP)/SDTTP
PDF1X  = EXP(-(1/2)*((X1X)**2))/SQRT(2*PI)
LOGPFX = ((1/(TIM*SDTTP))*PDF1X/(1-PHI(X1X))) ; lognormal baseline hazard
DADT(1)= LOGPFX*HZCOV         

$ERROR
DELX = 1E-8 ; to avoid value zero of time

;----------TTE Model------------------------------
  IF(NEWIND.NE.2) OLDCHZ=0  ; Reset the cumulative hazard
  CHZ = A(1)-OLDCHZ         ; Cumulative hazard 
  OLDCHZ = A(1)             ; Rename old cumulative hazard
  SUR = EXP(-CHZ)           ; Survival probability

TIMX    = TIME+DELX
LTIMX   = LOG(TIMX)

X1     = (LTIMX-MUTTP)/SDTTP
PDF1   = EXP(-(1/2)*(X1**2))/SQRT(2*PI)
LOGPFS = ((1/(TIMX*SDTTP))*PDF1/(1-PHI(X1)))

HAZNOW = LOGPFS*HZCOV  ; hazard for PFS

  PDF = SUR*HAZNOW          ; Probability density function

  IF(DV.EQ.0) Y=SUR         ; Censored event (prob of not porgressing)
  IF(DV.NE.0) Y=SUR*HAZNOW  ; Probability density function of event

  IF(ICALL.EQ.4) THEN  ; For simulation
   CALL RANDOM (2,R)
   DV=0                ; Event censeored or not
   RTTE = 0            ; Flag to tell whether there is an event or not
   IF(TIME.EQ.20832) RTTE = 1 ; for the censored observation at 124 weeks
   IF(R.GT.SUR) THEN
      DV=1
      RTTE = 1
   ENDIF
  ENDIF

$THETA  (0,0.9)     ; 1. Standard deviation of lognormal hazard model
$THETA  (0,9)       ; 2. Mean of lognormal hazard fucntion
$THETA  (0,0.1)     ; 3. Covariate: CRP cycle 3
$THETA  (-0.26)     ; 4. Covariate: Difference in CRP concentration between cycle 3 and cycle 2

$OMEGA  0  FIX  ; 1. IIV standard deviation of lognormal hazard model

;Sim_start : add/remove for simulation
;$SIMULATION (5988566) (39978 UNIFORM) ONLYSIM NOPREDICTION SUB=100

$ESTIM MAXEVAL=9990 METHOD=0 LIKE PRINT=1 MSFO=msfb1 SIGL=9 NSIG=3
$COV PRINT=E
;Sim_end
