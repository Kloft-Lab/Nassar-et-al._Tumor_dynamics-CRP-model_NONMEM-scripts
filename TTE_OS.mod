
$PROBLEM    

$INPUT      ID    ; Patient number
            TIME  ; Time of event
            EVID  ; Event identifier
            DV    ; Death/Censored event
            SIM   ; Flag for simulation records
            BLSD  ; Baseline tumor size
            C3D1  ; CRP concentration at cycle 3
	    ABC3C2 ; Difference in CRP concentration between cycle 3 and 2	
            LIVERLES  ; Liver lesions: Presence/absence
            LM   ; Flag for landmark time (day 42)

$DATA 4-CRP_SUR_COV_LM.csv IGNORE=@  
IGNORE=(ID.EQ.13028)         ;Patient with no CRP (all are BLQ)
IGNORE=(LM.EQ.1)             ;Patients with events before cycle 3

;Sim_start : Add/remove for simulation
      IGNORE=(SIM.EQ.1)
;Sim_end

$ABBR PROTECT

$SUBR ADVAN =13 TOL = 9

$MODEL COMP = (HAZARD)

$PK
;**************Covariates*************

;[1] CRP cycle 3
HZC3D1= ((C3D1)**THETA(3))
;--------------------------------

;[2] Liver lesions
IF(LIVERLES.EQ.0) HZLIVERLES = 1  ; Most common
IF(LIVERLES.EQ.1) HZLIVERLES = ( 1 + THETA(4))
;--------------------------------

;[3] Difference in CRP concentration between cycle 3 and cycle 2
HZABC3C2= ((ABC3C2)**THETA(5))
;--------------------------------

;[4] Baseline tumor size
HZBLSD= ((BLSD)**THETA(6))
;--------------------------------
HZCOV=HZC3D1*HZLIVERLES*HZABC3C2*HZBLSD


;**************Hazard function: Weibull************* 

HZLAM  = THETA(1)*EXP(ETA(1)) ; Scale factor
HZALPH = THETA(2)*EXP(ETA(2)) ; Shape factor

$DES
DEL = 1E-8 ; To avoid value zero of time

DADT(1) = (HZLAM*HZALPH*(HZLAM*(T+DEL))**(HZALPH-1))*HZCOV ; Weibull 

$ERROR
DELX = 1E-8 ; to avoid value zero of time

;----------TTE Model------------------------------
  IF(NEWIND.NE.2) OLDCHZ=0  ; Reset the cumulative hazard
  CHZ = A(1)-OLDCHZ         ; Cumulative hazard from previous time point in data set
  OLDCHZ = A(1)             ; Rename old cumulative hazard
  SUR = EXP(-CHZ)           ; Survival probability

 HAZNOW = (HZLAM*HZALPH*(HZLAM*(TIME+DELX))**(HZALPH-1))*HZCOV  ;weibull
                            ; Rate of event each time point

  PDF = SUR*HAZNOW          ; Probability density function

  IF(DV.EQ.0) Y=SUR         ; Censored event (prob of survival)
  IF(DV.NE.0) Y=SUR*HAZNOW  ; Probability density function of event

  IF(ICALL.EQ.4) THEN  ; For simulation
   CALL RANDOM (2,R)
   DV=0                ; Event censeored or not
   RTTE = 0            ; Flag to tell whether there is an event or not
   IF(TIME.EQ.23520) RTTE = 1 ; For the censored observation at 140 weeks
   IF(R.GT.SUR) THEN
      DV=1
      RTTE = 1
   ENDIF
  ENDIF

$THETA  (0.0000000001,0.00009)   ; 1. LAM: Scale factor
$THETA  (0, 1.1)                 ; 2. ALPHA: Shape factor
$THETA  (0,0.4)                  ; 3. Covariate: CRP cycle 3
$THETA  (0,0.5)                  ; 4. Covariate: Liver lesions
$THETA  (-0.3)                   ; 5. Covariate: Difference in CRP concentration between cycle 3 and cycle 2
$THETA  (0, 0.4)                 ; 6. Covariate: Baseline tumor size


$OMEGA  0  FIX       ; IIV LAM
$OMEGA  0  FIX       ; IIV SHP

;Sim_start : add/remove for simulation
;$SIMULATION (5988566) (39978 UNIFORM) ONLYSIM NOPREDICTION SUB=100

$ESTIM MAXEVAL=9990 METHOD=0 LIKE PRINT=1 MSFO=msfb1 SIGL=9 NSIG=3
$COV PRINT=E
;Sim_end


