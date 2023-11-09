
$SIZES     PD=-70
$PROBLEM   Coupled tumor size-CRP turnover model

$INPUT  ID       ; Patient number	
        OCC      ; Treatment cycle	
        DAY      ; Day within cycle
        TIME     ; Time of CRP measurement		
        EVID     ; Event identifier	
        MDV      ; Missing dependent variable	
        DV	 ; CRP (log-transformed)
        CMT	 ; Compartment identifier (CRP, Tumour)
        FLAGCRP  ; FLag for BLQ and ULOQ value for CRP
        BLIL6	 ; Baseline interleukin-6
        NEWBLIL6 ; Baseline interleukin-6, missing value imputed by median		
        SMOK	 ; Smoking status		
        STAGE	 ; Disease stage					
        BLSD     ; Baseline tumor size
        IGR	 ; Individual tumour growth rate
        IBETA	 ; Individual drug effect
        ITBASE	 ; Individual estimated baseline tumor size
        ILAMBDA	 ; Individual drug resistance rate constant
        SEGR	 ; Individual standard error of tumor growth rate
        SEBETA	 ; Individual standard error of drug effect
        SELAMBDA ; Individual standard error of resistance
        SEERR	 ; Individual standard error of baseline tumour size
        MAUC     ; Paclitaxel AUC per cycle

$DATA  2-CRP_TGI_Model.csv IGNORE=@

      IGNORE=(FLAGCRP.EQ.2)        ;Exclude flagged items (BLQ)
      IGNORE=(ID.EQ.13028)         ;Patient with no CRP (all are BLQ)

$ABBREV PROTECT
$SUBROUTINE   ADVAN13 TOL9  

$MODEL      NCOMP=3
            COMP=(TS)              ;1 Tumour size 
            COMP=(CRP)             ;2 CRP
            COMP=(SIZE8)           ;3 Tumour size at week 8
            
$PK
;-----------------------------------------------------------
;------ CRP model parameters -------------------------------
;-----------------------------------------------------------

;Covariate function: sum of diameters-exponential
KINSUMDIA = EXP(THETA(5)*(BLSD-8.25))

;Covariate function: IL6-linear
KINIL6 = ( 1 + THETA(4)*(NEWBLIL6-2.57))

;Covariate function: disease stage
IF(STAGE.EQ.1)  KINSTAGE = 1 ; Most common
IF(STAGE.EQ.0)  KINSTAGE = ( 1 + THETA(6))

;Covariate function: smoking status
IF(SMOK.EQ.1)  KINSMOK = 1 ; Most common
IF(SMOK.EQ.2)  KINSMOK = ( 1 + THETA(7))
IF(SMOK.EQ.3)  KINSMOK = ( 1 + THETA(8))

; Production constant: Kin, Zero-order

TVKIN = THETA(1)*KINIL6*KINSUMDIA*KINSTAGE*KINSMOK
KIN   = TVKIN*EXP(ETA(1))	

; Basal production constant
KINB =THETA(9)

KINT=KINB+KIN

; Degradation constant: Kout, 1st-order
TVKOUT = THETA(2)
KOUT   = TVKOUT*EXP(ETA(2)) 		

; Baseline
BASE = KINT/KOUT

; Relationship between tumor size and Kin 			              
TVSLP =THETA(10)
SLP   =TVSLP*EXP(ETA(7))

; Tumor growth inhibition model parameters
GR     = IGR*EXP(ETA(3)*SEGR)
BETA   = IBETA*EXP(ETA(4)*SEBETA)
LAMBDA = ILAMBDA*EXP(ETA(5)*SELAMBDA)
TBASE  = BLSD*EXP(ETA(6)*SEERR)

; Intialize compartments
A_0(1)       = TBASE   ; Baseline tumor size
A_0(2)       = BASE    ; CRP (mg/L)
A_0(3)       = TBASE   ; Baseline tumor size

$DES

; TGI model
EFF = BETA * MAUC*EXP(-LAMBDA*T)
DADT(1) = GR - EFF * A(1)

; CRP
KINP=KIN*SLP*(A(1)/TBASE)
DADT(2) = (KINB+KINP)-KOUT*A(2)		; Time-course of CRP

;Output of tumor size week 8 
 IF(T.LT.1344) SLOP2 = GR - EFF * A(1)           ; 8*7*24=1344
 IF(T.GE.1344) SLOP2 = 0                         ; 8*7*24=1344
 DADT(3) = SLOP2

$ERROR

;Output of amounts
TSIZE = A(1)

CRP = A(2)
IPRED=0
IF(CRP.NE.0) IPRED=LOG(CRP)

RWSIZE= A(3)

;Relatve change of size at week 8
RS8 = (RWSIZE/BLSD)*100

; CRP RUV model: additive
W    = THETA(3)                         ; Sigma fixed to 1
IRES = DV-IPRED
IWRES= IRES/W
Y    = IPRED+W*EPS(1)

$THETA
    (0,0.4)                ;1. KIN
    0.036473684 FIX        ;2. KOUT
    (0,0.8)                ;3. SIGMA SD
    (0.3)                  ;4. COV IL6 on Kin
    0 FIX                  ;5. COV Sumdia on Kin
    (-0.4)                 ;6. Stage 0: IIIB
    (0.6)                  ;7. Smok 2: former
    (1.3)                  ;8. Smok 3: current
    0.01094211 FIX         ;9. KINB set to equivalent to LLOQ
    (0.8)                  ;10. TS:linear

$OMEGA  0.8                  ;1. IIV-KIN
$OMEGA  0 FIX                ;2. IIV-KOUT
$OMEGA  1 FIX                ;3. IGR
$OMEGA  1 FIX                ;4. IBETA
$OMEGA  1 FIX                ;5. ILAMBDA
$OMEGA  1 FIX                ;6. IIV BSL, residual
$OMEGA  0.36                 ;7. IIV slope

$SIGMA  1 FIX                ; CRP ERROR

$ESTIMATION METHOD=1 INTER MAX=9999 PRINT=1 NOABORT

$COV






