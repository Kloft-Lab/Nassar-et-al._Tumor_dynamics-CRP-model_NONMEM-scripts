$SIZES PD=-100

$PROBLEM Paclitaxel AUC-driven drug-induced tumour decay

$INPUT  ID       ; Patient number
        OCC      ; Treatment cycle
        DAY      ; Day within cycle
        TIME     ; Time of tumour size measurment
        DV       ; Tumour size (log-transformed)
        BSL      ; Baseline tumour size
        EVID     ; Even tidentifier
        MDV      ; Missing dependent variable
        FLGIG    ; Flag for tumour size records 
        AUC      ; Paclitaxel AUC per cycle
        SimNo    ; Simulation number (idnetifier with respect to multiple imputation)

$DATA data1.csv IGNORE=I  
IGNORE(ID.EQ.16037)
IGNORE(FLGIG.EQ.1)    ; exclude non-tumour size records
IGNORE(TIME.GT.5040)

$SUBROUTINES ADVAN6 TOL=5

$MODEL COMP=(SD)              ;1 Sum of diameter
       COMP=(SIZE8)           ;2 Tumour size at week 8
       
$PK 
" FIRST
" COMMON/PRCOMG/IDUM1,IDUM2,IMAX,IDUM4,IDUM5
" INTEGER IDUM1,IDUM2,IMAX,IDUM4,IDUM5
" IMAX=10000000 

; Sum of diameters     
 GR     = THETA(1) * EXP(ETA(1))         ; Growth rate
 BETA   = THETA(2) * EXP(ETA(2))         ; Drug induced decay
 BASE   = EXP(LOG(BSL)+THETA(3)*ETA(3))  ; Baseline tumour size accroding to B2 method
 LAMBDA = THETA(4)*EXP(ETA(4))           ; Decline in drug effect in a cycle
 
 A_0(1) = BASE
 A_0(2) = BASE

$DES

; Sum of diameters
 EFF = BETA * AUC*EXP(-LAMBDA*TIME)
 DADT(1) = GR - EFF * A(1)
 CT = A(1)

;Output of week 8 size 
 IF(T.LT.1344) SLOP2 = GR - EFF * A(1)           ; 8*7*24=1344
 IF(T.GE.1344) SLOP2 = 0                         ; 8*7*24=1344
 DADT(2) = SLOP2
   
$ERROR 

;Output of amounts 
 AA1 = A(1)
 AA2 = A(2)

;Relatve change of size at week 8
TSIZE = AA2
RS8 = 100/BSL*TSIZE
                
;Residual error  
 IPRED=0.0001
 IF(A(1).GT.0) IPRED=LOG(A(1))
 W=THETA(3)
 Y=IPRED + W *EPS(1)
 IRES=DV-IPRED
 IWRES=IRES/W
 
$THETA
(0, 0.002)    ;1. GR 
(0, 0.0003)   ;2. BETA
(0, 0.12)     ;3. ERR
(0, 0.0004)   ;4. LAMBDA

$OMEGA 
 0.4       ;1. IIV GR
 0.6       ;2. IIV BETA
 1 FIX     ;3. IIV BSL, residual
 0.5       ;4. IIV LAMBDA, residual
 
$SIGMA 
 1 FIX     ;1.  RES                                              

$ESTIMATION PRINT=5 MAXEVAL=99999 METHOD=1 INTER NOABORT POSTHOC SIGDIG=3 MSFO=MSF001 NOTHETABOUNDTEST NOOMEGABOUNDTEST NOSIGMABOUNDTEST
$COVARIANCE PRINT=E
