      SUBROUTINE SSPACE (A,B,MAXA,R,EIGV,TT,W,AR,BR,VEC,D,RTOLV,BUP,BLO,SSP00001
     1 BUPC,NN,NNM,NWK,NWM,NROOT,RTOL,NC,NNC,NITEM,IFSS,IFPR,NSTIF,IOUT)SSP00002
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . SSP00003
C .                                                                   . SSP00004
C .   P R O G R A M                                                   . SSP00005
C .        TO SOLVE FOR THE SMALLEST EIGENVALUES-- ASSUMED .GT. 0 --  . SSP00006
C .        AND CORRESPONDING EIGENVECTORS IN THE GENERALIZED          . SSP00007
C .        EIGENPROBLEM USING THE SUBSPACE ITERATION METHOD           . SSP00008
C .                                                                   . SSP00009
C .  - - INPUT VARIABLES - -                                          . SSP00010
C .        A(NWK)    = STIFFNESS MATRIX IN COMPACTED FORM (ASSUMED    . SSP00011
C .                    POSITIVE DEFINITE)                             . SSP00012
C .        B(NWM)    = MASS MATRIX IN COMPACTED FORM                  . SSP00013
C .        MAXA(NNM) = VECTOR CONTAINING ADDRESSES OF DIAGONAL        . SSP00014
C .                    ELEMENTS OF STIFFNESS MATRIX A                 . SSP00015
C .        R(NN,NC)  = STORAGE FOR EIGENVECTORS                       . SSP00016
C .        EIGV(NC)  = STORAGE FOR EIGENVALUES                        . SSP00017
C .        TT(NN)    = WORKING VECTOR                                 . SSP00018
C .        W(NN)     = WORKING VECTOR                                 . SSP00019
C .        AR(NNC)   = WORKING MATRIX STORING PROJECTION OF K         . SSP00020
C .        BR(NNC)   = WORKING MATRIX STORING PROJECTION OF M         . SSP00021
C .        VEC(NC,NC)= WORKING MATRIX                                 . SSP00022
C .        D(NC)     = WORKING VECTOR                                 . SSP00023
C .        RTOLV(NC) = WORKING VECTOR                                 . SSP00024
C .        BUP(NC)   = WORKING VECTOR                                 . SSP00025
C .        BLO(NC)   = WORKING VECTOR                                 . SSP00026
C .        BUPC(NC)  = WORKING VECTOR                                 . SSP00027
C .        NN        = ORDER OF STIFFNESS AND MASS MATRICES           . SSP00028
C .        NNM       = NN + 1                                         . SSP00029
C .        NWK       = NUMBER OF ELEMENTS BELOW SKYLINE OF            . SSP00030
C .                    STIFFNESS MATRIX                               . SSP00031
C .        NWM       = NUMBER OF ELEMENTS BELOW SKYLINE OF            . SSP00032
C .                    MASS MATRIX                                    . SSP00033
C .                      I. E. NWM=NWK FOR CONSISTENT MASS MATRIX     . SSP00034
C .                            NWM=NN  FOR LUMPED MASS MATRIX         . SSP00035
C .        NROOT     = NUMBER OF REQUIRED EIGENVALUES AND EIGENVECTORS. SSP00036
C .        RTOL      = CONVERGENCE TOLERANCE ON EIGENVALUES           . SSP00037
C .                    ( 1.E-06 OR SMALLER )                          . SSP00038
C .        NC        = NUMBER OF ITERATION VECTORS USED               . SSP00039
C .                    (USUALLY SET TO MIN(2*NROOT, NROOT+8), BUT NC  . SSP00040
C .                    CANNOT BE LARGER THAN THE NUMBER OF MASS       . SSP00041
C .                    DEGREES OF FREEDOM)                            . SSP00042
C .        NNC       = NC*(NC+1)/2 DIMENSION OF STORAGE VECTORS AR,BR . SSP00043
C .        NITEM     = MAXIMUM NUMBER OF SUBSPACE ITERATIONS PERMITTED. SSP00044
C .                    (USUALLY SET TO 16)                            . SSP00045
C .                    THE PARAMETERS NC AND/OR NITEM MUST BE         . SSP00046
C .                    INCREASED IF A SOLUTION HAS NOT CONVERGED      . SSP00047
C .        IFSS      = FLAG FOR STURM SEQUENCE CHECK                  . SSP00048
C .                      EQ.0  NO CHECK                               . SSP00049
C .                      EQ.1  CHECK                                  . SSP00050
C .        IFPR      = FLAG FOR PRINTING DURING ITERATION             . SSP00051
C .                      EQ.0  NO PRINTING                            . SSP00052
C .                      EQ.1  PRINT                                  . SSP00053
C .        NSTIF     = SCRATCH FILE                                   . SSP00054
C .        IOUT      = UNIT USED FOR OUTPUT                           . SSP00055
C .                                                                   . SSP00056
C .  - - OUTPUT - -                                                   . SSP00057
C .        EIGV(NROOT) = EIGENVALUES                                  . SSP00058
C .        R(NN,NROOT) = EIGENVECTORS                                 . SSP00059
C .                                                                   . SSP00060
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . SSP00061
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               SSP00062
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . SSP00063
C .   THIS PROGRAM IS USED IN SINGLE PRECISION ARITHMETIC ON CRAY     . SSP00064
C .   EQUIPMENT AND DOUBLE PRECISION ARITHMETIC ON IBM MACHINES,      . SSP00065
C .   ENGINEERING WORKSTATIONS AND PCS. DEACTIVATE ABOVE LINE FOR     . SSP00066
C .   SINGLE PRECISION ARITHMETIC.                                    . SSP00067
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . SSP00068
      INTEGER MAXA(NNM)                                                 SSP00069
      DIMENSION A(NWK),B(NWM),R(NN,NC),TT(NN),W(NN),EIGV(NC),           SSP00070
     1          D(NC),VEC(NC,NC),AR(NNC),BR(NNC),RTOLV(NC),BUP(NC),     SSP00071
     2          BLO(NC),BUPC(NC)                                        SSP00072
C                                                                       SSP00073
C     SET TOLERANCE FOR JACOBI ITERATION                                SSP00074
      TOLJ=1.0D-12                                                      SSP00075
C                                                                       SSP00076
C     INITIALIZATION                                                    SSP00077
C                                                                       SSP00078
      ICONV=0                                                           SSP00079
      NSCH=0                                                            SSP00080
      NSMAX=12                                                          SSP00081
      N1=NC + 1                                                         SSP00082
      NC1=NC - 1                                                        SSP00083
      REWIND NSTIF                                                      SSP00084
      WRITE (NSTIF) A                                                   SSP00085
      DO 2 I=1,NC                                                       SSP00086
    2 D(I)=0.                                                           SSP00087
C                                                                       SSP00088
C     ESTABLISH STARTING ITERATION VECTORS                              SSP00089
C                                                                       SSP00090
      ND=NN/NC                                                          SSP00091
      IF (NWM.GT.NN) GO TO 4                                            SSP00092
      J=0                                                               SSP00093
      DO 6 I=1,NN                                                       SSP00094
      II=MAXA(I)                                                        SSP00095
      R(I,1)=B(I)                                                       SSP00096
      IF (B(I).GT.0) J=J + 1                                            SSP00097
    6 W(I)=B(I)/A(II)                                                   SSP00098
      IF (NC.LE.J) GO TO 16                                             SSP00099
      WRITE (IOUT,1007)                                                 SSP00100
      GO TO 800                                                         SSP00101
    4 DO 10 I=1,NN                                                      SSP00102
      II=MAXA(I)                                                        SSP00103
      R(I,1)=B(II)                                                      SSP00104
   10 W(I)=B(II)/A(II)                                                  SSP00105
   16 DO 20 J=2,NC                                                      SSP00106
      DO 20 I=1,NN                                                      SSP00107
   20 R(I,J)=0.                                                         SSP00108
C                                                                       SSP00109
      L=NN - ND                                                         SSP00110
      DO 30 J=2,NC                                                      SSP00111
      RT=0.                                                             SSP00112
      DO 40 I=1,L                                                       SSP00113
      IF (W(I).LT.RT) GO TO 40                                          SSP00114
      RT=W(I)                                                           SSP00115
      IJ=I                                                              SSP00116
   40 CONTINUE                                                          SSP00117
      DO 50 I=L,NN                                                      SSP00118
      IF (W(I).LE.RT) GO TO 50                                          SSP00119
      RT=W(I)                                                           SSP00120
      IJ=I                                                              SSP00121
   50 CONTINUE                                                          SSP00122
      TT(J)=FLOAT(IJ)                                                   SSP00123
      W(IJ)=0.                                                          SSP00124
      L=L - ND                                                          SSP00125
   30 R(IJ,J)=1.                                                        SSP00126
C                                                                       SSP00127
      WRITE (IOUT,1008)                                                 SSP00128
      WRITE (IOUT,1002) (TT(J),J=2,NC)                                  SSP00129
C                                                                       SSP00130
C     A RANDOM VECTOR IS ADDED TO THE LAST VECTOR                       SSP00131
C                                                                       SSP00132
      PI=3.141592654D0                                                  SSP00133
      XX=0.5D0                                                          SSP00134
      DO 60 K=1,NN                                                      SSP00135
      XX=(PI + XX)**5                                                   SSP00136
      IX=INT(XX)                                                        SSP00137
      XX=XX - FLOAT(IX)                                                 SSP00138
   60 R(K,NC)=R(K,NC) + XX                                              SSP00139
C                                                                       SSP00140
C     FACTORIZE MATRIX A INTO (L)*(D)*(L(T))                            SSP00141
C                                                                       SSP00142
      ISH=0                                                             SSP00143
      CALL DECOMP (A,MAXA,NN,ISH,IOUT)                                  SSP00144
C                                                                       SSP00145
C - - - S T A R T   O F   I T E R A T I O N   L O O P                   SSP00146
C                                                                       SSP00147
      NITE=0                                                            SSP00148
      TOLJ2=1.0D-24                                                     SSP00149
  100 NITE=NITE + 1                                                     SSP00150
      IF (IFPR.EQ.0) GO TO 90                                           SSP00151
      WRITE (IOUT,1010) NITE                                            SSP00152
C                                                                       SSP00153
C     CALCULATE THE PROJECTIONS OF A AND B                              SSP00154
C                                                                       SSP00155
   90 IJ=0                                                              SSP00156
      DO 110 J=1,NC                                                     SSP00157
      DO 120 K=1,NN                                                     SSP00158
  120 TT(K)=R(K,J)                                                      SSP00159
      CALL REDBAK (A,TT,MAXA,NN)                                        SSP00160
      DO 130 I=J,NC                                                     SSP00161
      ART=0.                                                            SSP00162
      DO 140 K=1,NN                                                     SSP00163
  140 ART=ART + R(K,I)*TT(K)                                            SSP00164
      IJ=IJ + 1                                                         SSP00165
  130 AR(IJ)=ART                                                        SSP00166
      DO 150 K=1,NN                                                     SSP00167
  150 R(K,J)=TT(K)                                                      SSP00168
  110 CONTINUE                                                          SSP00169
      IJ=0                                                              SSP00170
      DO 160 J=1,NC                                                     SSP00171
      CALL MULT (TT,B,R(1,J),MAXA,NN,NWM)                               SSP00172
      DO 180 I=J,NC                                                     SSP00173
      BRT=0.                                                            SSP00174
      DO 190 K=1,NN                                                     SSP00175
  190 BRT=BRT + R(K,I)*TT(K)                                            SSP00176
      IJ=IJ + 1                                                         SSP00177
  180 BR(IJ)=BRT                                                        SSP00178
      IF (ICONV.GT.0) GO TO 160                                         SSP00179
      DO 200 K=1,NN                                                     SSP00180
  200 R(K,J)=TT(K)                                                      SSP00181
  160 CONTINUE                                                          SSP00182
C                                                                       SSP00183
C     SOLVE FOR EIGENSYSTEM OF SUBSPACE OPERATORS                       SSP00184
C                                                                       SSP00185
      IF (IFPR.EQ.0) GO TO 320                                          SSP00186
      IND=1                                                             SSP00187
  210 WRITE (IOUT,1020)                                                 SSP00188
      II=1                                                              SSP00189
      DO 300 I=1,NC                                                     SSP00190
      ITEMP=II + NC - I                                                 SSP00191
      WRITE (IOUT,1005) (AR(J),J=II,ITEMP)                              SSP00192
  300 II=II + N1 - I                                                    SSP00193
      WRITE (IOUT,1030)                                                 SSP00194
      II=1                                                              SSP00195
      DO 310 I=1,NC                                                     SSP00196
      ITEMP=II + NC - I                                                 SSP00197
      WRITE (IOUT,1005) (BR(J),J=II,ITEMP)                              SSP00198
  310 II=II + N1 - I                                                    SSP00199
      IF (IND.EQ.2) GO TO 350                                           SSP00200
C                                                                       SSP00201
  320 CALL JACOBI (AR,BR,VEC,EIGV,W,NC,NNC,TOLJ,NSMAX,IFPR,IOUT)        SSP00202
C                                                                       SSP00203
      IF (IFPR.EQ.0) GO TO 350                                          SSP00204
      WRITE (IOUT,1040)                                                 SSP00205
      IND=2                                                             SSP00206
      GO TO 210                                                         SSP00207
C                                                                       SSP00208
C     ARRANGE EIGENVALUES IN ASCENDING ORDER                            SSP00209
C                                                                       SSP00210
  350 IS=0                                                              SSP00211
      II=1                                                              SSP00212
      DO 360 I=1,NC1                                                    SSP00213
      ITEMP=II + N1 - I                                                 SSP00214
      IF (EIGV(I+1).GE.EIGV(I)) GO TO 360                               SSP00215
      IS=IS + 1                                                         SSP00216
      EIGVT=EIGV(I+1)                                                   SSP00217
      EIGV(I+1)=EIGV(I)                                                 SSP00218
      EIGV(I)=EIGVT                                                     SSP00219
      BT=BR(ITEMP)                                                      SSP00220
      BR(ITEMP)=BR(II)                                                  SSP00221
      BR(II)=BT                                                         SSP00222
      DO 370 K=1,NC                                                     SSP00223
      RT=VEC(K,I+1)                                                     SSP00224
      VEC(K,I+1)=VEC(K,I)                                               SSP00225
  370 VEC(K,I)=RT                                                       SSP00226
  360 II=ITEMP                                                          SSP00227
      IF (IS.GT.0) GO TO 350                                            SSP00228
      IF (IFPR.EQ.0) GO TO 375                                          SSP00229
      WRITE (IOUT,1035)                                                 SSP00230
      WRITE (IOUT,1006) (EIGV(I),I=1,NC)                                SSP00231
C                                                                       SSP00232
C     CALCULATE B TIMES APPROXIMATE EIGENVECTORS (ICONV.EQ.0)           SSP00233
C        OR     FINAL EIGENVECTOR APPROXIMATIONS (ICONV.GT.0)           SSP00234
C                                                                       SSP00235
  375 DO 420 I=1,NN                                                     SSP00236
      DO 422 J=1,NC                                                     SSP00237
  422 TT(J)=R(I,J)                                                      SSP00238
      DO 424 K=1,NC                                                     SSP00239
      RT=0.                                                             SSP00240
      DO 430 L=1,NC                                                     SSP00241
  430 RT=RT + TT(L)*VEC(L,K)                                            SSP00242
  424 R(I,K)=RT                                                         SSP00243
  420 CONTINUE                                                          SSP00244
C                                                                       SSP00245
C     CALCULATE ERROR BOUNDS AND CHECK FOR CONVERGENCE OF EIGENVALUES   SSP00246
C                                                                       SSP00247
      DO 380 I=1,NC                                                     SSP00248
      VDOT=0.                                                           SSP00249
      DO 382 J=1,NC                                                     SSP00250
  382 VDOT=VDOT + VEC(I,J)*VEC(I,J)                                     SSP00251
      EIGV2=EIGV(I)*EIGV(I)                                             SSP00252
      DIF=VDOT - EIGV2                                                  SSP00253
      RDIF=MAX(DIF,TOLJ2*EIGV2)/EIGV2                                   SSP00254
      RDIF=SQRT(RDIF)                                                   SSP00255
      RTOLV(I)=RDIF                                                     SSP00256
  380 CONTINUE                                                          SSP00257
      IF (IFPR.EQ.0 .AND. ICONV.EQ.0) GO TO 385                         SSP00258
      WRITE (IOUT,1050)                                                 SSP00259
      WRITE (IOUT,1005) (RTOLV(I),I=1,NC)                               SSP00260
  385 IF (ICONV.GT.0) GO TO 500                                         SSP00261
C                                                                       SSP00262
      DO 390 I=1,NROOT                                                  SSP00263
      IF (RTOLV(I).GT.RTOL) GO TO 400                                   SSP00264
  390 CONTINUE                                                          SSP00265
      WRITE (IOUT,1060) RTOL                                            SSP00266
      ICONV=1                                                           SSP00267
      GO TO 100                                                         SSP00268
  400 IF (NITE.LT.NITEM) GO TO 100                                      SSP00269
      WRITE (IOUT,1070)                                                 SSP00270
      ICONV=2                                                           SSP00271
      IFSS=0                                                            SSP00272
      GO TO 100                                                         SSP00273
C                                                                       SSP00274
C - - - E N D   O F   I T E R A T I O N   L O O P                       SSP00275
C                                                                       SSP00276
  500 WRITE (IOUT,1100)                                                 SSP00277
      WRITE (IOUT,1006) (EIGV(I),I=1,NROOT)                             SSP00278
      WRITE (IOUT,1110)                                                 SSP00279
      DO 530 J=1,NROOT                                                  SSP00280
  530 WRITE (IOUT,1005) (R(K,J),K=1,NN)                                 SSP00281
C                                                                       SSP00282
C     CALCULATE AND PRINT ERROR MEASURES                                SSP00283
C                                                                       SSP00284
      REWIND NSTIF                                                      SSP00285
      READ (NSTIF) A                                                    SSP00286
C                                                                       SSP00287
      DO 580 L=1,NROOT                                                  SSP00288
      RT=EIGV(L)                                                        SSP00289
      CALL MULT(TT,A,R(1,L),MAXA,NN,NWK)                                SSP00290
      VNORM=0.                                                          SSP00291
      DO 590 I=1,NN                                                     SSP00292
  590 VNORM=VNORM + TT(I)*TT(I)                                         SSP00293
      CALL MULT(W,B,R(1,L),MAXA,NN,NWM)                                 SSP00294
      WNORM=0.                                                          SSP00295
      DO 600 I=1,NN                                                     SSP00296
      TT(I)=TT(I) - RT*W(I)                                             SSP00297
  600 WNORM=WNORM + TT(I)*TT(I)                                         SSP00298
      VNORM=SQRT(VNORM)                                                 SSP00299
      WNORM=SQRT(WNORM)                                                 SSP00300
      D(L)=WNORM/VNORM                                                  SSP00301
  580 CONTINUE                                                          SSP00302
      WRITE (IOUT,1115)                                                 SSP00303
      WRITE (IOUT,1005) (D(I),I=1,NROOT)                                SSP00304
C                                                                       SSP00305
C     APPLY STURM SEQUENCE CHECK                                        SSP00306
C                                                                       SSP00307
      IF (IFSS.EQ.0) GO TO 900                                          SSP00308
      CALL SCHECK (EIGV,RTOLV,BUP,BLO,BUPC,D,NC,NEI,RTOL,SHIFT,IOUT)    SSP00309
C                                                                       SSP00310
      WRITE (IOUT,1120) SHIFT                                           SSP00311
C                                                                       SSP00312
C     SHIFT MATRIX A                                                    SSP00313
C                                                                       SSP00314
      REWIND NSTIF                                                      SSP00315
      READ (NSTIF) A                                                    SSP00316
      IF (NWM.GT.NN) GO TO 645                                          SSP00317
      DO 640 I=1,NN                                                     SSP00318
      II=MAXA(I)                                                        SSP00319
  640 A(II)=A(II) - B(I)*SHIFT                                          SSP00320
      GO TO 660                                                         SSP00321
  645 DO 650 I=1,NWK                                                    SSP00322
  650 A(I)=A(I) - B(I)*SHIFT                                            SSP00323
C                                                                       SSP00324
C     FACTORIZE SHIFTED MATRIX                                          SSP00325
C                                                                       SSP00326
  660 ISH=1                                                             SSP00327
      CALL DECOMP (A,MAXA,NN,ISH,IOUT)                                  SSP00328
C                                                                       SSP00329
C     COUNT NUMBER OF NEGATIVE DIAGONAL ELEMENTS                        SSP00330
C                                                                       SSP00331
      NSCH=0                                                            SSP00332
      DO 664 I=1,NN                                                     SSP00333
      II=MAXA(I)                                                        SSP00334
      IF (A(II).LT.0.) NSCH=NSCH + 1                                    SSP00335
  664 CONTINUE                                                          SSP00336
      IF (NSCH.EQ.NEI) GO TO 670                                        SSP00337
      NMIS=NSCH - NEI                                                   SSP00338
      WRITE (IOUT,1130) NMIS                                            SSP00339
      GO TO 900                                                         SSP00340
  670 WRITE (IOUT,1140) NSCH                                            SSP00341
      GO TO 900                                                         SSP00342
C                                                                       SSP00343
  800 STOP                                                              SSP00344
  900 RETURN                                                            SSP00345
C                                                                       SSP00346
 1002 FORMAT (' ',10F10.0)                                              SSP00347
 1005 FORMAT (' ',12E11.4)                                              SSP00348
 1006 FORMAT (' ',6E22.14)                                              SSP00349
 1007 FORMAT (///,' STOP, NC IS LARGER THAN THE NUMBER OF MASS ',       SSP00350
     1        'DEGREES OF FREEDOM')                                     SSP00351
 1008 FORMAT (///,' DEGREES OF FREEDOM EXCITED BY UNIT STARTING ',      SSP00352
     1        'ITERATION VECTORS')                                      SSP00353
 1010 FORMAT (//,' I T E R A T I O N   N U M B E R ',I8)                SSP00354
 1020 FORMAT (/,' PROJECTION OF A (MATRIX AR)')                         SSP00355
 1030 FORMAT (/,' PROJECTION OF B (MATRIX BR)')                         SSP00356
 1035 FORMAT (/,' EIGENVALUES OF AR-LAMBDA*BR')                         SSP00357
 1040 FORMAT (//,' AR AND BR AFTER JACOBI DIAGONALIZATION')             SSP00358
 1050 FORMAT (/,' ERROR BOUNDS REACHED ON EIGENVALUES')                 SSP00359
 1060 FORMAT (///,' CONVERGENCE REACHED FOR RTOL ',E10.4)               SSP00360
 1070 FORMAT (' *** NO CONVERGENCE IN MAXIMUM NUMBER OF ITERATIONS',    SSP00361
     1        ' PERMITTED',/,                                           SSP00362
     2        ' WE ACCEPT CURRENT ITERATION VALUES',/,                  SSP00363
     3        ' THE STURM SEQUENCE CHECK IS NOT PERFORMED')             SSP00364
 1100 FORMAT (///,' THE CALCULATED EIGENVALUES ARE')                    SSP00365
 1115 FORMAT (//,' ERROR MEASURES ON THE EIGENVALUES')                  SSP00366
 1110 FORMAT (//,' THE CALCULATED EIGENVECTORS ARE',/)                  SSP00367
 1120 FORMAT (///,' CHECK APPLIED AT SHIFT ',E22.14)                    SSP00368
 1130 FORMAT (//,' THERE ARE ',I8,' EIGENVALUES MISSING')               SSP00369
 1140 FORMAT (//,' WE FOUND THE LOWEST ',I8,' EIGENVALUES')             SSP00370
C                                                                       SSP00371
      END                                                               SSP00372
      SUBROUTINE DECOMP (A,MAXA,NN,ISH,IOUT)                            SSP00373
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . SSP00374
C .                                                                   . SSP00375
C .   P R O G R A M                                                   . SSP00376
C .        TO CALCULATE (L)*(D)*(L)(T) FACTORIZATION OF               . SSP00377
C .        STIFFNESS MATRIX                                           . SSP00378
C .                                                                   . SSP00379
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . SSP00380
C                                                                       SSP00381
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               SSP00382
      DIMENSION A(1),MAXA(1)                                            SSP00383
      IF (NN.EQ.1) GO TO 900                                            SSP00384
C                                                                       SSP00385
      DO 200 N=1,NN                                                     SSP00386
      KN=MAXA(N)                                                        SSP00387
      KL=KN + 1                                                         SSP00388
      KU=MAXA(N+1) - 1                                                  SSP00389
      KH=KU - KL                                                        SSP00390
      IF (KH) 304,240,210                                               SSP00391
  210 K=N - KH                                                          SSP00392
      IC=0                                                              SSP00393
      KLT=KU                                                            SSP00394
      DO 260 J=1,KH                                                     SSP00395
      IC=IC + 1                                                         SSP00396
      KLT=KLT - 1                                                       SSP00397
      KI=MAXA(K)                                                        SSP00398
      ND=MAXA(K+1) - KI - 1                                             SSP00399
      IF (ND) 260,260,270                                               SSP00400
  270 KK=MIN0(IC,ND)                                                    SSP00401
      C=0.                                                              SSP00402
      DO 280 L=1,KK                                                     SSP00403
  280 C=C + A(KI+L)*A(KLT+L)                                            SSP00404
      A(KLT)=A(KLT) - C                                                 SSP00405
  260 K=K + 1                                                           SSP00406
  240 K=N                                                               SSP00407
      B=0.                                                              SSP00408
      DO 300 KK=KL,KU                                                   SSP00409
      K=K - 1                                                           SSP00410
      KI=MAXA(K)                                                        SSP00411
      C=A(KK)/A(KI)                                                     SSP00412
      IF (ABS(C).LT.1.E07) GO TO 290                                    SSP00413
      WRITE (IOUT,2010) N,C                                             SSP00414
      GO TO 800                                                         SSP00415
  290 B=B + C*A(KK)                                                     SSP00416
  300 A(KK)=C                                                           SSP00417
      A(KN)=A(KN) - B                                                   SSP00418
  304 IF (A(KN)) 310,310,200                                            SSP00419
  310 IF (ISH.EQ.0) GO TO 320                                           SSP00420
      IF (A(KN).EQ.0.) A(KN)=-1.E-16                                    SSP00421
      GO TO 200                                                         SSP00422
  320 WRITE (IOUT,2000) N,A(KN)                                         SSP00423
      GO TO 800                                                         SSP00424
  200 CONTINUE                                                          SSP00425
      GO TO 900                                                         SSP00426
C                                                                       SSP00427
  800 STOP                                                              SSP00428
  900 RETURN                                                            SSP00429
C                                                                       SSP00430
 2000 FORMAT (//' STOP - STIFFNESS MATRIX NOT POSITIVE DEFINITE',//,    SSP00431
     1          ' NONPOSITIVE PIVOT FOR EQUATION ',I8,//,               SSP00432
     2          ' PIVOT = ',E20.12)                                     SSP00433
 2010 FORMAT (//' STOP - STURM SEQUENCE CHECK FAILED BECAUSE OF',       SSP00434
     1          ' MULTIPLIER GROWTH FOR COLUMN NUMBER ',I8,//,          SSP00435
     2          ' MULTIPLIER = ',E20.8)                                 SSP00436
      END                                                               SSP00437
      SUBROUTINE REDBAK (A,V,MAXA,NN)                                   SSP00438
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . SSP00439
C .                                                                   . SSP00440
C .   P R O G R A M                                                   . SSP00441
C .        TO REDUCE AND BACK-SUBSTITUTE ITERATION VECTORS            . SSP00442
C .                                                                   . SSP00443
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . SSP00444
C                                                                       SSP00445
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               SSP00446
      DIMENSION A(1),V(1),MAXA(1)                                       SSP00447
C                                                                       SSP00448
      DO 400 N=1,NN                                                     SSP00449
      KL=MAXA(N) + 1                                                    SSP00450
      KU=MAXA(N+1) - 1                                                  SSP00451
      IF (KU-KL) 400,410,410                                            SSP00452
  410 K=N                                                               SSP00453
      C=0.                                                              SSP00454
      DO 420 KK=KL,KU                                                   SSP00455
      K=K - 1                                                           SSP00456
  420 C=C + A(KK)*V(K)                                                  SSP00457
      V(N)=V(N) - C                                                     SSP00458
  400 CONTINUE                                                          SSP00459
C                                                                       SSP00460
      DO 480 N=1,NN                                                     SSP00461
      K=MAXA(N)                                                         SSP00462
  480 V(N)=V(N)/A(K)                                                    SSP00463
      IF (NN.EQ.1) GO TO 900                                            SSP00464
      N=NN                                                              SSP00465
      DO 500 L=2,NN                                                     SSP00466
      KL=MAXA(N) + 1                                                    SSP00467
      KU=MAXA(N+1) - 1                                                  SSP00468
      IF (KU-KL) 500,510,510                                            SSP00469
  510 K=N                                                               SSP00470
      DO 520 KK=KL,KU                                                   SSP00471
      K=K - 1                                                           SSP00472
  520 V(K)=V(K) - A(KK)*V(N)                                            SSP00473
  500 N=N - 1                                                           SSP00474
C                                                                       SSP00475
  900 RETURN                                                            SSP00476
      END                                                               SSP00477
      SUBROUTINE MULT (TT,B,RR,MAXA,NN,NWM)                             SSP00478
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . SSP00479
C .                                                                   . SSP00480
C .   P R O G R A M                                                   . SSP00481
C .        TO EVALUATE PRODUCT OF B TIMES RR AND STORE RESULT IN TT   . SSP00482
C .                                                                   . SSP00483
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . SSP00484
C                                                                       SSP00485
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               SSP00486
      DIMENSION TT(1),B(1),RR(1),MAXA(1)                                SSP00487
C                                                                       SSP00488
      IF (NWM.GT.NN) GO TO 20                                           SSP00489
      DO 10 I=1,NN                                                      SSP00490
   10 TT(I)=B(I)*RR(I)                                                  SSP00491
      GO TO 900                                                         SSP00492
C                                                                       SSP00493
   20 DO 40 I=1,NN                                                      SSP00494
   40 TT(I)=0.                                                          SSP00495
      DO 100 I=1,NN                                                     SSP00496
      KL=MAXA(I)                                                        SSP00497
      KU=MAXA(I+1) - 1                                                  SSP00498
      II=I + 1                                                          SSP00499
      CC=RR(I)                                                          SSP00500
      DO 100 KK=KL,KU                                                   SSP00501
      II=II - 1                                                         SSP00502
  100 TT(II)=TT(II) + B(KK)*CC                                          SSP00503
      IF (NN.EQ.1) GO TO 900                                            SSP00504
      DO 200 I=2,NN                                                     SSP00505
      KL=MAXA(I) + 1                                                    SSP00506
      KU=MAXA(I+1) - 1                                                  SSP00507
      IF (KU-KL) 200,210,210                                            SSP00508
  210 II=I                                                              SSP00509
      AA=0.                                                             SSP00510
      DO 220 KK=KL,KU                                                   SSP00511
      II=II - 1                                                         SSP00512
  220 AA=AA + B(KK)*RR(II)                                              SSP00513
      TT(I)=TT(I) + AA                                                  SSP00514
  200 CONTINUE                                                          SSP00515
C                                                                       SSP00516
  900 RETURN                                                            SSP00517
      END                                                               SSP00518
      SUBROUTINE SCHECK (EIGV,RTOLV,BUP,BLO,BUPC,NEIV,NC,NEI,RTOL,      SSP00519
     1                   SHIFT,IOUT)                                    SSP00520
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . SSP00521
C .                                                                   . SSP00522
C .   P R O G R A M                                                   . SSP00523
C .        TO EVALUATE SHIFT FOR STURM SEQUENCE CHECK                 . SSP00524
C .                                                                   . SSP00525
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . SSP00526
C                                                                       SSP00527
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               SSP00528
      DIMENSION EIGV(NC),RTOLV(NC),BUP(NC),BLO(NC),BUPC(NC),NEIV(NC)    SSP00529
C                                                                       SSP00530
      FTOL=0.01                                                         SSP00531
C                                                                       SSP00532
      DO 100 I=1,NC                                                     SSP00533
      BUP(I)=EIGV(I)*(1.+FTOL)                                          SSP00534
  100 BLO(I)=EIGV(I)*(1.-FTOL)                                          SSP00535
      NROOT=0                                                           SSP00536
      DO 120 I=1,NC                                                     SSP00537
  120 IF (RTOLV(I).LT.RTOL) NROOT=NROOT + 1                             SSP00538
      IF (NROOT.GE.1) GO TO 200                                         SSP00539
      WRITE (IOUT,1010)                                                 SSP00540
      GO TO 800                                                         SSP00541
C                                                                       SSP00542
C      FIND UPPER BOUNDS ON EIGENVALUE CLUSTERS                         SSP00543
C                                                                       SSP00544
  200 DO 240 I=1,NROOT                                                  SSP00545
  240 NEIV(I)=1                                                         SSP00546
      IF (NROOT.NE.1) GO TO 260                                         SSP00547
      BUPC(1)=BUP(1)                                                    SSP00548
      LM=1                                                              SSP00549
      L=1                                                               SSP00550
      I=2                                                               SSP00551
      GO TO 295                                                         SSP00552
  260 L=1                                                               SSP00553
      I=2                                                               SSP00554
  270 IF (BUP(I-1).LE.BLO(I)) GO TO 280                                 SSP00555
      NEIV(L)=NEIV(L) + 1                                               SSP00556
      I=I + 1                                                           SSP00557
      IF (I.LE.NROOT) GO TO 270                                         SSP00558
  280 BUPC(L)=BUP(I-1)                                                  SSP00559
      IF (I.GT.NROOT) GO TO 290                                         SSP00560
      L=L + 1                                                           SSP00561
      I=I + 1                                                           SSP00562
      IF (I.LE.NROOT) GO TO 270                                         SSP00563
      BUPC(L)=BUP(I-1)                                                  SSP00564
  290 LM=L                                                              SSP00565
      IF (NROOT.EQ.NC) GO TO 300                                        SSP00566
  295 IF (BUP(I-1).LE.BLO(I)) GO TO 300                                 SSP00567
      IF (RTOLV(I).GT.RTOL) GO TO 300                                   SSP00568
      BUPC(L)=BUP(I)                                                    SSP00569
      NEIV(L)=NEIV(L) + 1                                               SSP00570
      NROOT=NROOT + 1                                                   SSP00571
      IF (NROOT.EQ.NC) GO TO 300                                        SSP00572
      I=I + 1                                                           SSP00573
      GO TO 295                                                         SSP00574
C                                                                       SSP00575
C      FIND SHIFT                                                       SSP00576
C                                                                       SSP00577
  300 WRITE (IOUT,1020)                                                 SSP00578
      WRITE (IOUT,1005) (BUPC(I),I=1,LM)                                SSP00579
      WRITE (IOUT,1030)                                                 SSP00580
      WRITE (IOUT,1006) (NEIV(I),I=1,LM)                                SSP00581
      LL=LM - 1                                                         SSP00582
      IF (LM.EQ.1) GO TO 310                                            SSP00583
  330 DO 320 I=1,LL                                                     SSP00584
  320 NEIV(L)=NEIV(L) + NEIV(I)                                         SSP00585
      L=L - 1                                                           SSP00586
      LL=LL - 1                                                         SSP00587
      IF (L.NE.1) GO TO 330                                             SSP00588
  310 WRITE (IOUT,1040)                                                 SSP00589
      WRITE (IOUT,1006) (NEIV(I),I=1,LM)                                SSP00590
      L=0                                                               SSP00591
      DO 340 I=1,LM                                                     SSP00592
      L=L + 1                                                           SSP00593
      IF (NEIV(I).GE.NROOT) GO TO 350                                   SSP00594
  340 CONTINUE                                                          SSP00595
  350 SHIFT=BUPC(L)                                                     SSP00596
      NEI=NEIV(L)                                                       SSP00597
      GO TO 900                                                         SSP00598
C                                                                       SSP00599
  800 STOP                                                              SSP00600
  900 RETURN                                                            SSP00601
C                                                                       SSP00602
 1005 FORMAT (' ',6E22.14)                                              SSP00603
 1006 FORMAT (' ',6I22)                                                 SSP00604
 1010 FORMAT (' *** ERROR ***  SOLUTION STOP IN *SCHECK*',/,            SSP00605
     1        ' NO EIGENVALUES FOUND',/)                                SSP00606
 1020 FORMAT (///,' UPPER BOUNDS ON EIGENVALUE CLUSTERS')               SSP00607
 1030 FORMAT (//,' NO. OF EIGENVALUES IN EACH CLUSTER')                 SSP00608
 1040 FORMAT (' NO. OF EIGENVALUES LESS THAN UPPER BOUNDS')             SSP00609
      END                                                               SSP00610
      SUBROUTINE JACOBI (A,B,X,EIGV,D,N,NWA,RTOL,NSMAX,IFPR,IOUT)       SSP00611
C ..................................................................... SSP00612
C .                                                                   . SSP00613
C .   P R O G R A M                                                   . SSP00614
C .        TO SOLVE THE GENERALIZED EIGENPROBLEM USING THE            . SSP00615
C .        GENERALIZED JACOBI ITERATION                               . SSP00616
C ..................................................................... SSP00617
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               SSP00618
      DIMENSION A(NWA),B(NWA),X(N,N),EIGV(N),D(N)                       SSP00619
C                                                                       SSP00620
C     INITIALIZE EIGENVALUE AND EIGENVECTOR MATRICES                    SSP00621
C                                                                       SSP00622
      N1=N + 1                                                          SSP00623
      II=1                                                              SSP00624
      DO 10 I=1,N                                                       SSP00625
      IF (A(II).GT.0. .AND. B(II).GT.0.) GO TO 4                        SSP00626
      WRITE (IOUT,2020) II,A(II),B(II)                                  SSP00627
      GO TO 800                                                         SSP00628
    4 D(I)=A(II)/B(II)                                                  SSP00629
      EIGV(I)=D(I)                                                      SSP00630
   10 II=II + N1 - I                                                    SSP00631
      DO 30 I=1,N                                                       SSP00632
      DO 20 J=1,N                                                       SSP00633
   20 X(I,J)=0.                                                         SSP00634
   30 X(I,I)=1.                                                         SSP00635
      IF (N.EQ.1) GO TO 900                                             SSP00636
C                                                                       SSP00637
C     INITIALIZE SWEEP COUNTER AND BEGIN ITERATION                      SSP00638
C                                                                       SSP00639
      NSWEEP=0                                                          SSP00640
      NR=N - 1                                                          SSP00641
   40 NSWEEP=NSWEEP + 1                                                 SSP00642
      IF (IFPR.EQ.1) WRITE (IOUT,2000) NSWEEP                           SSP00643
C                                                                       SSP00644
C     CHECK IF PRESENT OFF-DIAGONAL ELEMENT IS LARGE ENOUGH TO REQUIRE  SSP00645
C     ZEROING                                                           SSP00646
C                                                                       SSP00647
      EPS=(.01)**(NSWEEP*2)                                             SSP00648
      DO 210 J=1,NR                                                     SSP00649
      JP1=J + 1                                                         SSP00650
      JM1=J - 1                                                         SSP00651
      LJK=JM1*N - JM1*J/2                                               SSP00652
      JJ=LJK + J                                                        SSP00653
      DO 210 K=JP1,N                                                    SSP00654
      KP1=K + 1                                                         SSP00655
      KM1=K - 1                                                         SSP00656
      JK=LJK + K                                                        SSP00657
      KK=KM1*N - KM1*K/2 + K                                            SSP00658
      EPTOLA=(A(JK)/A(JJ))*(A(JK)/A(KK))                                SSP00659
      EPTOLB=(B(JK)/B(JJ))*(B(JK)/B(KK))                                SSP00660
      IF (EPTOLA.LT.EPS .AND. EPTOLB.LT.EPS) GO TO 210                  SSP00661
C                                                                       SSP00662
C     IF ZEROING IS REQUIRED, CALCULATE THE ROTATION MATRIX ELEMENTS CA SSP00663
C     AND CG                                                            SSP00664
C                                                                       SSP00665
      AKK=A(KK)*B(JK) - B(KK)*A(JK)                                     SSP00666
      AJJ=A(JJ)*B(JK) - B(JJ)*A(JK)                                     SSP00667
      AB=A(JJ)*B(KK) - A(KK)*B(JJ)                                      SSP00668
      SCALE=A(KK)*B(KK)                                                 SSP00669
      ABCH=AB/SCALE                                                     SSP00670
      AKKCH=AKK/SCALE                                                   SSP00671
      AJJCH=AJJ/SCALE                                                   SSP00672
      CHECK=(ABCH*ABCH+4.0*AKKCH*AJJCH)/4.0                             SSP00673
      IF (CHECK) 50,60,60                                               SSP00674
   50 WRITE (IOUT,2020) JJ,A(JJ),B(JJ)                                  SSP00675
      GO TO 800                                                         SSP00676
   60 SQCH=SCALE*SQRT(CHECK)                                            SSP00677
      D1=AB/2. + SQCH                                                   SSP00678
      D2=AB/2. - SQCH                                                   SSP00679
      DEN=D1                                                            SSP00680
      IF (ABS(D2).GT.ABS(D1)) DEN=D2                                    SSP00681
      IF (DEN) 80,70,80                                                 SSP00682
   70 CA=0.                                                             SSP00683
      CG=-A(JK)/A(KK)                                                   SSP00684
      GO TO 90                                                          SSP00685
   80 CA=AKK/DEN                                                        SSP00686
      CG=-AJJ/DEN                                                       SSP00687
C                                                                       SSP00688
C     PERFORM THE GENERALIZED ROTATION TO ZERO THE PRESENT OFF-DIAGONAL SSP00689
C     ELEMENT                                                           SSP00690
C                                                                       SSP00691
   90 IF (N-2) 100,190,100                                              SSP00692
  100 IF (JM1-1) 130,110,110                                            SSP00693
  110 DO 120 I=1,JM1                                                    SSP00694
      IM1=I - 1                                                         SSP00695
      IJ=IM1*N - IM1*I/2 + J                                            SSP00696
      IK=IM1*N - IM1*I/2 + K                                            SSP00697
      AJ=A(IJ)                                                          SSP00698
      BJ=B(IJ)                                                          SSP00699
      AK=A(IK)                                                          SSP00700
      BK=B(IK)                                                          SSP00701
      A(IJ)=AJ + CG*AK                                                  SSP00702
      B(IJ)=BJ + CG*BK                                                  SSP00703
      A(IK)=AK + CA*AJ                                                  SSP00704
  120 B(IK)=BK + CA*BJ                                                  SSP00705
  130 IF (KP1-N) 140,140,160                                            SSP00706
  140 LJI=JM1*N - JM1*J/2                                               SSP00707
      LKI=KM1*N - KM1*K/2                                               SSP00708
      DO 150 I=KP1,N                                                    SSP00709
      JI=LJI + I                                                        SSP00710
      KI=LKI + I                                                        SSP00711
      AJ=A(JI)                                                          SSP00712
      BJ=B(JI)                                                          SSP00713
      AK=A(KI)                                                          SSP00714
      BK=B(KI)                                                          SSP00715
      A(JI)=AJ + CG*AK                                                  SSP00716
      B(JI)=BJ + CG*BK                                                  SSP00717
      A(KI)=AK + CA*AJ                                                  SSP00718
  150 B(KI)=BK + CA*BJ                                                  SSP00719
  160 IF (JP1-KM1) 170,170,190                                          SSP00720
  170 LJI=JM1*N - JM1*J/2                                               SSP00721
      DO 180 I=JP1,KM1                                                  SSP00722
      JI=LJI + I                                                        SSP00723
      IM1=I - 1                                                         SSP00724
      IK=IM1*N - IM1*I/2 + K                                            SSP00725
      AJ=A(JI)                                                          SSP00726
      BJ=B(JI)                                                          SSP00727
      AK=A(IK)                                                          SSP00728
      BK=B(IK)                                                          SSP00729
      A(JI)=AJ + CG*AK                                                  SSP00730
      B(JI)=BJ + CG*BK                                                  SSP00731
      A(IK)=AK + CA*AJ                                                  SSP00732
  180 B(IK)=BK + CA*BJ                                                  SSP00733
  190 AK=A(KK)                                                          SSP00734
      BK=B(KK)                                                          SSP00735
      A(KK)=AK + 2.*CA*A(JK) + CA*CA*A(JJ)                              SSP00736
      B(KK)=BK + 2.*CA*B(JK) + CA*CA*B(JJ)                              SSP00737
      A(JJ)=A(JJ) + 2.*CG*A(JK) + CG*CG*AK                              SSP00738
      B(JJ)=B(JJ) + 2.*CG*B(JK) + CG*CG*BK                              SSP00739
      A(JK)=0.                                                          SSP00740
      B(JK)=0.                                                          SSP00741
C                                                                       SSP00742
C     UPDATE THE EIGENVECTOR MATRIX AFTER EACH ROTATION                 SSP00743
C                                                                       SSP00744
      DO 200 I=1,N                                                      SSP00745
      XJ=X(I,J)                                                         SSP00746
      XK=X(I,K)                                                         SSP00747
      X(I,J)=XJ + CG*XK                                                 SSP00748
  200 X(I,K)=XK + CA*XJ                                                 SSP00749
  210 CONTINUE                                                          SSP00750
C                                                                       SSP00751
C     UPDATE THE EIGENVALUES AFTER EACH SWEEP                           SSP00752
C                                                                       SSP00753
      II=1                                                              SSP00754
      DO 220 I=1,N                                                      SSP00755
      IF (A(II).GT.0. .AND. B(II).GT.0.) GO TO 215                      SSP00756
      WRITE (IOUT,2020) II,A(II),B(II)                                  SSP00757
      GO TO 800                                                         SSP00758
  215 EIGV(I)=A(II)/B(II)                                               SSP00759
  220 II=II + N1 - I                                                    SSP00760
      IF (IFPR.EQ.0) GO TO 230                                          SSP00761
      WRITE (IOUT,2030)                                                 SSP00762
      WRITE (IOUT,2010) (EIGV(I),I=1,N)                                 SSP00763
C                                                                       SSP00764
C     CHECK FOR CONVERGENCE                                             SSP00765
C                                                                       SSP00766
  230 DO 240 I=1,N                                                      SSP00767
      TOL=RTOL*D(I)                                                     SSP00768
      DIF=ABS(EIGV(I)-D(I))                                             SSP00769
      IF (DIF.GT.TOL) GO TO 280                                         SSP00770
  240 CONTINUE                                                          SSP00771
C                                                                       SSP00772
C     CHECK ALL OFF-DIAGONAL ELEMENTS TO SEE IF ANOTHER SWEEP IS        SSP00773
C     REQUIRED                                                          SSP00774
C                                                                       SSP00775
      EPS=RTOL**2                                                       SSP00776
      DO 250 J=1,NR                                                     SSP00777
      JM1=J - 1                                                         SSP00778
      JP1=J + 1                                                         SSP00779
      LJK=JM1*N - JM1*J/2                                               SSP00780
      JJ=LJK + J                                                        SSP00781
      DO 250 K=JP1,N                                                    SSP00782
      KM1=K - 1                                                         SSP00783
      JK=LJK + K                                                        SSP00784
      KK=KM1*N - KM1*K/2 + K                                            SSP00785
      EPSA=(A(JK)/A(JJ))*(A(JK)/A(KK))                                  SSP00786
      EPSB=(B(JK)/B(JJ))*(B(JK)/B(KK))                                  SSP00787
      IF (EPSA.LT.EPS .AND. EPSB.LT.EPS) GO TO 250                      SSP00788
      GO TO 280                                                         SSP00789
  250 CONTINUE                                                          SSP00790
C                                                                       SSP00791
C     SCALE EIGENVECTORS                                                SSP00792
C                                                                       SSP00793
  255 II=1                                                              SSP00794
      DO 275 I=1,N                                                      SSP00795
      BB=SQRT(B(II))                                                    SSP00796
      DO 270 K=1,N                                                      SSP00797
  270 X(K,I)=X(K,I)/BB                                                  SSP00798
  275 II=II + N1 - I                                                    SSP00799
      GO TO 900                                                         SSP00800
C                                                                       SSP00801
C     UPDATE  D  MATRIX AND START NEW SWEEP, IF ALLOWED                 SSP00802
C                                                                       SSP00803
  280 DO 290 I=1,N                                                      SSP00804
  290 D(I)=EIGV(I)                                                      SSP00805
      IF (NSWEEP.LT.NSMAX) GO TO 40                                     SSP00806
      GO TO 255                                                         SSP00807
C                                                                       SSP00808
  800 STOP                                                              SSP00809
  900 RETURN                                                            SSP00810
C                                                                       SSP00811
 2000 FORMAT (//,' SWEEP NUMBER IN *JACOBI* = ',I8)                     SSP00812
 2010 FORMAT (' ',6E20.12)                                              SSP00813
 2020 FORMAT (' *** ERROR *** SOLUTION STOP',/,                         SSP00814
     1        ' MATRICES NOT POSITIVE DEFINITE',/,                      SSP00815
     2        ' II = ',I8,' A(II) = ',E20.12,' B(II) = ',E20.12)        SSP00816
 2030 FORMAT (/,' CURRENT EIGENVALUES IN *JACOBI* ARE',/)               SSP00817
      END                                                               SSP00818
