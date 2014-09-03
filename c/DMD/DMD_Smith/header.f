<PRE>
	COMMON /PHSSPC / RX(NOP),RY(NOP),RZ(NOP),
     &                   VX(NOP),VY(NOP),VZ(NOP)
	COMMON /INTGRS / IN(4)
	COMMON /REALS /  RL(29)

C** Make this a separate array for those lucky 64 bit machines
C** and their users who want to exceed the 2 billion limit of 32 bit integers
        INTEGER INDBLE
	COMMON /INTDBL  / INDBLE(2)
	INTEGER NCOLL, COLL

C** Chain number for each segment
        COMMON /CHNUM / ICHN(NOP)

C** Using equivalence statements to reduce variable passing
	EQUIVALENCE (INDBLE(1),NCOLL), (INDBLE(2),COLL)

	INTEGER CHNLEN
	EQUIVALENCE (IN(1),CHNLEN), (IN(2),NCHN)
	EQUIVALENCE (IN(3),NABLIM), (IN(4), NLIMIT)

	EQUIVALENCE (RL(1),SIGMA), (RL(2),SIGSQ), (RL(3),DELFAC)
	EQUIVALENCE (RL(4),BONDL), (RL(5),BNDLSQ), (RL(6),SIGSQIJ)
	EQUIVALENCE (RL(7),RLIST), (RL(8),ETA), (RL(9),TEMP)
	EQUIVALENCE (RL(10),T), (RL(11),TNEXT), (RL(12),W) 
	EQUIVALENCE (RL(13),BOXL), (RL(14),VOL), (RL(15),BOXLM1)
	EQUIVALENCE (RL(16),RLSTSQ), (RL(17),DSPLIM) 
	EQUIVALENCE (RL(18),RCHLEN), (RL(19), CELLI), (RL(20), BOXHLF)
C** Square well items
        EQUIVALENCE (RL(21),PE), (RL(22),TOTE),(RL(23),TSTAR) 
	EQUIVALENCE (RL(24),V0), (RL(25),DELPE), (RL(26),CIJ2)
	EQUIVALENCE (RL(27),SIG2SQIJ),(RL(28),SIG2SQ),(RL(29),SIG2)
</PRE>
