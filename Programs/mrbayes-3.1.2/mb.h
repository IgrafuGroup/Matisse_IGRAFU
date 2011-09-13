#if !defined(UNIX_VERSION) && !defined(WIN_VERSION) && !defined(MAC_VERSION)
#ifdef __MWERKS__
#define MAC_VERSION
#else
#define WIN_VERSION
#endif
#endif

/* definition of UNIX_VERSION, WIN_VERSION or MAC_VERSION is now set in the
   Makefile file; for compilers not using a project file or Makefile, use the #defines above to select version */

typedef	double MrBFlt;		/* double used for parameter values and generally for floating point values */
typedef float CLFlt;		/* single-precision float used for cond likes (CLFlt) to increase speed and reduce memory requirement */
							/* set CLFlt to double if you want increased precision */

#undef FAST_LOG			/* define this to get a fast polynomial approximation to the log for 4by4 nucleotide data */
						/* this results in speed gains of 10% to 20% for typical data sets */
						/* however, it does result in less precision and the code does not work with all compilers */

#undef SSE		/* DO NOT #define THESE: ONLY FOR TESTING PURPOSES FOR NOW,    */
#undef ALTIVEC  /* NOT YET FUNCTIONAL											*/

/* These statements are needed for the SIMD code tessting ONLY, DO NOT #DEFINE THESE */
/*
#undef MS_VCPP
#undef GCC
#undef ICC
*/

#if defined GCC			/* gcc compiler */
#define ALIGNEDMALLOC memalign
#elif defined ICC		/* icc compiler */
#define ALIGNEDMALLOC _mm_malloc
#elif defined MS_VCPP   /* Visual .Net */
#define ALIGNEDMALLOC _aligned_malloc
#include <xmmintrin.h>
#endif

#ifdef PAUL /* for debugging */
#define DEBUG(fmt, arg) printf("%s:%d ",__FILE__,__LINE__);printf(fmt,arg);
#define free(a) free(a);a=NULL
#else
#define DEBUG(a,b) 
#endif
extern void *SafeMalloc(size_t);
extern int SafeFclose(FILE **);

/* For comparing floating points: two values are the same if the absolute difference is less then 
   this value. 
*/
#ifndef ETA
#define ETA (1E-30)
#endif

#if defined (MPI_ENABLED)
#	if defined (__MWERKS__) & defined (MAC_VERSION)
#		include "macmpi.h"
#	else
#		include "mpi.h"
#	endif
#endif

#undef ASYMMETRY

#define RELEASE
#ifdef RELEASE
#define	VERSION_NUMBER			"3.1.2"
#else
#define VERSION_NUMBER                  "3.1.2-cvs"
#endif

#undef NO_ERROR
#undef ERROR
#define NO_ERROR				0
#define ERROR					1
#define	NO_ERROR_QUIT			2
#define ABORT					3

#undef FALSE
#undef TRUE
#define FALSE					0
#define TRUE					1

#define NO						0
#define YES						1

#define UP						0
#define DOWN					1

#define	UPPER					0
#define	MIDDLE					1
#define	LOWER					2

#define NONINTERACTIVE			0
#define INTERACTIVE				1

#define STANDARD_USER           1
#define DEVELOPER               3

#define	DIFFERENT				0
#define	SAME					1
#define	CONSISTENT_WITH			2

#define	LINETERM_UNIX			0
#define	LINETERM_MAC			1
#define	LINETERM_DOS			2

#define	SCREENWIDTH				60
#define	SCREENWIDTH2			61

#define	NONE					0
#define	DNA						1
#define	RNA						2
#define	PROTEIN					3
#define	RESTRICTION				4
#define	STANDARD				5
#define	MIXED					6
#define	CONTINUOUS				7

#define	AAMODEL_POISSON			0
#define	AAMODEL_JONES			1
#define	AAMODEL_DAY				2
#define	AAMODEL_MTREV			3
#define	AAMODEL_MTMAM			4
#define	AAMODEL_WAG				5
#define	AAMODEL_RTREV			6
#define	AAMODEL_CPREV			7
#define	AAMODEL_VT				8
#define	AAMODEL_BLOSUM			9
#define	AAMODEL_EQ				10
#define	AAMODEL_GTR				11 /* aa models with free parameters must be listed last */

#define	MISSING					10000000
#define	GAP						10000001

#define	UNORD					0
#define	ORD						1
#define	DOLLO					2
#define	IRREV					3

#define	IN_CMD					0
#define	IN_FILE					1

#define	NOTHING					0
#define	COMMAND					1
#define	PARAMETER				2
#define	EQUALSIGN				3
#define	COLON					4
#define	SEMICOLON				5
#define	COMMA					6
#define	POUNDSIGN				7
#define	QUESTIONMARK			8
#define	DASH					9
#define	LEFTPAR					10
#define	RIGHTPAR				11
#define	LEFTCOMMENT				12
#define	RIGHTCOMMENT			13
#define	ALPHA					14
#define	NUMBER					15
#define	RETURNSYMBOL			16
#define	ASTERISK				17
#define	BACKSLASH				18
#define	FORWARDSLASH			19
#define	EXCLAMATIONMARK			20
#define	PERCENT					21
#define	QUOTATIONMARK			22
#define	WEIRD					23
#define	UNKNOWN_TOKEN_TYPE		24
#define LEFTCURL				25
#define RIGHTCURL				26

#define	MAX_Q_RATE				100.0f
#define	MIN_SHAPE_PARAM			0.0001f
#define	MAX_SHAPE_PARAM			200.0f
#define	MAX_SITE_RATE			10.0f
#define	MAX_GAMMA_CATS			20
#define	MAX_GAMMA_CATS_SQUARED	400
#define	BRLENS_MIN				0.000001f
#define	BRLENS_MAX				100.0f
#define TREEHEIGHT_MAX			100000.0f
#define KAPPA_MIN				0.01f
#define	KAPPA_MAX				1000.0f
#define	GROWTH_MIN				-1000000.0f
#define	GROWTH_MAX				1000000.0f
#define RATE_MIN				0.000001f
#define RATE_MAX				100.0f
#define SYMPI_MIN				0.000001f
#define	SYMPI_MAX				100.0f
#define ALPHA_MIN				0.0001f
#define ALPHA_MAX				10000.0f
#define DIR_MIN					0.000001f
#define MIN_OFFSET_EXP_LAMBDA   0.000001f
#define MAX_OFFSET_EXP_LAMBDA   100000.0f

#define NEG_INFINITY			-1000000.0f

#define	CMD_STRING_LENGTH		20000

#define	pos(i,j,n)				((i)*(n)+(j))

#define	NUM_ALLOCS				 83

#define	ALLOC_MATRIX			 0
#define	ALLOC_CHARINFO			 2
#define	ALLOC_CHARSETNAMES		 3
#define	ALLOC_TAXANAMES			 4
#define	ALLOC_TMPSET			 5
#define	ALLOC_PARTITIONNAMES	 6
#define	ALLOC_TAXAINFO			 7
#define	ALLOC_TAXASETNAMES		 8
#define	ALLOC_CONSTRAINTNAMES	 9
#define	ALLOC_USERTREE			 10
#define	ALLOC_SUMTTREE			 11
#define	ALLOC_SUMTSTRING		 12
#define	ALLOC_SUMTDP			 13
#define	ALLOC_TREEBITS			 14
#define	ALLOC_TREEPARTS			 15
#define	ALLOC_NUMOFPART			 16
#define	ALLOC_ABRLENS			 17
#define	ALLOC_SUMB      		 18
#define	ALLOC_TAXONMASK			 19
#define	ALLOC_SBRLENS			 20
#define	ALLOC_CONNODES			 21
#define	ALLOC_OUTPART			 22
#define	ALLOC_TRANSFROM			 23
#define	ALLOC_TRANSTO			 24
#define	ALLOC_AVAILNODES		 25
#define	ALLOC_AVAILINDICES		 26
#define ALLOC_SCALERS			 27
#define	ALLOC_CURLNL			 28
#define	ALLOC_CURLNPR			 29
#define	ALLOC_CHAINID			 30
#define	ALLOC_PARAMS			 31
#define	ALLOC_TREE				 32
#define	ALLOC_NODES				 33
#define	ALLOC_LOCTAXANAMES		 34
#define	ALLOC_PARSMATRIX		 35
#define	ALLOC_TERMSTATE			 36
#define	ALLOC_ISPARTAMBIG		 37
#define	ALLOC_TERMCONDLIKES		 38
#define	ALLOC_COMPMATRIX		 39
#define	ALLOC_NUMSITESOFPAT		 40
#define ALLOC_COMPCOLPOS		 41
#define	ALLOC_COMPCHARPOS		 42
#define ALLOC_ORIGCHAR			 43
#define ALLOC_CHAINCONDLIKES	 44
#define ALLOC_NODEPTRS			 45
#define ALLOC_PARAMVALUES		 46
#define ALLOC_MCMCTREES			 47
#define ALLOC_MOVES				 48
#define ALLOC_CLSCALERS			 49
#define ALLOC_INVCONDLIKES		 50
#define ALLOC_TIPROBS			 51
#define	ALLOC_PRELIKES			 52
#define	ALLOC_CIJK				 53
#define ALLOC_SITEJUMP			 54
#define ALLOC_MARKOVTIS			 55
#define ALLOC_RATEPROBS			 56
#define ALLOC_STDTYPE			 57
#define	ALLOC_TAXAFOUND			 58
#define	ALLOC_FULLTREEINFO		 59
#define	ALLOC_PARTORIGORDER		 60
#define	ALLOC_PRUNEINFO			 61
#define	ALLOC_SUMPSTRING		 62
#define	ALLOC_SUMPINFO			 63
#define	ALLOC_SWAPINFO			 64
#define ALLOC_SYMPIINDEX		 65
#define	ALLOC_POSSELPROBS		 66
#define ALLOC_PARSSETS			 67
#define	ALLOC_PBF				 68
#define ALLOC_LOCALTAXONCALIBRATION		 69
#define	ALLOC_FULLCOMPTREEINFO	 70
#define	ALLOC_TOPO_DISTANCES	 71
#define	ALLOC_SPR_PARSSETS		 72
#define ALLOC_ANCSTATECONDLIKES  73
#define ALLOC_PFCOUNTERS         74
#define ALLOC_FILEPOINTERS       75
#define	ALLOC_STATS				 76
#define ALLOC_DIAGNUTREE         77
#define ALLOC_DIAGNRTREE         78
#define ALLOC_NUMINRUNOFPART     79
#define ALLOC_A_WITHIN_BRLENS    80
#define ALLOC_SUMSQB             81
#define ALLOC_S_WITHIN_BRLENS    82

#define	NUM_LINKED				18
#define	LINKED					0
#define	UNLINKED				1
#define	P_TRATIO				0
#define	P_REVMAT				1
#define	P_OMEGA					2
#define	P_PI					3
#define	P_SHAPE					4
#define	P_PINVAR				5
#define	P_CORREL				6
#define	P_SWITCH				7
#define	P_RATEMULT				8
#define	P_TOPOLOGY				9
#define	P_BRLENS				10
#define	P_SPECRATE				11
#define	P_EXTRATE				12
#define	P_THETA					13
#define	P_AAMODEL				14
#define	P_BRCORR				15
#define	P_BRSIGMA				16
#define	P_GROWTH				17 /* NOTE: If you add another parameter, change NUM_LINKED */

#define	MAX_NUM_DIVS			150
#define	MAX_NUM_CHARSETS		150
#define	MAX_NUM_PARTITIONS		30
#define	MAX_NUM_TREES			200
#define	MAX_NUM_DIV_BITS		150
#define MAX_NUM_DIV_LONGS       1 + MAX_NUM_DIV_BITS / 32
#define	MAX_NUM_TAXASETS		30
#define	MAX_NUM_CONSTRAINTS		30
#define	MAX_CHAINS				256

typedef void * VoidPtr;
typedef int (*CmdFxn)(void);
typedef int (*ParmFxn)(char *, char *);

typedef struct
	{
	MrBFlt			sum;
	MrBFlt			numPartitions;
	MrBFlt			numSamples;
	MrBFlt			avgStdDev;
	MrBFlt			**pair;
	} STATS;

enum CALPRIOR
	{
	unconstrained,
	fixed,
	offsetExponential,
	uniform
	};

typedef struct calibration
	{
	char			name[30];
	enum CALPRIOR   prior;
	MrBFlt			upper;
	MrBFlt			lower;
	MrBFlt			offset;
	MrBFlt			lambda;
	MrBFlt			age;
	}
	Calibration;
	
/* NOTE: Any variable added here must also be copied in CopyTrees */
typedef struct node
	{
	struct node		*left, *right, *anc;
	int				memoryIndex, index, upDateCl, upDateTi, marked, x, y,
					scalerNode, isLocked, lockID, uL, dL, mL, isDated;
	long int		scalersSet[MAX_NUM_DIV_LONGS], clSpace[MAX_NUM_DIV_LONGS], tiSpace[MAX_NUM_DIV_LONGS];
	long int 		*partition;
	char			label[100];
	MrBFlt			length, nodeDepth, d, age;
	Calibration		*calibration;
	}
	TreeNode;

typedef struct 
	{
	int				nNodes;
	int				nIntNodes;
	int				isRooted;
	int				isCalibrated;
	int				nRelParts;
	int				*relParts;
	int				checkConstraints;
	int				*constraints;
	int				nConstraints;
	int				nLocks;
	TreeNode		**allDownPass;
	TreeNode		**intDownPass;
	TreeNode		*root;
	TreeNode		*nodes;
	MrBFlt			clockRate;
	}
	Tree;

	/* typedef for nodes in polytomous tree */
typedef struct pNode
	{
	struct pNode	*left, *sib, *anc;
	int				x, y, mark, index, memoryIndex, isLocked, lockID, isDated;
	MrBFlt			length, support, f, age;
	char			label[100];
	long int		*partition;
	Calibration		*calibration;
	}
	PolyNode;


/* typedef for polytomous tree */
typedef struct 
	{
	int				nNodes;
	int				nIntNodes;
	PolyNode		**allDownPass;
	PolyNode		**intDownPass;
	PolyNode		*root;
	PolyNode		*nodes;
	}
	PolyTree;

/* struct for holding model parameter info for the mcmc run */
typedef struct param
	{
	int				index;			/* index to the parameter (0, 1, 2, ...)        */
	int				paramType;		/* the type of the parameter					*/
	int				paramId;		/* unique ID for parameter x prior combination	*/
	MrBFlt			*values;		/* main values of parameter						*/
	MrBFlt			*subValues;		/* subvalues of parameter						*/
	int				nValues;		/* number of values								*/
	int				nSubValues;		/* number of subvalues							*/	
	int				*relParts;		/* pointer to relevant divisions				*/
	int				nRelParts;		/* number of relevant divisions					*/
	int				upDate;			/* update flag (for copying)					*/
	struct param	**subParams;	/* pointers to subparams (for topology)			*/
	int				nSubParams;		/* number of subparams							*/
	Tree			*tree;			/* pointer to first tree (for brlens & topology)*/
	int				treeIndex;		/* index to first tree in mcmcTree				*/
	MrBFlt			relProposalProb;/* the relProposalProb for the parameter		*/
	int				*sympiBsIndex;	/* pointer to sympi bsIndex (std chars)			*/
	int				*sympinStates;	/* pointer to sympi nStates (std chars)			*/
	int				*sympiCType;	/* pointer to sympi cType (std chars)			*/
	int				nSympi;			/* number of sympis								*/
	int				printParam;     /* whether parameter should be printed          */
	char			paramName[1000];/* a string holding the name of the parameter(s)*/
	} Param;

/* parameter ID values */
/* identifies unique model parameter x prior combinations */
#define TRATIO_DIR						1
#define TRATIO_FIX						2
#define REVMAT_DIR						3
#define REVMAT_FIX						4
#define	OMEGA_DIR						5
#define	OMEGA_FIX						6
#define SYMPI_UNI						7
#define	SYMPI_UNI_MS					8
#define	SYMPI_EXP						9
#define SYMPI_EXP_MS					10
#define	SYMPI_FIX						11
#define	SYMPI_FIX_MS					12
#define SYMPI_EQUAL						13
#define	PI_DIR							14
#define	PI_USER							15
#define	PI_EMPIRICAL					16
#define	PI_EQUAL						17
#define	PI_FIXED						18
#define	SHAPE_UNI						19
#define SHAPE_EXP						20
#define SHAPE_FIX						21
#define PINVAR_UNI						22
#define PINVAR_FIX						23
#define	CORREL_UNI						24
#define CORREL_FIX						25
#define SWITCH_UNI						26
#define SWITCH_EXP						27
#define SWITCH_FIX						28
#define RATEMULT_DIR					29
#define RATEMULT_FIX					30
#define TOPOLOGY_NCL_UNIFORM			31
#define TOPOLOGY_NCL_CONSTRAINED		32
#define TOPOLOGY_CL_UNIFORM				33
#define TOPOLOGY_CL_CONSTRAINED			34
#define TOPOLOGY_CCL_UNIFORM			35
#define TOPOLOGY_CCL_CONSTRAINED		36
#define TOPOLOGY_NCL_UNIFORM_HOMO		37
#define TOPOLOGY_NCL_UNIFORM_HETERO		38
#define TOPOLOGY_NCL_CONSTRAINED_HOMO	39
#define TOPOLOGY_NCL_CONSTRAINED_HETERO 40
#define TOPOLOGY_CL_UNIFORM_HOMO		41
#define TOPOLOGY_CL_UNIFORM_HETERO		42
#define TOPOLOGY_CL_CONSTRAINED_HOMO	43
#define TOPOLOGY_CL_CONSTRAINED_HETERO	44
#define TOPOLOGY_CCL_UNIFORM_HOMO		45
#define TOPOLOGY_CCL_UNIFORM_HETERO		46
#define TOPOLOGY_CCL_CONSTRAINED_HOMO	47
#define TOPOLOGY_CCL_CONSTRAINED_HETERO	48
#define	TOPOLOGY_PARSIMONY_UNIFORM		49
#define	TOPOLOGY_PARSIMONY_CONSTRAINED	50
#define BRLENS_CLOCK_UNI				51
#define BRLENS_CLOCK_COAL				52
#define BRLENS_CLOCK_BD					53
#define BRLENS_CCLOCK_UNI				54
#define BRLENS_CCLOCK_COAL				55
#define BRLENS_CCLOCK_BD				56
#define BRLENS_UNI						57
#define BRLENS_EXP						58
#define	BRLENS_PARSIMONY				59
#define SPECRATE_UNI					60
#define SPECRATE_EXP					61
#define SPECRATE_FIX					62
#define EXTRATE_UNI						63
#define EXTRATE_EXP						64
#define EXTRATE_FIX						65
#define THETA_UNI						66
#define THETA_EXP						67
#define THETA_FIX						68
#define	AAMODEL_FIX						69
#define	AAMODEL_MIX						70
#define GROWTH_UNI						71
#define GROWTH_EXP						72
#define GROWTH_FIX						73
#define	GROWTH_NORMAL					74
#define	OMEGA_BUD						75
#define	OMEGA_BUF						76
#define	OMEGA_BED						77
#define	OMEGA_BEF						78
#define	OMEGA_BFD						79
#define	OMEGA_BFF						80
#define	OMEGA_FUD						81
#define	OMEGA_FUF						82
#define	OMEGA_FED						83
#define	OMEGA_FEF						84
#define	OMEGA_FFD						85
#define	OMEGA_FFF						86
#define	OMEGA_ED						87
#define	OMEGA_EF						88
#define	OMEGA_FD						89
#define	OMEGA_FF						90
#define	OMEGA_10UUB						91
#define	OMEGA_10UUF						92
#define	OMEGA_10UEB						93
#define	OMEGA_10UEF						94
#define	OMEGA_10UFB						95
#define	OMEGA_10UFF						96
#define	OMEGA_10EUB						97
#define	OMEGA_10EUF						98
#define	OMEGA_10EEB						99
#define	OMEGA_10EEF						100
#define	OMEGA_10EFB						101
#define	OMEGA_10EFF						102
#define	OMEGA_10FUB						103
#define	OMEGA_10FUF						104
#define	OMEGA_10FEB						105
#define	OMEGA_10FEF						106
#define	OMEGA_10FFB						107
#define	OMEGA_10FFF						108

/* typedef for a MoveFxn */
typedef int (MoveFxn)(Param *param, int chain, long int *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);

/* struct holding info on each move type that the program handles */
typedef struct
	{
	MoveFxn		*moveFxn;			/* pointer to the move function					*/
	int			nApplicable;		/* number of relevant params					*/
	int			applicableTo[40];	/* pointer to ID of relevant params				*/
	char		*name;				/* name of the move type						*/
	char		*nameTuning[2];		/* name of tuning params                        */
	char		*shortName;         /* abbreviated name of the move type            */
	MrBFlt		relProposalProb;	/* this holds the set proposal probability      */
	int         numTuningParams;    /* number of tuning parameters                  */
	MrBFlt		tuningParam[2];	    /* tuning parameters for the proposal           */
	MrBFlt		minimum[2];         /* minimum values for tuning params             */
	MrBFlt		maximum[2];         /* maximum values for tuning params             */
	int         parsimonyBased;     /* this move is based on parsimony (YES/NO)     */
	int         level;              /* user level of this move                      */
	} MoveType;


/* max number of move types */
#define NUM_MOVE_TYPES 100


/* struct holding info on each move */
/* Note: This allows several different moves to affect the same parameter */
/* It also allows the same move to affect different parameters as before */
/* This is also a good place to keep the proposal probs */
typedef struct
	{
	char		*name;				/* pointer to the name of the move type         */
	char		*shortName;	        /* pointer to the short name of the move        */
	MoveFxn		*moveFxn;			/* pointer to the move function					*/
	Param		*parm;				/* ptr to parameter the move applies to			*/
	MrBFlt		relProposalProb;	/* the actual proposal probability used			*/
	MrBFlt		cumProposalProb;	/* the cumulative proposal probability			*/
	int			*nAccepted;			/* number of accepted moves						*/
	int			*nTried;			/* number of tried moves						*/
	MrBFlt		tuningParam[2];     /* tuning parameters for the move               */
	} MCMCMove;

typedef int (*LikeDownFxn)(TreeNode *, int, int);
typedef int (*LikeRootFxn)(TreeNode *, int, int);
typedef int (*LikeScalerFxn)(TreeNode *, int, int);
typedef int (*LikeFxn)(TreeNode *, int, int, MrBFlt *, int);
typedef int (*TiProbFxn)(TreeNode *, int, int);
typedef int (*LikeUpFxn)(TreeNode *, int, int);
typedef int (*PrintAncStFxn)(TreeNode *, int, int);
typedef int (*StateCodeFxn) (int);
typedef int (*PrintSiteRateFxn) (TreeNode *, int, int);


typedef struct cmdtyp			
	{
	int			cmdNumber;
	char		*string;
	int			specialCmd;
	CmdFxn		cmdFxnPtr;
	short		numParms;
	short		parmList[40];
	int			expect;
	char		*cmdDescription;
	int			cmdUse;
	int			hiding;
	}
	CmdType;
	
typedef struct parm
	{
	char		*string;	/* parameter name */
	char		*valueList;	/* list of values that could be input */
	ParmFxn 	fp;	        /* function pointer */
	}
	ParmInfo, *ParmInfoPtr;

typedef struct model
	{
	int			dataType;          /* data type for partition                      */
	int			nStates;           /* number of states for this type of data       */
	int			codon[64];         /* gives protein ID for each codon              */
	int			codonNucs[64][3];  /* gives triplet for each codon                 */
	int			codonAAs[64];      /* gives protein ID for implemented code        */
	
	char		nucModel[100];     /* nucleotide model used                        */
	char		nst[100];          /* number of substitution types                 */
	char		parsModel[100];    /* use the (so-called) parsimony model          */
	char		geneticCode[100];  /* genetic code used                            */
	char		coding[100];       /* type of patterns encoded                     */
	char		ploidy[100];       /* ploidy level                                 */
	char		omegaVar[100];     /* type of omega variation model                */
	char		ratesModel[100];   /* rates across sites model                     */
	int 		numGammaCats;      /* number of categories for gamma approximation */
	int			numBetaCats;       /* number of categories for beta approximation  */
	int 		numM10GammaCats;   /* number of cats for gamma approx (M10 model)  */
	int			numM10BetaCats;    /* number of cats for beta approx (M10 model)   */
	char		covarionModel[100];/* use covarion model? (yes/no)                 */
	char		augmentData[100];  /* should data be augmented                     */

	char		tRatioPr[100];     /* prior for ti/tv rate ratio                   */
	MrBFlt		tRatioFix;   
	MrBFlt		tRatioDir[2];      
	char		revMatPr[100];     /* prior for GTR model                          */
	MrBFlt		revMatFix[6];
	MrBFlt		revMatDir[6];
	char		aaModelPr[100];    /* prior for amino acid model                   */
	char		aaModel[100];
	MrBFlt		aaModelPrProbs[10];
	char		aaRevMatPr[100];     /* prior for aa GTR model                     */
	MrBFlt		aaRevMatFix[190];
	MrBFlt		aaRevMatDir[190];
	char		omegaPr[100];      /* prior for omega                              */
	MrBFlt		omegaFix;
	MrBFlt		omegaDir[2];
	char		ny98omega1pr[100]; /* prior for class 1 omega (Ny98 model)         */
	MrBFlt		ny98omega1Fixed;
	MrBFlt		ny98omega1Beta[2];
	char		ny98omega3pr[100]; /* prior for class 3 omega (Ny98 model)         */
	MrBFlt		ny98omega3Fixed;
	MrBFlt		ny98omega3Uni[2];
	MrBFlt		ny98omega3Exp;
	char		m3omegapr[100];    /* prior for all three omegas (M3 model)        */
	MrBFlt		m3omegaFixed[3];
	char		m10betapr[100];    /* prior for omega variation (M10 model)        */
	char		m10gammapr[100];
	MrBFlt		m10betaExp;
	MrBFlt		m10betaUni[2];
	MrBFlt		m10betaFix[2];
	MrBFlt		m10gammaExp;
	MrBFlt		m10gammaUni[2];
	MrBFlt		m10gammaFix[2];
	char		codonCatFreqPr[100];/* prior for selection cat frequencies         */
	MrBFlt		codonCatFreqFix[3];
	MrBFlt		codonCatDir[3];
	char		stateFreqPr[100];  /* prior for character state frequencies        */
	MrBFlt		stateFreqsFix[200];
	MrBFlt		stateFreqsDir[200];
	char		stateFreqsFixType[100];
	int			numDirParams;
	char		shapePr[100];      /* prior for gamma shape parameter              */
	MrBFlt		shapeFix;
	MrBFlt		shapeUni[2];
	MrBFlt		shapeExp;
	char		pInvarPr[100];     /* prior for proportion of invariable sites     */
	MrBFlt		pInvarFix;
	MrBFlt		pInvarUni[2];
	char		adGammaCorPr[100]; /* prior for correlation param of adGamma model */
	MrBFlt		corrFix;
	MrBFlt		corrUni[2];
	char		covSwitchPr[100];  /* prior for switching rates of covarion model  */
	MrBFlt		covswitchFix[2];
	MrBFlt		covswitchUni[2];
	MrBFlt		covswitchExp;
	char		symPiPr[100];      /* prior for pi when unidentifiable states used */
	MrBFlt		symBetaFix;
	MrBFlt		symBetaUni[2];
	MrBFlt		symBetaExp;
	char		ratePr[100];       /* prior on rate for a partition                */
	MrBFlt		ratePrDir;			
	char		brownCorPr[100];   /* prior for correlation of Brownian model      */
	MrBFlt		brownCorrFix;
	MrBFlt		brownCorrUni[2];
	char		brownScalesPr[100];/* prior for scales of Brownian model           */
	MrBFlt		brownScalesFix;
	MrBFlt		brownScalesUni[2];
	MrBFlt		brownScalesGamma[2];
	MrBFlt		brownScalesGammaMean;

	char		topologyPr[100];   /* prior for tree topology                      */
	int			activeConstraints[30]; 
	int			numActiveConstraints;
	int			numActiveLocks;
	char		brlensPr[100];     /* prior on branch lengths                      */
	MrBFlt		brlensUni[2];
	MrBFlt		brlensExp;
	char		unconstrainedPr[100]; /* prior on branch if unconstrained          */
	char		clockPr[100];      /* prior on branch if clock enforced            */
	char		speciationPr[100]; /* prior on speciation rate                     */
	MrBFlt		speciationFix;
	MrBFlt		speciationUni[2];
	MrBFlt		speciationExp;
	char		extinctionPr[100]; /* prior on extinction rate                     */
	MrBFlt		extinctionFix;
	MrBFlt		extinctionUni[2];
	MrBFlt		extinctionExp;
	MrBFlt		sampleProb;        /* taxon sampling fraction (for b-d process)    */
	char		treeHeightPr[100]; /* prior on tree height for clock models        */
	MrBFlt		treeHeightGamma[2];
	MrBFlt		treeHeightExp;
	char		thetaPr[100];      /* prior on coalescence                         */
	MrBFlt		thetaFix;
	MrBFlt		thetaUni[2];
	MrBFlt		thetaExp;
	char		growthPr[100];      /* prior on coalescence growth rate            */
	MrBFlt		growthFix;
	MrBFlt		growthUni[2];
	MrBFlt		growthExp;
	MrBFlt		growthNorm[2];
	char		calWaitPr[100];    /* prior on clock waiting time                  */
	MrBFlt		calWaitFix;
	MrBFlt		calWaitUni[2];
	MrBFlt		calWaitExp;

	char		tratioFormat[30];      /* format used to report tratio				   */
	char		revmatFormat[30];      /* format used to report revmat				   */
	char		ratemultFormat[30];    /* format used to report ratemult     		   */
	char	    inferAncStates[5];     /* should ancestral states be inferred (Yes/No)?*/
	char	    inferSiteRates[5];     /* should site rates be inferred (Yes/No)?      */
	char	    inferPosSel[5];        /* should site selection be inferred (Yes/No)?  */
	} Model, ModelParams;

typedef struct chain
	{
	int			numGen;                /* number of MCMC cycles                         */
	int			sampleFreq;            /* frequency to sample chain                     */
	int			printFreq;             /* frequency to print chain                      */
	int			swapFreq;              /* frequency to attempt swap of states           */
	int			numRuns;               /* number of runs                                */
	int			numChains;             /* number of chains                              */
	MrBFlt		chainTemp;             /* chain temperature                             */
	int			userDefinedTemps;      /* should we use the users temperatures?         */
	MrBFlt		userTemps[MAX_CHAINS]; /* user-defined chain temperatures               */
	char		chainFileName[100];    /* chain file name for output                    */
	int			chainBurnIn;           /* chain burn in length                          */
	int			numStartPerts;         /* number of perturbations to starting tree      */
	char		chainStartTree[100];   /* starting tree for chain (random/user)         */
	int			saveBrlens;            /* should branch lengths be saved                */
	long int	chainSeed;             /* random seed for chain                         */
	MrBFlt		weightScheme[3];       /* percent chars to increase/decrease in weight  */
	int			calcPbf;               /* should we calculate the pseudo Bayes factor   */
	int			pbfInitBurnin;         /* initial burnin when calculating pseudo BF     */
	int			pbfSampleFreq;         /* sample frequency for pseudo BF                */
	int			pbfSampleTime;         /* how many cycles to calcualate site prob.      */
	int			pbfSampleBurnin;       /* burnin period for each site for pseudo BF     */
	int			swapAdjacentOnly;      /* whether we only swap adjacent temperatures    */
	int			redirect;              /* should output go to stdout                    */
	int			allChains;             /* should stats be output for all chains?        */
	int         numSwaps;              /* number of swaps to try each time              */
	int         mcmcDiagn;             /* should mcmc diagnostics be output?            */
	int         diagnFreq;             /* mcmc diagnostics frequency                    */
	int         relativeBurnin;        /* should a relative burnin be used ?            */
	MrBFlt      burninFraction;        /* the sample fraction to discard as burnin      */
	int         allComps;              /* top conv diagnosis for all pairs?             */
	MrBFlt      minPartFreq;           /* minimum partition frequency for conv diagn    */
	MrBFlt      stopVal;               /* top conv diagn value to reach before stopping */
	int		    stopRule;              /* use stop rule?                                */
	STATS		*stat;				   /* ptr to structs with mcmc diagnostics info     */
	Tree		*utree;				   /* pointing to utree used for conv diagnostics   */
	Tree		*rtree;				   /* pointing to rtree used for conv diagnostics   */
	int			printMax;			   /* maximum number of chains to print             */
	int			printAll;			   /* whether to print all or only cold chains      */
	int			runWithData;		   /* should we run with or without data?           */
	int			orderTaxa;		       /* order taxa before printing tree to file?      */
	char	mpiPath[100];				/* path for temp files in mpi version */
	} Chain;

typedef struct modelinfo
	{
	int			dataType;          			/* data type for partition                  */
	int			nucModelId;					/* structure of nucleotide model            */
	int			nst;						/* # substitution types                     */
	int 		aaModelId;					/* amino acid model type                    */
	int			parsModelId;				/* is parsimony model used YES/NO           */

	Param		*tRatio;					/* ptr to tRatio used in model				*/
	Param		*revMat;					/* ptr to revMat used in model				*/
	Param		*omega;						/* ptr to omega used in model				*/
	Param		*stateFreq;					/* ptr to statFreq used in model			*/
	Param		*shape;						/* ptr to shape used in model				*/
	Param		*pInvar;					/* ptr to pInvar used in model				*/
	Param		*correlation;				/* ptr to correlation used in model			*/
	Param		*switchRates;				/* ptr to switchRates (off<->on)            */
	Param		*rateMult;					/* ptr to rateMult used in model			*/
	Param		*speciationRates;			/* ptr to speciationRates used in model		*/
	Param		*extinctionRates;			/* ptr to extinctionRates used in model		*/
	Param		*theta;						/* ptr to theta used in model			 	*/
	Param		*growthRate;				/* ptr to theta used in model			 	*/
	Param		*topology;					/* ptr to topology used in model			*/
	Param		*brlens;					/* ptr to brlens (and tree) used in model	*/
	Param		*aaModel;					/* ptr to amino acid matrix used            */
	
	int			numChars;					/* number of compressed characters			*/
	int			numUncompressedChars;		/* number of uncompressed characters		*/
	int			numDummyChars;				/* number of dummy characters				*/
	int			compMatrixStart;			/* start column in compressed matrix		*/
	int			compMatrixStop;				/* stop column in compressed matrix 		*/
	int			compCharStart;				/* start char among compressed chars		*/
	int			compCharStop;				/* stop char among compressed chars			*/
	int			parsMatrixStart;			/* start column in parsimony matrix			*/
	int			parsMatrixStop;				/* stop collumn in parsimony matrix			*/
	int			nParsIntsPerSite;			/* # parsimony ints per character			*/	
	int			nCharsPerSite;				/* number chars per site (eg 3 for codon)	*/
	int			condLikeStart;				/* start column	in cond like matrix			*/
	int			condLikeStop;				/* stop column in cond like matrix			*/
	int			rateProbStart;				/* start of rate probs (for adgamma)		*/
	int			sprParsMatrixStart;			/* start column in SPR parsimony matrix		*/
	int			sprParsMatrixStop;			/* stop collumn in SPR parsimony matrix		*/
				
	int			numGammaCats;				/* number of gamma cats (1 if inapplic.)	*/
	int			numBetaCats;				/* number of beta cats (1 if inapplic.)	    */
	int			numOmegaCats;				/* number of omega cats	(1 if inapplic.)    */
	int			numTiCats;					/* number of cats needing different tis     */
	int			numModelStates;				/* number of states including hidden ones	*/
	int			numStates;					/* number of observable discrete states		*/

	int			tiProbStart;
	int			upDateCl;

	int			nCijk;						/* stores length of cijk vector */
	int			nCijkParts;					/* stores number of cijk partitions (for omega/covarion models) */
	int			cijkStart;
	int 		cijkBits[MAX_CHAINS][2];	/* bits indicate which cijk to use */
	int			upDateCijk[MAX_CHAINS][2];	/* whether cijk vector needs to be updated */

	int			*tiIndex;					/* index to trans probs for each compressed char */
	int			*bsIndex;					/* index to stat freqs for each compressed char  */
	int			*nStates;					/* # states of each compressed char */
	int			*cType;						/* whether char is ord, unord or irrev */
	int			*weight;					/* prior weight of each compressed char */
	int			isTiNeeded[20];				/* marks whether a trans prob matrix is needed */

	MrBFlt		parsTreeLength[MAX_CHAINS * 2];/* parsimony tree lengths for chains        */
	int			parsNodeLenStart;			/* start column in parsimony node length array */
	
	int			tiProbsId;
	int			*isPartAmbig;
	int			mark;
	MrBFlt		*invCondLikes;

	MrBFlt		lnLike[MAX_CHAINS * 2];

	LikeDownFxn		CondLikeDown;
	LikeRootFxn		CondLikeRoot;
	LikeScalerFxn	CondLikeScaler;
	LikeFxn			Likelihood;
	TiProbFxn		TiProbs;
	LikeUpFxn		CondLikeUp;
	PrintAncStFxn   PrintAncStates;
	StateCodeFxn	StateCode;
	PrintSiteRateFxn PrintSiteRates;

	int			printAncStates;				/* should ancestral states be printed (YES/NO)  */
	int			printSiteRates;             /* should site rates be printed (YES/NO)        */
	int			printPosSel;                /* should selection be printed (YES/NO)         */

	int         parsimonyBasedMove;         /* is parsimony-based move used (YES/NO)        */

	} ModelInfo;

typedef struct sumt
	{
	char		sumtFileName[100];     /* name of input file                            */
	int			sumtBurnIn;            /* burn in                                       */
	char		sumtConType[100];      /* consensus tree type                           */
	int			calcTrprobs;           /* should the individual tree probs be calculated*/
	int			showSumtTrees;         /* should the individual tree probs be shown     */
	MrBFlt		freqDisplay;           /* threshold for printing partitions to screen   */
	int			printBrlensToFile;     /* should branch lengths be printed to file      */
	MrBFlt		brlensFreqDisplay;     /* threshold for printing branch lengths to file */
	int			numRuns;			   /* number of independent analyses to summarize   */
	int			numTrees;              /* number of trees to summarize                  */
	int			orderTaxa;             /* order taxa in trees?                          */
	} Sumt;

typedef struct comptree
	{
	char		comptFileName1[100];    /* name of first input file                      */
	char		comptFileName2[100];    /* name of second input file                     */
	char		comptOutfile[100];      /* name of output file                           */
	int			comptBurnIn;            /* burn in                                       */
	} Comptree;

typedef struct sump
	{
	char		sumpFileName[100];     /* name of input file                            */
	int			sumpBurnIn;            /* burn in                                       */
	int			plot;                  /* output plot (y/n)?                            */
	int			table;                 /* output table (y/n)?                           */
	int			margLike;              /* output marginal likelihood (y/n)?             */
	int			printToFile;           /* print output to file (y/n)?                   */
	char	    sumpOutfile[100];      /* name of output file                           */
	int			numRuns;			   /* number of independent analyses to summarize   */
	int			allRuns;			   /* should data for all runs be printed (yes/no)? */
	int			overlayPlot;		   /* should plots from several runs be overlaid?   */
	} Sump;

typedef struct plot
	{
	char		plotFileName[100];     /* name of input file                            */
	char		parameter[100];        /* parameter(s) to be plotted                    */
	char		match[100];            /* whether the match needs to be perfect         */
	int			plotBurnIn;            /* burn in                                       */
	} Plot;

typedef struct
	{
	int		    numTrees;		       /* number of trees to reassemble                 */
	int		    numRuns;		       /* number of runs to reassemble                  */
	} ReassembleInfo;

typedef struct doublet
	{
	long int	first, second;
	} Doublet;

typedef struct matrix
	{
	long *origin;
	int rowSize;
	int nRows;
	int column;
	int row;
	} Matrix;

typedef struct charinfo
	{
	int dType;
	int cType;
	int nStates;
	int constant[10];
	int variable;
	int informative;
	} CharInfo;
	
typedef struct 
	{
	int			isExcluded;            /* is the character excluded                     */
	int			numStates;             /* number of observed states for the character   */
	int			charType;              /* type of character                             */
	int			isMissAmbig;           /* is the character missing or ambiguous         */
	int			ctype;                 /* ordering of character                         */
	int			charId;                /* char ID index for doublet and codon models    */
	int			pairsId;               /* char ID for doublets                          */
	int			bigBreakAfter;         /* is there a large break after this character   */
	int			charSet[MAX_NUM_CHARSETS]; /* holds defined character sets              */
	int			partitionId[MAX_NUM_PARTITIONS];/* the partitions                       */
	}
	CharInformation;

typedef struct 
	{
	int			isDeleted;             /* is the taxon deleted                          */
	Calibration	calibration;           /* the age of the taxon                          */
	int			charCount;             /* count holder                                  */
	int			taxaSet[MAX_NUM_TAXASETS]; /* holds defined taxon sets                  */
	int			constraints[MAX_NUM_CONSTRAINTS];/* the constraints                     */
	}
	TaxaInformation;

typedef struct
	{
	MrBFlt		curScore;
	MrBFlt		minScore;
	MrBFlt		totalScore;
	MrBFlt		stopScore;
	MrBFlt		warp;
	TreeNode	**leaf;
	TreeNode	**vertex;
	}
	TreeInfo;
