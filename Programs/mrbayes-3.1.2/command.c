/*
 *  MrBayes 3.1.2
 *
 *  copyright 2002-2005
 *
 *  John P. Huelsenbeck
 *  Section of Ecology, Behavior and Evolution
 *  Division of Biological Sciences
 *  University of California, San Diego
 *  La Jolla, CA 92093-0116
 *
 *  johnh@biomail.ucsd.edu
 *
 *	Fredrik Ronquist
 *  School of Computational Science
 *  Florida State University
 *  Tallahassee, FL 32306-4120
 *
 *  ronquist@csit.fsu.edu
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details (www.gnu.org).
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include "mb.h"
#include "globals.h"
#include "command.h"
#include "bayes.h"
#include "model.h"
#include "mcmc.h"
#include "plot.h"
#include "sump.h"
#include "sumt.h"
#if defined(__MWERKS__)
#include "SIOUX.h"
#endif


#	define	NUMCOMMANDS					53  /* Note: NUMCOMMANDS gives the total number  */
											/*       of commands in the program           */
#define	NUMPARAMS						200
#define PARAM(i, s, f, l)				p->string = s;    \
										p->fp = f;        \
										p->valueList = l; \
										p++;
#define	HIDE							0
#define	SHOW							1
#undef	SHOW_TOKENS

int      AddToSet (int i, int j, int k, int id);
int      AllocMatrix (void);
char     ChangeCase (char c);
int      CharacterCode (char ch, int *charCode, int chType);
int		 CharacterNumber (int charCode, int chType);
int      CheckInitialPartitions (void);
int      Dex (TreeNode *p);
int      DoAbout (void);
int      DoAcknowledgments (void);
int      DoBeginParm (char *parmName, char *tkn);
int      DoBreaks (void);
int      DoBreaksParm (char *parmName, char *tkn);
int      DoCalibrate (void);
int      DoCalibrateParm (char *parmName, char *tkn);
int      DoCharset (void);
int      DoCharsetParm (char *parmName, char *tkn);
int      DoCharStat (void);
int      DoCitations (void);
int      DoConstraints (void);
int      DoConstraintsParm (char *parmName, char *tkn);
int      DoCtype (void);
int      DoCtypeParm (char *parmName, char *tkn);
int      DoDelete (void);
int      DoDeleteParm (char *parmName, char *tkn);
int      DoDeroot (void);
int      DoDimensions (void);
int      DoDimensionsParm (char *parmName, char *tkn);
int      DoDisclaimer (void);
int      DoEndBlock (void);
int      DoExecute (void);
int      DoExecuteParm (char *parmName, char *tkn);
int      DoExclude (void);
int      DoExcludeParm (char *parmName, char *tkn);
int      DoFormat (void);
int      DoFormatParm (char *parmName, char *tkn);
int      DoHelp (void);
int      DoHelpParm (char *parmName, char *tkn);
int      DoInclude (void);
int      DoIncludeParm (char *parmName, char *tkn);
int      DoLog (void);
int      DoLogParm (char *parmName, char *tkn);
int      DoManual (void);
int      DoManualParm (char *parmName, char *tkn);
int      DoMatrix (void);
int      DoMatrixParm (char *parmName, char *tkn);
int      DoNexusParm (char *parmName, char *tkn);
int      DoOutgroup (void);
int      DoOutgroupParm (char *parmName, char *tkn);
int      DoPairs (void);
int      DoPairsParm (char *parmName, char *tkn);
int      DoPartition (void);
int      DoPartitionParm (char *parmName, char *tkn);
int      DoProps (void);
int      DoRestore (void);
int      DoRestoreParm (char *parmName, char *tkn);
int      DoRoot (void);
int      DoSet (void);
int      DoSetParm (char *parmName, char *tkn);
int      DoShowMatrix (void);
int      DoShowtree (void);
int      DoTaxaset (void);
int      DoTaxasetParm (char *parmName, char *tkn);
int      DoTaxaStat (void);
int      DoUserTree (void);
int      DoUserTreeParm (char *parmName, char *tkn);
int      DoVersion (void);
int      FindValidParam (char *tk, int *numMatches);
int      GetNumPartDivisions (int n);
int      GetUserHelp (char *helpTkn);
int		 IsAmbig (int charCode, int dType);
int      IsMissing (int charCode, int dType);
int      NBits (int x);
int      NucID (char nuc);
void     PrintYesNo (int yn, char s[4]);
int      ProtID (char aa);
int      RemoveLastFromString (char *s1);
int      MBResID (char nuc);
int      SetPartitionInfo (int n);
int      StandID (char nuc);
int      StateCode_AA (int n);
int      StateCode_NUC4 (int n);
int      StateCode_Std (int n);
void     WhatVariableExp (unsigned long int exp, char *st);
char     WhichAA (int x);
MrBFlt   WhichCont (int x);
char     WhichRes (int x);
char     WhichStand (int x);



/* globals */
int				autoClose;             /* autoclose                                     */
int 			autoOverwrite;                          /* Overwrite or append outputfiles when nowarnings=yes */
Calibration		*calibrationPtr;       /* ptr to calibration being set                  */
CharInformation *charInfo;             /* holds critical information about characters   */
char			*charSetNames;         /* holds names of character sets                 */
Comptree		comptreeParams;        /* holds parameters for comparetree command      */
Calibration		constraintCalibration[MAX_NUM_CONSTRAINTS];/* holds calibration of constraints */
char			*constraintNames;      /* holds names of constraints                    */
int				dataType;              /* type of data                                  */
int				echoMB;				   /* flag used by Manual to prevent echoing        */
unsigned long int expecting;           /* variable denoting expected token type         */
int				foundNewLine;          /* whether a new line has been found             */
int 			inComment;             /* flag for whether input stream is commented    */
int				inferAncStates;		   /* should ancestral states be inferred (y/n)     */
int				inferSiteRates;		   /* should site rates be inferred (y/n)           */
int				inMrbayesBlock;        /* flag for whether we are in a mrbayes block    */
int				inSumtBlock;           /* are we in the sumt block                      */
int				inValidCommand;        /* a useful flag set whenever you enter a cmd    */
int  			isInAmbig, isInPoly;   /* flags whether we are within () or {}          */
int				isTranslateDef;        /* is a translation block defined                */
int				isUserTreeDefined;     /* flag indicating whether user tree is found    */
char			logFileName[100];      /* name of the log file                          */
int				logToFile;             /* should screen output be logged to a file      */
FILE			*logFileFp;            /* file pointer to log file                      */
char			manFileName[100];      /* name of the file for the command help info    */
int				*matrix;               /* matrix containing original data               */
int  			matrixHasPoly;		   /* flag for whether matrix has polymorphisms     */
int				memAllocs[NUM_ALLOCS]; /* allocated memory flags                        */
int				mode;                  /* mode of program (interactive/noninteractive)  */
int				noWarn;                /* no warnings on overwriting files              */
int				numCharSets;           /* holds number of character sets                */
int				numComments;		   /* counts how deeply nested a comment is         */
int				numDefinedConstraints; /* number of constraints defined                 */
int				numDefinedPartitions;  /* number of partitions defined                  */
int				numOpenExeFiles;       /* number of execute files open                  */
int				numTaxaSets;           /* holds number of taxa sets                     */
int				outGroupNum;           /* number of outgroup taxon                      */
ParmInfo		paramTable[NUMPARAMS]; /* information on parameters                     */
char			*partitionNames;       /* hold names of partitions (first is "default") */
int				partitionNum;          /* number of current partition                   */
Plot			plotParams;            /* holds parameters for plot command             */
int				quitOnError;		   /* quit on error?					            */
MrBFlt			relConstraintProbs[MAX_NUM_CONSTRAINTS];/* rel. probs. of constraint    */
int				replaceLogFile;        /* should logfile be replace/appended to         */
char			spacer[10];            /* holds blanks for printing indentations        */
Sump			sumpParams;            /* holds parameters for sump command             */
Sumt			sumtParams;            /* holds parameters for sumt command             */
TaxaInformation *taxaInfo;             /* holds critical information about taxa         */
char			*taxaNames;            /* holds name of taxa                            */
int				*tempSet;              /* temporarily holds defined character set       */
char			*taxaSetNames;         /* holds names of taxa sets                      */
int  			theAmbigChar;          /* int containing ambiguous character            */
char			*transFrom;            /* translation block information                 */
char			*transTo;              /* translation block information                 */
int				userBrlensDef;         /* are the branch lengths on user tree defined   */
Tree			*userTree;             /* user tree                                     */
int				longIntegerSize;       /* size of an unsigned integer                   */

/* local (to this file) */
char			*tokenP, token[CMD_STRING_LENGTH];
CmdType			commands[] =
				{
			/*	Information on commands initialization:
			 
			 		1 = Command number (cmdNumber)
			 		2 = Command name (string)
			 		3 = Special command (YES/NO) (specialCmd) 
			 		4 = Pointer to finishing function (fp)
			 		5 = Number of valid parameters (numParms)
			 		6 = List of valid parameters (parmList) 
			 		7 = Expecting (2^TokenType) (expect) (PARAMETER = 4; SEMICOLON = 32; ALPHA = 16384; 
			 		    ALPHA | QUESTIONMARK | DASH | NUMBER | ASTERISK | EXCLAMATIONMARK | PERCENT | WEIRD | SEMICOLON = 11715360;
			 		    ALPHA | QUESTIONMARK | DASH | NUMBER | ASTERISK | EXCLAMATIONMARK | PERCENT | WEIRD | SEMICOLON | LEFTPAR | RIGHTPAR | LEFTCURL | RIGHTCURL = 112381728;
			 		    PARAMETER | SEMICOLON = 36; NUMBER | ALPHA = 49152; ALPHA | SEMICOLON = 16416; EQUALSIGN = 8; NUMBER = 32768)
			 		8 = Description of the command (cmdDescription)
			 		9 = Where should the command be used (cmdUse) (IN_CMD = used from command line or mrbayes block; IN_FILE = used in data block or in tree block)
			 	   10 = Should the command be shown when "help" is typed (hiding).
			 
			  #1                 #2   #3                 #4  #5                                                                                                #6        #7                                                            #8       #9   #10
			 -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- */
			{  0,               "#",  NO,              NULL,  1,                                                                                              {0},        4,                                                           "", IN_FILE, HIDE },
			{  1,           "About",  NO,           DoAbout,  0,                                                                                             {-1},       32,                                      "Describes the program",  IN_CMD, SHOW },
			{  2, "Acknowledgments",  NO, DoAcknowledgments,  0,                                                                                             {-1},       32,                              "Shows program acknowledgments",  IN_CMD, SHOW },
			{  3,           "Begin",  NO,              NULL,  3,                                                                                          {1,2,3},        4,                         "Denotes beginning of block in file", IN_FILE, SHOW },
			{  4,       "Calibrate",  NO,       DoCalibrate,  1,                                                                                            {120},        4,            "Assigns dates to terminals or constrained nodes",  IN_CMD, HIDE },
			{  5,         "Charset",  NO,         DoCharset,  1,                                                                                             {15},        4,                          "Assigns a group of sites to a set",  IN_CMD, SHOW },
			{  6,        "Charstat",  NO,        DoCharStat,  0,                                                                                             {-1},       32,                                 "Shows status of characters",  IN_CMD, SHOW },
			{  7,       "Citations",  NO,       DoCitations,  0,                                                                                             {-1},       32,                            "Appropriate citation of program",  IN_CMD, SHOW },
			{  8,     "Comparetree",  NO,     DoCompareTree,  4,																				{127,128,129,130},       36,                     "Compares the trees from two tree files",  IN_CMD, SHOW },
			{  9,      "Constraint",  NO,     DoConstraints,  1,                                                                                             {66},        4,                      "Defines a constraint on tree topology",  IN_CMD, SHOW },
			{ 10,           "Ctype",  NO,           DoCtype,  1,                                                                                             {65},        4,                        "Assigns ordering for the characters",  IN_CMD, SHOW },
			{ 11,      "Databreaks", YES,          DoBreaks,  1,                                                                                             {94},    32768,        "Defines nucleotide pairs (doublets) for stem models",  IN_CMD, SHOW },
			{ 12,          "Delete", YES,          DoDelete,  1,                                                                                             {47},    49152,                             "Deletes taxa from the analysis",  IN_CMD, SHOW },
			{ 13,          "Deroot",  NO,          DoDeroot,  0,                                                                                             {-1},       32,                                          "Deroots user tree",  IN_CMD, SHOW },
			{ 14,      "Dimensions",  NO,      DoDimensions,  2,                                                                                            {4,5},        4,                           "Defines size of character matrix", IN_FILE, SHOW },
			{ 15,      "Disclaimer",  NO,      DoDisclaimer,  0,                                                                                             {-1},       32,                               "Describes program disclaimer",  IN_CMD, SHOW },
			{ 16,             "End",  NO,        DoEndBlock,  0,                                                                                             {-1},       32,                             "Denotes end of a block in file", IN_FILE, SHOW },
			{ 17,        "Endblock",  NO,        DoEndBlock,  0,                                                                                             {-1},       32,                 "Alternative way of denoting end of a block", IN_FILE, SHOW },
			{ 18,         "Exclude", YES,         DoExclude,  1,                                                                                             {45},    49152,                           "Excludes sites from the analysis",  IN_CMD, SHOW },
			{ 19,         "Execute", YES,         DoExecute,  1,                                                                                             {12},    16384,                                            "Executes a file",  IN_CMD, SHOW },
			{ 20,          "Format",  NO,          DoFormat,  5,                                                                                     {6,7,8,9,10},        4,                     "Defines character format in data block", IN_FILE, SHOW },
			{ 21,            "Help", YES,            DoHelp,  1,                                                                                             {50},    16416,                  "Provides detailed description of commands",  IN_CMD, SHOW },
			{ 22,         "Include", YES,         DoInclude,  1,                                                                                             {46},    49152,                                             "Includes sites",  IN_CMD, SHOW },
			{ 23,            "Link",  NO,            DoLink, 16,                                              {55,56,57,58,59,60,61,62,63,72,73,74,75,76,106,119},        4,               "Links parameters across character partitions",  IN_CMD, SHOW },
			{ 24,             "Log",  NO,             DoLog,  5,                                                                                 {86,87,88,89,90},        4,                               "Logs screen output to a file",  IN_CMD, SHOW },
			{ 25,            "Lset",  NO,            DoLset, 14,                                                     {28,29,30,31,32,33,34,40,51,52,53,91,92,131},        4,                "Sets the parameters of the likelihood model",  IN_CMD, SHOW },
			{ 26,	       "Manual",  NO,          DoManual,  1,								    														{127},       36,				  "Prints a command reference to a text file",  IN_CMD, SHOW },
			{ 27,          "Matrix", YES,          DoMatrix,  1,                                                                                             {11},112381728,                 "Defines matrix of characters in data block", IN_FILE, SHOW },
			{ 28,            "Mcmc",  NO,            DoMcmc, 37,  {17,18,19,20,21,22,23,24,25,26,27,85,99,113,114,115,116,117,132,143,144,145,149,150,151,152,153,
																														  154,155,156,157,158,159,160,161,167,170},    36,                   "Starts Markov chain Monte Carlo analysis",  IN_CMD, SHOW },
			{ 29,           "Mcmcp",  NO,           DoMcmcp, 37,  {17,18,19,20,21,22,23,24,25,26,27,85,99,113,114,115,116,117,132,143,144,145,149,150,151,152,153,
																														  154,155,156,157,158,159,160,161,167,170},        4, "Sets the parameters of a chain (without starting analysis)",  IN_CMD, SHOW },
			{ 30,        "Outgroup", YES,        DoOutgroup,  1,                                                                                             {79},    49152,                                     "Changes outgroup taxon",  IN_CMD, SHOW },
			{ 31,           "Pairs", YES,           DoPairs,  1,                                                                                             {93},    32768,        "Defines nucleotide pairs (doublets) for stem models",  IN_CMD, SHOW },
			{ 32,       "Partition",  NO,       DoPartition,  1,                                                                                             {16},        4,                              "Assigns a character partition",  IN_CMD, SHOW },
			{ 33,            "Plot",  NO,            DoPlot,  4,                                                                                {107,108,109,110},       36,                        "Plots parameters from MCMC analysis",  IN_CMD, SHOW },
			{ 34,           "Prset",  NO,           DoPrset, 30,  {35,36,37,38,39,41,42,43,44,54,64,67,68,69,70,71,77,101,102,103,104,105,111,112,118,121,122,133,
																																						 134,169},        4,                         "Sets the priors for the parameters",  IN_CMD, SHOW },
			{ 35,           "Props",  NO,           DoProps,  0,                                                                                             {-1},       32,                                 "Set proposal probabilities",  IN_CMD, SHOW },
			{ 36,            "Quit",  NO,            DoQuit,  0,                                                                                             {-1},       32,                                          "Quits the program",  IN_CMD, SHOW },
			{ 37,          "Report",  NO,          DoReport,  7,															        {123,124,125,134,135,136,137},        4,                 "Controls how model parameters are reported",  IN_CMD, SHOW },
			{ 38,         "Restore", YES,         DoRestore,  1,                                                                                             {48},    49152,                                              "Restores taxa",  IN_CMD, SHOW },
			{ 39,            "Root",  NO,            DoRoot,  0,                                                                                             {-1},       32,                                            "Roots user tree",  IN_CMD, SHOW },
			{ 40,             "Set",  NO,             DoSet,  5,                                                                               {13,14,95,146,171},        4,      "Sets run conditions and defines active data partition",  IN_CMD, SHOW },
			{ 41,      "Showmatrix",  NO,      DoShowMatrix,  0,                                                                                             {-1},       32,                             "Shows current character matrix",  IN_CMD, SHOW },
			{ 42,       "Showmodel",  NO,       DoShowModel,  0,                                                                                             {-1},       32,                                       "Shows model settings",  IN_CMD, SHOW },
			{ 43,        "Showtree",  NO,        DoShowtree,  0,                                                                                             {-1},       32,                                            "Shows user tree",  IN_CMD, SHOW },
			{ 44,            "Sump",  NO,            DoSump,  9,                                                              {97,98,138,139,140,141,142,162,163},       36,                   "Summarizes parameters from MCMC analysis",  IN_CMD, SHOW },
			{ 45,            "Sumt",  NO,            DoSumt, 11,                                                        {81,82,83,96,100,147,148,164,165,166,168},       36,                        "Summarizes trees from MCMC analysis",  IN_CMD, SHOW },
			{ 46,        "Taxastat",  NO,        DoTaxaStat,  0,                                                                                             {-1},       32,                                       "Shows status of taxa",  IN_CMD, SHOW },
			{ 47,          "Taxset",  NO,         DoTaxaset,  1,                                                                                             {49},        4,                           "Assigns a group of taxa to a set",  IN_CMD, SHOW },
			{ 48,       "Translate", YES,       DoTranslate,  1,                                                                                             {84},    49152,                         "Defines alternative names for taxa", IN_FILE, SHOW },
			{ 49,            "Tree",  NO,            DoTree,  1,                                                                                             {80},        4,                          "Defines a tree from MCMC analysis", IN_FILE, SHOW },
			{ 50,          "Unlink",  NO,          DoUnlink, 16,                                              {55,56,57,58,59,60,61,62,63,72,73,74,75,76,106,119},        4,             "Unlinks parameters across character partitions",  IN_CMD, SHOW },
			{ 51,        "Usertree", YES,        DoUserTree,  1,                                                                                             {78},        8,                                 "Defines a single user tree",  IN_CMD, SHOW },
			{ 52,         "Version",  NO,         DoVersion,  0,                                                                                             {-1},       32,                                      "Shows program version",  IN_CMD, SHOW },
		    
		/* NOTE: If you add a command here, make certain to change NUMCOMMANDS (above, in this file) appropriately! */
		    { 999,             NULL,  NO,              NULL,  0,                                                                                             {-1},       32,                                                           "",  IN_CMD, HIDE }  
		};
int					inDataBlock, inForeignBlock, isInterleaved, isFirstMatrixRead, isFirstInterleavedBlock, 
					taxonCount, fromI, toJ, everyK, foundDash, foundSlash, foundFirst, isMixed, whichPartition,
					isNegative, numPartitions, charOrdering, foundExp, foundColon, isFirstNode, nextAvailableNode,
					pairId, firstPair;
char				gapId, missingId, matchId, tempSetName[100];
CmdType 			*commandPtr;
ParmInfoPtr			paramPtr;
TreeNode			*pPtr, *qPtr;




int AddToSet (int i, int j, int k, int id)

{

	int		m, n;
	
	if (id <= 0)
		{
		MrBayesPrint ("%s   The id for a temporary set should be greater than 0\n", spacer);
		return (ERROR);
		}
	
	if (i < 0 && j < 0)
		return (ERROR);
	else if (i < 0 && j >= 0)
		return (ERROR);
	else if (i >= 0 && j < 0)
		{
		if (k >= 0)
			return (ERROR);
		else
			{
			if (tempSet[i] != 0)
				{
				MrBayesPrint ("%s   Character %d defined more than once\n", spacer, i+1);
				return (ERROR);
				}
			tempSet[i] = id;
			}
		}
	else if (i >= 0 && j >= 0)
		{
		if (k < 0)
			{
			for (m=i; m<=j; m++)
				{
				if (tempSet[m] != 0)
					{
					MrBayesPrint ("%s   Character %d defined more than once\n", spacer, m+1);
					return (ERROR);
					}
				tempSet[m] = id;
				}
			}
		else
			{
			n = k;
			for (m=i; m<=j; m++)	
				{
				if (n % k == 0)
					{
					if (tempSet[m] != 0)
						{
						MrBayesPrint ("%s   Character %d defined more than once\n", spacer, m+1);
						return (ERROR);
						}
					tempSet[m] = id;
					}
				n++;
				}
			}
		}


	return (NO_ERROR);
	
}





int AddToString (char *s1, char *s2, int *x)

{

	int		i, j, startI=0, numPrev;
	
	/* s1 is the character string that should be appended to character string s2 
	   We assume that the end of each character string has a "|". */
	   
	i = numPrev = 0;
	while (s2[i] != '\0')
		{
		if (s2[i] == '|')
			numPrev++;
		i++;
		}
		
	if (numPrev == 0)
		startI = 0;
	else
		{
		i = j = 0;
		while (s2[i] != '\0')
			{
			if (s2[i] == '|')
				j++;
			i++;
			if (j == numPrev)
				{
				startI = i;
				break;
				}
			}
		}
		
	if (s2[startI] == '\0')
		{
		MrBayesPrint ("%s   String is too full (1 %d)\n", spacer, startI);
		return (ERROR);
		}
	
	i = startI;
	j = 0;
	while(s1[j] != '\0')
		{
		s2[i++] = s1[j++];
		if (s2[i] == '\0')
			{
			MrBayesPrint ("%s   String is too full (2 %d)\n", spacer, i);
			return (ERROR);
			}
		}
	s2[i] = '|';
	
	*x = numPrev + 1;

	return (NO_ERROR);
	
}





int AllocMatrix (void)

{

	int		i, j, tempSetSize;

	if (memAllocs[ALLOC_MATRIX] == YES)
		FreeMatrix();
	matrix = (int *)SafeMalloc((size_t) (numTaxa * numChar * sizeof(int)));
	if (!matrix)
		{
		MrBayesPrint ("%s   Problem allocating matrix (%d)\n", spacer, numTaxa * numChar * sizeof(int));
		goto errorExit;
		}
	for (i=0; i<numTaxa * numChar; i++)
		matrix[i] = 0;
	memAllocs[ALLOC_MATRIX] = YES;

	if (memAllocs[ALLOC_CHARINFO] == YES)
		goto errorExit;
	charInfo = (CharInformation *)SafeMalloc((size_t) (numChar * sizeof(CharInformation)));
	if (!charInfo)
		{
		MrBayesPrint ("%s   Problem allocating charInfo (%d)\n", spacer, numChar * sizeof(CharInformation));
		goto errorExit;
		}
	for (i=0; i<numChar; i++)
		{
		charInfo[i].isExcluded = NO;
		charInfo[i].numStates = 0;
		charInfo[i].charType = 0;
		charInfo[i].isMissAmbig = NO;
		charInfo[i].ctype = UNORD;
		charInfo[i].charId = 0;
		charInfo[i].pairsId = 0;
		charInfo[i].bigBreakAfter = NO;
		for (j=0; j<MAX_NUM_CHARSETS; j++)
			charInfo[i].charSet[j] = NO;
		for (j=0; j<MAX_NUM_PARTITIONS; j++)
			charInfo[i].partitionId[j] = NO;
		}
	numTaxaSets = 0;
	memAllocs[ALLOC_CHARINFO] = YES;

	if (memAllocs[ALLOC_TAXAINFO] == YES)
		goto errorExit;
	taxaInfo = (TaxaInformation *)SafeMalloc((size_t) (numTaxa * sizeof(TaxaInformation)));
	if (!taxaInfo)
		{
		MrBayesPrint ("%s   Problem allocating taxaInfo (%d)\n", spacer, numTaxa * sizeof(TaxaInformation));
		goto errorExit;
		}
	for (i=0; i<numTaxa; i++)
		{
		taxaInfo[i].isDeleted = NO;
		taxaInfo[i].calibration.prior = unconstrained; /* no calibration set */
		strcpy (taxaInfo[i].calibration.name, "Unconstrained");
		taxaInfo[i].calibration.age = -1.0;
		taxaInfo[i].calibration.lower = -1.0;
		taxaInfo[i].calibration.upper = -1.0;
		taxaInfo[i].calibration.offset = -1.0;
		taxaInfo[i].calibration.lambda = -1.0;
		taxaInfo[i].charCount = 0;
		for (j=0; j<MAX_NUM_TAXASETS; j++)
			taxaInfo[i].taxaSet[j] = 0;
		for (j=0; j<MAX_NUM_CONSTRAINTS; j++)
			taxaInfo[i].constraints[j] = 0;
		}
	for (i=0; i<MAX_NUM_CONSTRAINTS; i++)
		{
		relConstraintProbs[i] = 0.0;
		constraintCalibration[i].prior = unconstrained; /* no calibration set */
		strcpy (constraintCalibration[i].name, "Unconstrained");
		constraintCalibration[i].age = -1.0;
		constraintCalibration[i].lower = -1.0;
		constraintCalibration[i].upper = -1.0;
		constraintCalibration[i].offset = -1.0;
		constraintCalibration[i].lambda = -1.0;
		}
	numDefinedConstraints = 0;
	numTaxaSets = 0;
	memAllocs[ALLOC_TAXAINFO] = YES;

	if (memAllocs[ALLOC_PARTITIONNAMES] == YES)
		goto errorExit;
	partitionNames = (char *)SafeMalloc((size_t) (MAX_NUM_PARTITIONS * 100 * sizeof(char)));
	if (!partitionNames)
		{
		MrBayesPrint ("%s   Problem allocating partitionNames (%d)\n", spacer, MAX_NUM_PARTITIONS * 100 * sizeof(char));
		goto errorExit;
		}
	for (i=0; i<100*MAX_NUM_PARTITIONS; i++)
		{
		partitionNames[i] = ' ';
		if (i == (MAX_NUM_PARTITIONS*100 - 1))
			partitionNames[i] = '\0';
		}
	memAllocs[ALLOC_PARTITIONNAMES] = YES;

	if (memAllocs[ALLOC_CHARSETNAMES] == YES)
		goto errorExit;
	charSetNames = (char *)SafeMalloc((size_t) (MAX_NUM_CHARSETS * 100 * sizeof(char)));
	if (!charSetNames)
		{
		MrBayesPrint ("%s   Problem allocating charSetNames (%d)\n", spacer, MAX_NUM_CHARSETS * 100 * sizeof(char));
		goto errorExit;
		}
	for (i=0; i<MAX_NUM_CHARSETS*100; i++)
		{
		charSetNames[i] = ' ';
		if (i == (MAX_NUM_CHARSETS*100 - 1))
			charSetNames[i] = '\0';
		}
	memAllocs[ALLOC_CHARSETNAMES] = YES;

	if (memAllocs[ALLOC_TMPSET] == YES)
		goto errorExit;
	if (numChar > numTaxa)
		tempSetSize = numChar;
	else
		tempSetSize = numTaxa;
	tempSet = (int *)SafeMalloc((size_t) (tempSetSize * sizeof(int)));
	if (!tempSet)
		{
		MrBayesPrint ("%s   Problem allocating tempSet (%d)\n", spacer, tempSetSize * sizeof(int));
		goto errorExit;
		}
	for (i=0; i<tempSetSize; i++)
		tempSet[i] = 0;
	memAllocs[ALLOC_TMPSET] = YES;

	if (memAllocs[ALLOC_TAXANAMES] == YES)
		goto errorExit;
	taxaNames = (char *)SafeMalloc((size_t) (numTaxa * 100 * sizeof(char)));
	if (!taxaNames)
		{
		MrBayesPrint ("%s   Problem allocating taxaNames (%d)\n", spacer, numTaxa * 100 * sizeof(char));
		goto errorExit;
		}
	for (i=0; i<numTaxa*100; i++)
		{
		taxaNames[i] = ' ';
		if (i == numTaxa*100 - 1)
			taxaNames[i] = '\0';
		}
	memAllocs[ALLOC_TAXANAMES] = YES;

	if (memAllocs[ALLOC_TAXASETNAMES] == YES)
		goto errorExit;
	taxaSetNames = (char *)SafeMalloc((size_t) (MAX_NUM_DIVS * 100 * sizeof(char)));
	if (!taxaSetNames)
		{
		MrBayesPrint ("%s   Problem allocating taxaSetNames (%d)\n", spacer, 3000 * sizeof(char));
		goto errorExit;
		}
	for (i=0; i<3000; i++)
		{
		taxaSetNames[i] = ' ';
		if (i == (3000 - 1))
			taxaSetNames[i] = '\0';
		}
	memAllocs[ALLOC_TAXASETNAMES] = YES;

	if (memAllocs[ALLOC_CONSTRAINTNAMES] == YES)
		goto errorExit;
	constraintNames = (char *)SafeMalloc((size_t) (30 * 100 * sizeof(char)));
	if (!constraintNames)
		{
		MrBayesPrint ("%s   Problem allocating constraintNames (%d)\n", spacer, 3000 * sizeof(char));
		goto errorExit;
		}
	for (i=0; i<3000; i++)
		{
		constraintNames[i] = ' ';
		if (i == (3000 - 1))
			constraintNames[i] = '\0';
		}
	memAllocs[ALLOC_CONSTRAINTNAMES] = YES;
	
	if (memAllocs[ALLOC_USERTREE] == YES)
		goto errorExit;
	isUserTreeDefined = NO;
	memAllocs[ALLOC_USERTREE] = NO;
	
	if (memAllocs[ALLOC_TRANSFROM] == YES)
		goto errorExit;
	transFrom = (char *)SafeMalloc((size_t) (numTaxa * 100 * sizeof(char)));
	if (!transFrom)
		{
		MrBayesPrint ("%s   Problem allocating transFrom (%d)\n", spacer, numTaxa * 100 * sizeof(char));
		goto errorExit;
		}
	for (i=0; i<numTaxa*100; i++)
		{
		transFrom[i] = ' ';
		if (i == numTaxa*100 - 1)
			transFrom[i] = '\0';
		}
	memAllocs[ALLOC_TRANSFROM] = YES;

	if (memAllocs[ALLOC_TRANSTO] == YES)
		goto errorExit;
	transTo = (char *)SafeMalloc((size_t) (numTaxa * 100 * sizeof(char)));
	if (!transTo)
		{
		MrBayesPrint ("%s   Problem allocating transTo (%d)\n", spacer, numTaxa * 100 * sizeof(char));
		goto errorExit;
		}
	for (i=0; i<numTaxa*100; i++)
		{
		transTo[i] = ' ';
		if (i == numTaxa*100 - 1)
			transTo[i] = '\0';
		}
	memAllocs[ALLOC_TRANSTO] = YES;

	MrBayesPrint ("%s   Allocated matrix\n", spacer);
	return (NO_ERROR);
	errorExit:
		MrBayesPrint ("%s   Problem allocating matrix\n", spacer);
		FreeMatrix();
		return (ERROR);

}





char ChangeCase (char c)

{

	int		x;
	
	x = tolower(c);
	return (x);
		
}





int CharacterCode (char ch, int *charCode, int chType)

{

	
	if (chType == DNA || chType == RNA)
		{
		if ((*charCode = NucID (ch)) == -1)
			{
			MrBayesPrint ("%s   Unrecognized DNA/RNA character '%c'\n", spacer, ch);
			return (ERROR);
			}
		}
	else if (chType == PROTEIN)
		{
		if ((*charCode = ProtID (ch)) == -1)
			{
			MrBayesPrint ("%s   Unrecognized Protein character '%c'\n", spacer, ch);
			return (ERROR);
			}
		}
	else if (chType == RESTRICTION)
		{
		if ((*charCode = MBResID (ch)) == -1)
			{
			MrBayesPrint ("%s   Unrecognized Restriction character '%c'\n", spacer, ch);
			return (ERROR);
			}
		}
	else if (chType == STANDARD)
		{
		if ((*charCode = StandID (ch)) == -1)
			{
			MrBayesPrint ("%s   Unrecognized Standard character '%c'\n", spacer, ch);
			return (ERROR);
			}
		}
	else if (chType == CONTINUOUS)
		{
		MrBayesPrint ("%s   CharacterCode function cannot check continuous characters\n", spacer);
		}
	else
		{
		MrBayesPrint ("%s   Unrecognized character type (%d)\n", spacer, chType);
		return (ERROR);
		}
		
	return (NO_ERROR);
	
}





int CharacterNumber (int charCode, int chType)

{

	int i, x = charCode;
	
	if (chType == CONTINUOUS)
		return 0;

	for (i=0; x!=0; i++)
		x >>= 1;

	return (i);

}





int CheckInitialPartitions (void)

{

	int		i;
	
	for (i=0; i<numChar; i++)
		{
		if (charInfo[i].partitionId[0] <= 0 || charInfo[i].partitionId[0] > numPartitions)
			{
			MrBayesPrint ("%s   The partition for site %d is incorrect\n", spacer, i+1); 
			return (ERROR);
			}
		}
		
	return (NO_ERROR);
	
}





int CheckStringValidity (char *s)

{

	int			i, numUnknownChars, tempNumComments, tempInComment;
	char		temp[100];

	i = 0;
	numUnknownChars = 0;
	tempNumComments = numComments;
	tempInComment = inComment;

	while(s[i] != '\0')
		{
		if (tempInComment == NO)
			{
			if ( !IsIn(s[i],"=abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ_0123456789.;:,#()[]?-*/'\\'!%\"&~+^$@|{}`>< ") )
				{
				if (IsWhite(s[i]) == 1 || IsWhite(s[i]) == 2)
					{
					
					}
				else
					{
					if ( commandPtr == NULL) 
						return (ERROR);
					MrBayesPrint ("%s   Unknown character \"%c\" (ASCII code %d)\n", spacer, s[i], s[i]);
					if (!strcmp(commandPtr->string,"Matrix"))
						{
						if (foundNewLine == NO)
							{
							MrBayesPrint ("%s   The error is in character %d for taxon %s\n", spacer, taxaInfo[taxonCount-1].charCount+i+1);
							}
						else
							{
							if (taxonCount == 0)
								MrBayesPrint ("%s   The error is in the first taxon name\n", spacer);
							else
								{
								GetNameFromString (taxaNames, temp, taxonCount);
								if (isInterleaved == NO)
									MrBayesPrint ("%s   The error is in the name of the taxon following taxon %s\n", spacer, temp);
								else
									{
									MrBayesPrint ("%s   The error is in the name of the taxon following taxon %s\n", spacer, temp);
									MrBayesPrint ("%s   in one of the interleaved data blocks\n", spacer);
									}
								}
							}
						}
					else if (!strcmp(commandPtr->string,"Execute"))
						{
						MrBayesPrint ("%s   Assuming irrelevant characters at beginning of file; processing continues\n", spacer);
						return (NO_ERROR);
						}
					return (ERROR);
					}
				}
			if (s[i]=='[')
				{
				tempInComment = YES;
				tempNumComments++;
				}
			}
		else if (tempInComment == YES)
			{
			if (s[i]==']')
				{
				tempNumComments--;
				if (tempNumComments == 0)
					tempInComment = NO;
				}
			}
		i++;
		}
		
	if (numUnknownChars > 0)
		return (ERROR);
	else
		return (NO_ERROR);

}





int CheckString (char *s1, char *s2, int *x)

{

	int			i, j, len1, len2, numOfStr, foundString, isIdentical, nDiff;
	char		tempStr[100];
		
	/* This function checks a string s1 against a list of names contained in s2. The
	   list of names in s2 is formatted like:
	   
	      "name1|name2|name3|                         \0"
	      
	   This string contains only three names and the names are separated by a "|".
	   When we go through the list of names, we know that we have found the end
	   of a name when we find "|" and that we have found the end of the list of
	   names when we find "\0" or " " (a blank). The command returns an ERROR if
	   a match was not found. Otherwise, a NO_ERROR is returned and *x is the
	   number in the list (s2) that matched s1. The number, x, is indexed 1, 2,
	   3, ... */
	
	len1 = (int) strlen(s1);
	i = j = numOfStr = 0;
	foundString = NO;
	while (s2[i] != '\0' && s2[i] != ' ')
		{
		if (j >= 100 - 2)
			{
			return (ERROR);
			}
		if (s2[i] != '|' && s2[i] != '\0' && s2[i] != ' ' && j < 100 - 2)
			{
			tempStr[j++] = s2[i];
			}
		else
			{
			tempStr[j++] = '\0';
			len2 = (int) strlen(tempStr);
			isIdentical = YES;
			if (len1 != len2)
				{
				isIdentical = NO;
				}
			else
				{
				nDiff = 0;
				for (j=0; j<len1; j++)
					if (tolower(s1[j]) != tolower(tempStr[j]))
						nDiff++;
				if (nDiff > 0)
					isIdentical = NO;
				}
			if (isIdentical == YES)
				{
				foundString = YES;
				*x = numOfStr + 1;
				break;
				}
			j = 0;
			numOfStr++;
			}
		i++;
		}
		
	if (foundString == NO)
		return (ERROR);
			
	return (NO_ERROR);
	
}





int Dex (TreeNode *p)

{

	return (p == NULL) ? -1 : p->index;

}





int DoAbout (void)

{

	MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
	MrBayesPrint ("   About the program                                                             \n");
    MrBayesPrint ("                                                                                 \n");
	MrBayesPrint ("   MrBayes is a program for the Bayesian estimation of phylogeny. Bayesian       \n");
	MrBayesPrint ("   inference of phylogeny is based upon the posterior probability distribution   \n");
	MrBayesPrint ("   of trees. Trees are labelled T1, T2, ..., Tn, where n is the number of        \n");
	MrBayesPrint ("   possible trees. The posterior probability of the i-th tree is calculated      \n");
	MrBayesPrint ("   using Bayes\'s formula as                                                     \n");
    MrBayesPrint ("                                                                                 \n");
	MrBayesPrint ("      Pr[Ti | X] = Pr[X | Ti] X Pr[Ti] / Pr[X]                                   \n");
    MrBayesPrint ("                                                                                 \n");
	MrBayesPrint ("   where X is a character matrix. Here, \"Pr[Ti | X]\" is the posterior          \n");
	MrBayesPrint ("   probability of the i-th tree, \"Pr[X | Ti]\" is the likelihood of the         \n");
	MrBayesPrint ("   i-th tree, and \"Pr[Ti]\" is the prior probability of the i-th tree. The      \n");
	MrBayesPrint ("   denominator of Bayes\'s formula (\"Pr[X]\") is a normalizing constant that    \n");
	MrBayesPrint ("   involves a summation over all possible trees. The likelihood, as described    \n");
	MrBayesPrint ("   above, cannot be calculated with knowledge of only the tree\'s topology. You  \n");
	MrBayesPrint ("   also need to have information on the lenths of the branches and on the        \n");
	MrBayesPrint ("   mechanism of character change. Hence, the likelihood (\"Pr[X | Ti]\")         \n");
	MrBayesPrint ("   involves a multidimensional integral over all possible combinations of        \n");
	MrBayesPrint ("   branch lengths and substitution model parameters.                             \n");
    MrBayesPrint ("                                                                                 \n");
	MrBayesPrint ("   In practice, it is impossible to calculate the posterior probability dist-    \n");
	MrBayesPrint ("   ribution of trees analytically. Instead, the posterior probability            \n");
	MrBayesPrint ("   of trees must be approximated. MrBayes uses a method called Markov chain      \n");
	MrBayesPrint ("   Monte Carlo (MCMC) to approximate the posterior probability of trees.         \n");
	MrBayesPrint ("   The object of MCMC is to construct a Markov chain that has as its state       \n");
	MrBayesPrint ("   space the parameters of the phylogenetic model and a stationary distribution  \n");
	MrBayesPrint ("   that is the posterior probability distribution of trees. MCMC takes valid,    \n");
	MrBayesPrint ("   albeit dependent, samples from the posterior probability distribution of      \n");
	MrBayesPrint ("   trees. The fraction of the time any tree appears in this sample is a          \n");
	MrBayesPrint ("   valid approximation of the posterior probability of the tree. MrBayes keeps   \n");
	MrBayesPrint ("   track of all the parameters of the phylogenetic model. The trees (with branch \n");
	MrBayesPrint ("   lengths) that were sampled by the MCMC procedure are saved in one file        \n");
	MrBayesPrint ("   (a file with a \".t\" extension) whereas the parameters of the model of       \n");
	MrBayesPrint ("   character change are saved in another file (a file with a \".p\" ext-         \n");
	MrBayesPrint ("   ension). You can summarize the results in the \".t\" and \".p\" files         \n");
	MrBayesPrint ("   using the \"sumt\" and \"sump\" commands, respectively.                       \n");
    MrBayesPrint ("                                                                                 \n");
	MrBayesPrint ("   MrBayes is cowritten by John Huelsenbeck and Fredrik Ronquist. It is          \n");
	MrBayesPrint ("   rather unusual to find a program of this sort that has multiple authors.      \n");
	MrBayesPrint ("   However, each author brings unique strengths to the development of the        \n");
	MrBayesPrint ("   program. Originally, the program was started by JH in August of               \n");
	MrBayesPrint ("   2000 and was intended to be distributed to a small number of people.          \n");
	MrBayesPrint ("   In March of 2001, Fredrik started making contributions to the program.        \n");
	MrBayesPrint ("   The contributions were of such a significant nature that he was made          \n");
	MrBayesPrint ("   a coauthor of the program. In particular, FR improved the speed of the        \n");
	MrBayesPrint ("   likelihood functions of the program and included new proposal mechanisms      \n");
	MrBayesPrint ("   for changing trees. The newest version of the program, v%s, is a              \n", VERSION_NUMBER);
	MrBayesPrint ("   completely rewritten version of the earlier program, and has more             \n");
	MrBayesPrint ("   integrally included FR\'s contributions.	                                    \n");
	MrBayesPrint ("   ---------------------------------------------------------------------------   \n");

	return (NO_ERROR);
	
}





int DoAcknowledgments (void)

{

	MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
	MrBayesPrint ("   Acknowledgments                                                               \n");
    MrBayesPrint ("                                                                                 \n");
	MrBayesPrint ("   JPH and FR would like to thank Gautam Altekar, Andrea Betancourt, Jon         \n");
    MrBayesPrint ("   Bollback, Barry Hall, Jimmy McGuire, Rasmus Nielsen, David Swofford,          \n");
    MrBayesPrint ("   Johan Nylander, Mikael Thollesson, and Derrick Zwickl for help during the     \n");
    MrBayesPrint ("   development of this program. Gautam Altekar, especially, was instrumental     \n");
    MrBayesPrint ("   in getting the parallel version of the program working.                       \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("   Our wives -- Edna Huelsenbeck and Eva Ronquist -- showed extraordinary        \n");
    MrBayesPrint ("   patience with us while we spent many late nights programming.                 \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("   JPH was supported by NSF grants DEB-007540 and MCB-0075404 and a Wenner-      \n");
    MrBayesPrint ("   Gren scholarship while writing this program. FR was supported by grants       \n");
    MrBayesPrint ("   from the Swedish Natural Science Research Council and the Swedish Research    \n");
    MrBayesPrint ("   Council.                                                                      \n");
	MrBayesPrint ("   ---------------------------------------------------------------------------   \n");

	return (NO_ERROR);
	
}





int DoBeginParm (char *parmName, char *tkn)

{
		
	if (expecting == Expecting(PARAMETER))
		{
		/* set Data (inDataBlock) *************************************************************/
		if (!strcmp(parmName, "Data"))
			{
			if (FreeMatrix () == ERROR)
				return (ERROR);
			MrBayesPrint ("   Reading data block\n");
			inDataBlock = YES;
			expecting = Expecting(SEMICOLON);
			strcpy (spacer, "   ");
			}
		/* set Mrbayes (inMrbayesBlock) *******************************************************/
		else if (!strcmp(parmName, "Mrbayes"))
			{
			MrBayesPrint ("   Reading MrBayes block\n");
			inMrbayesBlock = YES;
			expecting = Expecting(SEMICOLON);
			strcpy (spacer, "   ");
			}
		/* set Foreign (inForeignBlock) *******************************************************/
		else
			{
			MrBayesPrint ("   Skipping \"%s\" block\n", tkn);
			inForeignBlock = YES;
			expecting = Expecting(SEMICOLON);
			strcpy (spacer, "");
			}
		}
	else
		return (ERROR);

	return (NO_ERROR);

}





int DoBreaks (void)

{

	int			i, numBreaks;
	
	numBreaks = 0;
	for (i=0; i<numChar; i++)
		{
		if (charInfo[i].bigBreakAfter == YES)
			{
			numBreaks++;
			}
		}
		
	if (numBreaks > 0)
		{
		if (numBreaks == 1)
			MrBayesPrint ("%s   One data break found after character ", spacer, numBreaks);
		else
			MrBayesPrint ("%s   %d data breaks found after characters: ", spacer, numBreaks);
		for (i=0; i<numChar; i++)
			{
			if (charInfo[i].bigBreakAfter == YES)
				{
				MrBayesPrint ("%d ", i+1);
				}
			}
		MrBayesPrint ("\n");

		if (numBreaks == 1)
			MrBayesPrint ("%s   Successfully defined one break in data\n", spacer);
		else
			MrBayesPrint ("%s   Successfully defined %d breaks in data\n", spacer, numBreaks);
		}
	else
		{
		MrBayesPrint ("%s   No breaks in data found\n", spacer);
		}
		
	return (NO_ERROR);
	
}





int DoBreaksParm (char *parmName, char *tkn)

{

	int		i, tempInt;
		
	if (defMatrix == NO)
		{
		MrBayesPrint ("%s   A matrix must be specified before you can define breaks in the data\n", spacer);
		return (ERROR);
		}
			
	if (expecting == Expecting(NUMBER))
		{
		sscanf (tkn, "%d", &tempInt);
		if (tempInt <= 0 || tempInt > numChar)
			{
			MrBayesPrint ("%s   Character number %d is out of range (should be between %d and %d)\n", spacer, tempInt, 1, numChar);
			for (i=0; i<numChar; i++)
				charInfo[i].bigBreakAfter = NO;
			return (ERROR);
			}
		if (tempInt == numChar)
			{
			MrBayesPrint ("%s   Character number %d is the last character. MrBayes will define the\n", spacer, tempInt);
			MrBayesPrint ("%s   break, even though it doesn't make too much sense.\n", spacer);
			}
		tempInt--;
					
		charInfo[tempInt].bigBreakAfter = YES;
		
		expecting  = (Expecting(NUMBER) | Expecting(SEMICOLON));
		}
	else
		{
		for (i=0; i<numChar; i++)
			charInfo[i].bigBreakAfter = NO;
		return (ERROR);
		}

	return (NO_ERROR);
	MrBayesPrint ("%s", parmName); /* just because I am tired of seeing the unused parameter error msg */

}





int DoCalibrate (void)

{

#	if 0
	int			i;
#	endif

	/* show calibration times (for debugging) */
#	if 0
	MrBayesPrint ("Taxon ages\n");
	for (i=0; i<numTaxa; i++)
		MrBayesPrint ("%4d  --  %s\n", i+1, taxaInfo[i].calibration.name);
	MrBayesPrint ("Constraint ages\n");
	for (i=0; i<30; i++)
		MrBayesPrint ("%4d  --  %s\n", i+1, constraintCalibration[i].name);
#	endif

	return (NO_ERROR);

}





int DoCalibrateParm (char *parmName, char *tkn)

{

	static int		isTaxon;
	static char		nodeName[100];
	int				howMany, index;
	char			s[20], tempStr[100];
	MrBFlt			tempD;
		
	if (defMatrix == NO)
		{
		MrBayesPrint ("%s   A matrix must be specified before you can calibrate nodes\n", spacer);
		return (ERROR);
		}
		
	if (expecting == Expecting(PARAMETER))
		{
		if (strcmp(parmName, "Xxxxxxxxxx") != 0)
			{
			MrBayesPrint ("%s   Unexpected error - Wrong parmName in DoCalibrateParm\n", spacer);
			return (ERROR);
			}

		/* find taxon with this name */
		calibrationPtr = NULL;
		howMany = 0;

		/* first look in constraint names */
		if (CheckString (tkn, constraintNames, &index) != ERROR)
			{
			calibrationPtr = &constraintCalibration[index-1];
			howMany++;
			isTaxon = NO;
			strcpy (nodeName, tkn);
			}
		
		/* then look in terminal taxon names */
		if (CheckString (tkn, taxaNames, &index) != ERROR)
			{
			calibrationPtr = &taxaInfo[index-1].calibration;
			howMany++;
			isTaxon = YES;
			strcpy (nodeName, tkn);
			}

		/* return error if not found or ambiguous */
		if (howMany == 0)
			{
			MrBayesPrint ("%s   No taxon or constraint named ""%s"" found\n", spacer, tkn);
			return (ERROR);
			}
		else if (howMany > 1)
			{
			MrBayesPrint ("%s   Both a taxon and a constraint named ""%s"" encountered -- please rename one\n", spacer, tkn);
			return (ERROR);
			}

		if (calibrationPtr->prior != unconstrained)
			{
			MrBayesPrint ("%s   Resetting previous calibration for ""%s""\n", spacer, tkn);
			}

		/* reset the values of the calibration */
		strcpy (calibrationPtr->name, "Unconstrained");
		calibrationPtr->prior = unconstrained;
		calibrationPtr->lower = -1.0;
		calibrationPtr->upper = -1.0;
		calibrationPtr->age = -1.0;
		calibrationPtr->offset = -1.0;
		calibrationPtr->lambda = -1.0;

		/* get ready to find the equal sign */
		expecting = Expecting(EQUALSIGN);
		}

	else if (expecting == Expecting(EQUALSIGN))
		{
		/* get ready to find the calibration prior */
		expecting = Expecting(ALPHA);
		}

	else if (expecting == Expecting(ALPHA))
		{
		/* set the calibration prior type */
		if (IsArgValid(tkn,tempStr) == NO_ERROR)
			{
			if (!strcmp (tempStr, "Uniform"))
				calibrationPtr->prior = uniform;
			else if (!strcmp (tempStr, "Offsetexponential"))
				calibrationPtr->prior = offsetExponential;
			else if (!strcmp (tempStr, "Fixed"))
				calibrationPtr->prior = fixed;
			else /* if (!strcmp (tempStr, "Unconstrained")) */
				calibrationPtr->prior = unconstrained;
			strcpy (calibrationPtr->name, tempStr);
			if (calibrationPtr->prior == unconstrained)
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
			else
				expecting = Expecting(LEFTPAR);
			}
		else
			{
			MrBayesPrint ("%s   Invalid calibration prior argument \n", spacer);
			return (ERROR);
			}
		}
	else if (expecting == Expecting(LEFTPAR))
		{
		strcat (calibrationPtr->name, "(");
		expecting  = Expecting(NUMBER);
		}
	else if (expecting == Expecting(NUMBER))
		{
		if (calibrationPtr->prior == fixed)
			{
			sscanf (tkn, "%lf", &tempD);
			if (tempD <= 0.0)
				{
				MrBayesPrint ("%s   Age must be positive\n", spacer);
				calibrationPtr->prior = unconstrained;
				return (ERROR);
				}
			calibrationPtr->age = tempD;
			expecting = Expecting(RIGHTPAR);
			}
		else if (calibrationPtr->prior == uniform)
			{
			sscanf (tkn, "%lf", &tempD);
			if (tempD <= 0.0)
				{
				MrBayesPrint ("%s   Age must be positive\n", spacer);
				calibrationPtr->prior = unconstrained;
				return (ERROR);
				}
			if (calibrationPtr->lower < 0.0)
				{
				calibrationPtr->lower = tempD;
				expecting = Expecting(COMMA);
				}
			else
				{
				if (tempD <= calibrationPtr->lower)
					{
					MrBayesPrint ("%s   Maximum age must be larger than minimum age\n", spacer);
					calibrationPtr->prior = unconstrained;
					return (ERROR);
					}
				calibrationPtr->upper = tempD;
				expecting = Expecting(RIGHTPAR);
				}
			}
		else if (calibrationPtr->prior == offsetExponential)
			{
			sscanf (tkn, "%lf", &tempD);
			if (calibrationPtr->offset < 0.0)
				{
				if (tempD <= 0.0)
					{
					MrBayesPrint ("%s   Offset age must be positive\n", spacer);
					calibrationPtr->prior = unconstrained;
					return (ERROR);
					}
				calibrationPtr->offset = tempD;
				expecting = Expecting(COMMA);
				}
			else
				{
				if (tempD <= MIN_OFFSET_EXP_LAMBDA)
					{
					MrBayesPrint ("%s   Offset exponential lambda parameter must be larger than %lf\n", spacer, MIN_OFFSET_EXP_LAMBDA);
					calibrationPtr->prior = unconstrained;
					return (ERROR);
					}
				else if (tempD >= MAX_OFFSET_EXP_LAMBDA)
					{
					MrBayesPrint ("%s   Offset exponential lambda parameter must be smaller than %lf\n", spacer, MAX_OFFSET_EXP_LAMBDA);
					calibrationPtr->prior = unconstrained;
					return (ERROR);
					}
				calibrationPtr->lambda = tempD;
				expecting = Expecting(RIGHTPAR);
				}
			}
		sprintf (s, "%1.2lf", tempD);
		strcat (calibrationPtr->name, s);
		}
	else if (expecting == Expecting(COMMA))
		{
		strcat (calibrationPtr->name, ",");
		expecting  = Expecting(NUMBER);
		}
	else if (expecting == Expecting(RIGHTPAR))
		{
		strcat (calibrationPtr->name, ")");
		if (isTaxon == YES)
			MrBayesPrint ("%s   Setting age of taxon %s to %s\n", spacer, nodeName, calibrationPtr->name);
		else
			MrBayesPrint ("%s   Setting age of constraint node %s to %s\n", spacer, nodeName, calibrationPtr->name);
		/* get ready to find more calibrated nodes or taxa, if present */
		expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
		}
	else
		return (ERROR);

	return (NO_ERROR);
	
}




int DoCharset (void)

{

	int			i, howMany;

	/* first add name to charSetName */
	if (AddToString (tempSetName, charSetNames, &howMany) == ERROR)
		{
		MrBayesPrint ("%s   Problem adding charset %s to list\n", spacer, tempSetName);
		return (ERROR);
		}
	if (howMany != numCharSets + 1)
		{
		MrBayesPrint ("%s   Problem adding charset %s to list\n", spacer, tempSetName);
		if (RemoveLastFromString (charSetNames) == ERROR)
			return (ERROR);
		return (ERROR);
		}

	/* add set to tempSet */
	if (fromI >= 0 && toJ < 0)
		{
		if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
			{
			if (RemoveLastFromString (charSetNames) == ERROR)
				return (ERROR);
			return (ERROR);
			}
		}
	else if (fromI >= 0 && toJ >= 0 && everyK < 0)
		{
		if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
			{
			if (RemoveLastFromString (charSetNames) == ERROR)
				return (ERROR);
			return (ERROR);
			}
		}
	else if (fromI >= 0 && toJ >= 0 && everyK >= 0)
		{
		if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
			{
			if (RemoveLastFromString (charSetNames) == ERROR)
				return (ERROR);
			return (ERROR);
			}
		}
		
	/* merge tempSet with charset */
	for (i=0; i<numChar; i++)
		{
		if (tempSet[i] != 0)
			{
			charInfo[i].charSet[numCharSets] = tempSet[i];
			}
		}
	
	/* increment number of char sets */
	numCharSets++;

	return (NO_ERROR);
	
}





int DoCharsetParm (char *parmName, char *tkn)

{

	int		i, howMany, tempInt, allDigit;
	
	if (defMatrix == NO)
		{
		MrBayesPrint ("%s   A matrix must be specified before charsets can be defined\n", spacer);
		return (ERROR);
		}

	if (expecting == Expecting(PARAMETER))
		{
		if (!strcmp(parmName, "Xxxxxxxxxx"))
			{
			/* check that the name of the charset is not a number */
			allDigit = YES;
			for (i=0; i<(int)strlen(tkn); i++)
				{
				if (tkn[i] == '0' || tkn[i] == '1' || tkn[i] == '2' || tkn[i] == '3' || tkn[i] == '4' || 
				    tkn[i] == '5' || tkn[i] == '6' || tkn[i] == '7' || tkn[i] == '8' || tkn[i] == '9' || tkn[i] == '.')
				    {}
				else
					allDigit = NO;
				}
			if (allDigit == YES)
				{
				MrBayesPrint ("%s   Charset name may not be a number\n", spacer);
				return (ERROR);
				}
			
			/* check size of charset name */
			if (strlen(tkn) > 99)
				{
				MrBayesPrint ("%s   Charset name is too long\n", spacer);
				return (ERROR);
				}
				
			/* check to see if the name has already been used as a charset */
			if (numCharSets > 1)
				{
				if (CheckString (tkn, charSetNames, &howMany) == ERROR)
					{
					/* if the charset name has not been used, then we should have an ERROR returned */
					/* we _want_ to be here */

					}
				else
					{
					MrBayesPrint ("%s   Charset name has been used previously\n", spacer);
					return (ERROR);
					}
				}
			else if (numCharSets > MAX_NUM_CHARSETS)
				{
				MrBayesPrint ("%s   You cannot define more than %d charsets\n", spacer, MAX_NUM_CHARSETS);
				return (ERROR);
				}
				
			/* add the name to the character set */
			strcpy (tempSetName, tkn);
			
			/* clear tempSet */
			for (i=0; i<numChar; i++)
				tempSet[i] = 0;
			
			fromI = toJ = everyK = -1;
			foundDash = foundSlash = NO;
			MrBayesPrint ("%s   Defining charset called %s\n", spacer, tkn);
			expecting = Expecting(EQUALSIGN);
			}
		else
			return (ERROR);
		}
	else if (expecting == Expecting(EQUALSIGN))
		{
		expecting  = Expecting(ALPHA);
		expecting |= Expecting(NUMBER);
		}
	else if (expecting == Expecting(ALPHA))
		{
		/* We are defining a character set in terms of another (called tkn, here). We should be able
		   to find tkn in the list of character set names. If we cannot, then we have a problem and
		   return an error. */
		if (numCharSets < 1)
			{
			MrBayesPrint ("%s   Could not find a character set called %s\n", spacer, tkn);
			return (ERROR);
			}
		if (CheckString (tkn, charSetNames, &howMany) == ERROR)
			{
			MrBayesPrint ("%s   Could not find a character set called %s\n", spacer, tkn);
			return (ERROR);
			}
		/* add characters from charset "tkn" to new tempset */
		for (i=0; i<numChar; i++)
			{
			if (charInfo[i].charSet[howMany-1] == 1)
				tempSet[i] = 1;
			}		
		fromI = toJ = everyK = -1;

		expecting  = Expecting(ALPHA);
		expecting |= Expecting(NUMBER);
		expecting |= Expecting(SEMICOLON);
		}
	else if (expecting == Expecting(NUMBER))
		{
		if (strlen(tkn) == 1 && tkn[0] == '.')
			tempInt = numChar;
		else
			sscanf (tkn, "%d", &tempInt);
		if (tempInt <= 0 || tempInt > numChar)
			{
			MrBayesPrint ("%s   Character number %d is out of range (should be between %d and %d)\n", spacer, tempInt, 1, numChar);
			return (ERROR);
			}
		tempInt--;
		if (foundDash == YES)
			{
			if (fromI >= 0)
				toJ = tempInt;
			else
				{
				MrBayesPrint ("%s   Improperly formatted charset\n", spacer);
				return (ERROR);
				}
			foundDash = NO;
			}
		else if (foundSlash == YES)
			{
			tempInt++;
			if (tempInt <= 1)
				{
				MrBayesPrint ("%s   Improperly formatted charset\n", spacer);
				return (ERROR);
				}
			if (fromI >= 0 && toJ >= 0 && fromI < toJ)
				everyK = tempInt;
			else
				{
				MrBayesPrint ("%s   Improperly formatted charset\n", spacer);
				return (ERROR);
				}
			foundSlash = NO;
			}
		else
			{
			if (fromI >= 0 && toJ < 0)
				{
				if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
					return (ERROR);
				fromI = tempInt;
				}
			else if (fromI < 0 && toJ < 0)
				{
				fromI = tempInt;
				}
			else if (fromI >= 0 && toJ >= 0 && everyK < 0)
				{
				if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
					return (ERROR);
				fromI = tempInt;
				toJ = everyK = -1;
				}
			else if (fromI >= 0 && toJ >= 0 && everyK >= 0)
				{
				if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
					return (ERROR);
				fromI = tempInt;
				toJ = everyK = -1;
				}
			else
				{
				MrBayesPrint ("%s   Improperly formatted charset\n", spacer);
					{
					return (ERROR);
					}
				}
				
			}

		
		expecting  = Expecting(ALPHA);
		expecting |= Expecting(NUMBER);
		expecting |= Expecting(SEMICOLON);
		expecting |= Expecting(DASH);
		expecting |= Expecting(BACKSLASH);
		}
	else if (expecting == Expecting(DASH))
		{
		foundDash = YES;
		expecting = Expecting(NUMBER);
		}
	else if (expecting == Expecting(BACKSLASH))
		{
		foundSlash = YES;
		expecting = Expecting(NUMBER);
		}
	else
		return (ERROR);

	return (NO_ERROR);

}





int DoCharStat (void)

{

	int			i, j, numDivs;
	char		tempName[100];
	
	if (defMatrix == NO)
		{
		MrBayesPrint ("%s   A character matrix must be defined first\n", spacer);
		return (ERROR);
		}
			
	if (numDefinedPartitions == 1)
		MrBayesPrint ("%s   1 character partition defined:\n", spacer, numDefinedPartitions);
	else
		MrBayesPrint ("%s   %d character partitions defined:\n", spacer, numDefinedPartitions);
	for (i=0; i<numDefinedPartitions; i++)
		{
		if (GetNameFromString (partitionNames, tempName, i+1) == ERROR)
			{
			MrBayesPrint ("%s   Error getting partition names \n", spacer);
			return (ERROR);
			}
		numDivs = GetNumPartDivisions (i+1);
		if (numDivs == 1)
			MrBayesPrint ("%s      Partition %d (\"%s\") does not divide the characters\n", spacer, i+1, tempName);
		else
			MrBayesPrint ("%s      Partition %d (\"%s\") divides the characters into %d parts\n", spacer, i+1, tempName, numDivs);
		}
	if (GetNameFromString (partitionNames, tempName, partitionNum) == ERROR)
		{
		MrBayesPrint ("%s   Error getting current partition name \n", spacer);
		return (ERROR);
		}
	MrBayesPrint ("%s      Current partition is \"%s\"\n", spacer, tempName);
	MrBayesPrint ("\n");

	/* print out list of characters with information about each */
	MrBayesPrint ("%s   Showing character status:\n\n", spacer);
	MrBayesPrint ("%s                                                    Partition(s)\n", spacer);
	MrBayesPrint ("%s      #      Type      In/Out    Ambiguity Order  ", spacer);
	for (i=0; i<numDefinedPartitions; i++)
		MrBayesPrint (" %2d", i+1);
	MrBayesPrint ("\n");
	MrBayesPrint ("%s   -----------------------------------------------", spacer);
	for (i=0; i<numDefinedPartitions; i++)
		MrBayesPrint ("---");
	MrBayesPrint ("\n");
	for (i=0; i<numChar; i++)
		{
		MrBayesPrint ("%s   %4d -- ", spacer, i+1);
				
		if (charInfo[i].charType == DNA)
			MrBayesPrint ("   DNA");
		else if (charInfo[i].charType == RNA)
			MrBayesPrint ("   RNA");
		else if (charInfo[i].charType == PROTEIN)
			MrBayesPrint ("  Prot");
		else if (charInfo[i].charType == RESTRICTION)
			MrBayesPrint ("  Rest");
		else if (charInfo[i].charType == STANDARD)
			MrBayesPrint (" Stand");
		else if (charInfo[i].charType == CONTINUOUS)
			MrBayesPrint ("  Cont");
			
		if (charInfo[i].charType == DNA)
			MrBayesPrint ("   4");
		else if (charInfo[i].charType == RNA)
			MrBayesPrint ("   4");
		else if (charInfo[i].charType == PROTEIN)
			MrBayesPrint ("  20");
		else if (charInfo[i].charType == RESTRICTION)
			MrBayesPrint ("   2");
		else if (charInfo[i].charType == STANDARD)
			MrBayesPrint ("  %2d", charInfo[i].numStates);
		else if (charInfo[i].charType == CONTINUOUS)
			MrBayesPrint (" Inf");
			
		if (charInfo[i].isExcluded == NO)
			MrBayesPrint ("  Included");
		else
			MrBayesPrint ("  Excluded");
			
		if (charInfo[i].isMissAmbig == YES)
			MrBayesPrint ("  MissAmbig");
		else
			MrBayesPrint ("       None");
			
		if (charInfo[i].ctype == UNORD)
			MrBayesPrint (" Unord");
		else if (charInfo[i].ctype == ORD)
			MrBayesPrint ("   Ord");
		else if (charInfo[i].ctype == DOLLO)
			MrBayesPrint (" Dollo");
		else if (charInfo[i].ctype == IRREV)
			MrBayesPrint (" Irrev");

		MrBayesPrint ("  ");
			
		for (j=0; j<numDefinedPartitions; j++)
			MrBayesPrint (" %2d", charInfo[i].partitionId[j]);

		/* MrBayesPrint ("%4d   ", charSet[i]);*/
		
		if (charInfo[i].pairsId > 0)
			{
			/* find paired character */
			for (j=0; j<numChar; j++)
				{
				if (i != j && charInfo[j].pairsId == charInfo[i].pairsId)
					{
					MrBayesPrint (" (coupled with %d)", j+1);
					break;
					}
				}
			}
					
		MrBayesPrint ("\n");
		
		if (charInfo[i].bigBreakAfter == YES)
			{
			MrBayesPrint ("%s   ", spacer);
			MrBayesPrint ("     - - - - - - - - - - - - - - - - - - - -  \n");
			}
		
		/* we may want to pause */
		if (autoClose == NO)
			{
			if ((i+1) % 100 == 0)
				{
				MrBayesPrint ("%s   Hit return key to continue  ", spacer);
				fflush (stdin);
				fgets (tempName, 100, stdin);
				}
			}
		}

	return (NO_ERROR);

}





int DoCitations (void)

{

	MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
	MrBayesPrint ("   Citations                                                                     \n");
    MrBayesPrint ("                                                                                 \n");
	MrBayesPrint ("   If you publish results obtained using MrBayes you may want to cite the        \n");
    MrBayesPrint ("   program. The appropriate citation is:                                         \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Huelsenbeck, J. P. and F. Ronquist. 2001. MRBAYES: Bayesian                \n");
    MrBayesPrint ("         inference of phylogeny. Bioinformatics 17:754-755.                      \n");
    MrBayesPrint ("      Ronquist, F. and J. P. Huelsenbeck. 2003. MRBAYES 3: Bayesian phylogenetic \n");
    MrBayesPrint ("         inference under mixed models. Bioinformatics 19:1572-1574.              \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("   If you use the parallel abilities of the program, you may also want to cite   \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Altekar, G., S. Dwarkadas, J. P. Huelsenbeck, and F. Ronquist. 2004.       \n");
    MrBayesPrint ("         Parallel Metropolis-coupled Markov chain Monte Carlo for Bayesian       \n");
	MrBayesPrint ("         phylogenetic inference. Bioinformatics 20:407-415.                      \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("   You should also cite other papers for different ideas that are implemented    \n");
    MrBayesPrint ("   in the program. For example, the program performs Bayesian inference of       \n");
    MrBayesPrint ("   phylogeny, an idea that was first proposed in the following papers:           \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Larget, B., and D. Simon. 1999. Markov chain Monte Carlo                   \n");
    MrBayesPrint ("         algorithms for the Bayesian analysis of phylogenetic trees.             \n");
    MrBayesPrint ("         Mol. Biol. Evol. 16:750-759.                                            \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Li, S. 1996. Phylogenetic tree construction using Markov chain             \n");
    MrBayesPrint ("         Monte carlo. Ph. D. dissertation, Ohio State University, Columbus.      \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Mau, B. 1996. Bayesian phylogenetic inference via Markov chain             \n");
    MrBayesPrint ("         Monte carlo methods. Ph. D. dissertation, University of                 \n");
    MrBayesPrint ("         Wisconsin, Madison.                                                     \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Mau, B., and M. Newton. 1997. Phylogenetic inference for binary            \n");
    MrBayesPrint ("         data on dendrograms using Markov chain Monte Carlo. Journal of          \n");
    MrBayesPrint ("         Computational and Graphical Statistics 6:122-131.                       \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Mau, B., M. Newton, and B. Larget. 1999. Bayesian phylogenetic             \n");
    MrBayesPrint ("         inference via Markov chain Monte carlo methods. Biometrics. 55:1-12.    \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Newton, M., B. Mau, and B. Larget. 1999. Markov chain Monte Carlo          \n");
    MrBayesPrint ("         for the Bayesian analysis of evolutionary trees from aligned            \n");
    MrBayesPrint ("         molecular sequences. In Statistics in molecular biology (F. Seillier-   \n");
    MrBayesPrint ("         Moseiwitch, T. P. Speed, and M. Waterman, eds.). Monograph Series       \n");
    MrBayesPrint ("         of the Institute of Mathematical Statistics.                            \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Rannala, B., and Z. Yang. 1996. Probability distribution of                \n");
    MrBayesPrint ("         molecular evolutionary trees: a new method of phylogenetic              \n");
    MrBayesPrint ("         inference. J. Mol. Evol. 43:304-311.                                    \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Yang, Z., and B. Rannala. 1997. Bayesian phylogenetic inference            \n");
    MrBayesPrint ("         using DNA sequences: a Markov chain Monte carlo method. Molecular       \n");
    MrBayesPrint ("         Biology and Evolution. 14:717-724.                                      \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("   MrBayes uses Markov chain Monte Carlo (MCMC) to approximate the posterior     \n");
    MrBayesPrint ("   probability of trees. MCMC was developed in the following papers:             \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Metropolis, N., A. W. Rosenbluth, M. N. Rosenbluth, A. H. Teller,          \n");
    MrBayesPrint ("         and E. Teller. 1953. Equations of state calculations by fast            \n");
    MrBayesPrint ("         computing machines. J. Chem. Phys. 21:1087-1091.                        \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Hastings, W. K. 1970. Monte Carlo sampling methods using Markov            \n");
    MrBayesPrint ("         chains and their applications. Biometrika 57:97-109.                    \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("   In particular, MrBayes implements a variant of MCMC that was described by     \n");
    MrBayesPrint ("   Charles Geyer:                                                                \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Geyer, C. J. 1991. Markov chain Monte Carlo maximum likelihood.            \n");
    MrBayesPrint ("         Pages 156-163 in Computing Science and Statistics: Proceed-             \n");
    MrBayesPrint ("         ings of the 23rd Symposium on the Interface. (E. M. Keramidas,          \n");
    MrBayesPrint ("         ed.). Fairfax Station: Interface Foundation.                            \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("   MrBayes implements a large number of DNA substitution models. These models    \n");
    MrBayesPrint ("   are of three different structures. The \"4by4\" models are the usual          \n");
    MrBayesPrint ("   flavor of phylogenetic models. The \"Doublet\" model was first proposed       \n");
    MrBayesPrint ("   by                                                                            \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Schoniger, M., and A. von Haeseler. 1994. A stochastic model and the       \n");
    MrBayesPrint ("         evolution of autocorrelated DNA sequences. Molecular Phylogenetics      \n");
    MrBayesPrint ("         and Evolution 3:240-247.                                                \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("   The program also implements codon models. Two papers, published back-to-back  \n");
    MrBayesPrint ("   were the first to implement a codon model of DNA substitution in which the    \n");
    MrBayesPrint ("   substitution process is modelled on the codon, not on a site-by-site basis:   \n");
    MrBayesPrint ("                                                                                 \n");
	MrBayesPrint ("      Goldman, N., and Z. Yang. 1994. A codon-based model of nucleotide          \n");
	MrBayesPrint ("         substitution for protein coding DNA sequences. Molecular Biology        \n");
	MrBayesPrint ("         and Evolution. 11:725-736.                                              \n");
	MrBayesPrint ("                                                                                 \n");
	MrBayesPrint ("      Muse, S., and B. Gaut. 1994. A likelihood approach for comparing           \n");
	MrBayesPrint ("         synonymous and non-synonymous substitution rates, with application      \n");
	MrBayesPrint ("         to the chloroplast genome. Molecular Biology and Evolution.             \n");
	MrBayesPrint ("         11:715-724.                                                             \n");
	MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("   The program can be used to detect positively slected amino-acid sites using   \n");
    MrBayesPrint ("   a full hierarchical Bayes analysis. The method is based on the excellent paper\n");
    MrBayesPrint ("   by Nielsen and Yang:                                                          \n");
	MrBayesPrint ("                                                                                 \n");
	MrBayesPrint ("      Nielsen, R., and Z. Yang. 1998. Likelihood models for detecting            \n");
	MrBayesPrint ("         positively selected amino acid sites and applications to the HIV-1      \n");
	MrBayesPrint ("         envelope gene. Genetics. 148:929-936.                                   \n");
	MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("   The previous four papers describe three different stuctures for the nuc-      \n");
    MrBayesPrint ("   leotide models implemented in MrBayes--the four-by-four models, the           \n");
    MrBayesPrint ("   16-by-16 (doublet) models and the 64-by-64 (codon) models. The program        \n");
    MrBayesPrint ("   implements three different substitution models within each model structure.   \n");
    MrBayesPrint ("   These include the nst=1 models:                                               \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Jukes, T., and C. Cantor. 1969. Evolution of protein molecules.            \n");
    MrBayesPrint ("         Pages 21-132 in Mammalian Protein Metabolism. (H. Munro, ed.).          \n");
    MrBayesPrint ("         Academic Press, New York.                                               \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Felsenstein, J. 1981. Evolutionary trees from DNA sequences: A             \n");
    MrBayesPrint ("         maximum likelihood approach. Journal of Molecular Evolution             \n");
    MrBayesPrint ("         17:368-376.                                                             \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("   the nst=2 models:                                                             \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Kimura, M. 1980. A simple method for estimating evolutionary rates         \n");
    MrBayesPrint ("         of base substitutions through comparative studies of nucleotide         \n");
    MrBayesPrint ("         sequences. Journal of Molecular Evolution. 16:111-120.                  \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Hasegawa, M., T. Yano, and H. Kishino. 1984. A new molecular clock         \n");
    MrBayesPrint ("         of mitochondrial DNA and the evolution of Hominoids. Proc.              \n");
    MrBayesPrint ("         Japan Acad. Ser. B 60:95-98.                                            \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Hasegawa, M., H. Kishino, and T. Yano. 1985. Dating the human-ape          \n");
    MrBayesPrint ("         split by a molecular clock of mitochondrial DNA. Journal of             \n");
    MrBayesPrint ("         Molecular Evolution 22:160-174.                                         \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("   and the the nst=6 models:                                                     \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Tavare, S. 1986. Some probabilistic and statisical problems on the         \n");
    MrBayesPrint ("         analysis of DNA sequences. Lect. Math. Life Sci. 17:57-86.              \n");
    MrBayesPrint ("         17:368-376.                                                             \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("   MrBayes implements a large number of amino-acid models. These include:        \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Poisson --                                                                 \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Bishop, M.J., and A.E. Friday. 1987. Tetropad relationships: the           \n");
    MrBayesPrint ("         molecular evidence. Pp. 123�139 in Molecules and morphology in          \n");
    MrBayesPrint ("         evolution: conflict or compromise? (C. Patterson, ed.). Cambridge       \n");
    MrBayesPrint ("         University Press, Cambridge, England.                                   \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Jones --                                                                   \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Jones, D.T., W. R. Taylor, and J. M. Thornton. 1992. The rapid generation  \n");
    MrBayesPrint ("         of mutation data matrices from protein sequences. Comput. Appl.         \n");
    MrBayesPrint ("         Biosci. 8:275�282.                                                      \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Dayhoff --                                                                 \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Dayhoff, M.O., R.M. Schwartz, and B.C. Orcutt. 1978. A model of evol-      \n");
    MrBayesPrint ("         utionary change in proteins. Pp. 345�352 in Atlas of protein sequence   \n");
    MrBayesPrint ("         and structure. Vol. 5, Suppl. 3. National Biomedical Research           \n");
    MrBayesPrint ("          Foundation, Washington, D.C.                                           \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Mtrev --                                                                   \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Adachi, J. and M. Hasegawa. 1996. MOLPHY version 2.3: programs for         \n");
    MrBayesPrint ("         molecular phylogenetics based on maximum likelihood.  Computer Science  \n");
    MrBayesPrint ("         Monographs of Institute of Statistical Mathematics 28:1-150.            \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Mtmam --                                                                   \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Cao, Y., A. Janke, P.J. Waddell, M. Westerman, O. Takenaka, S. Murata,     \n");
    MrBayesPrint ("         N. Okada, S. Paabo, and M. Hasegawa. 1998. Conflict amongst individual  \n");
    MrBayesPrint ("         mitochondrial proteins in resolving the phylogeny of eutherian orders.  \n");
    MrBayesPrint ("         Journal of Molecular Evolution                                          \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Yang, Z., R. Nielsen, and M. Hasegawa. 1998.  Models of amino acid         \n");
    MrBayesPrint ("         substitution and applications to mitochondrial protein evolution        \n");
    MrBayesPrint ("         Molecular Biology and Evolution 15:1600�1611.                           \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      WAG --                                                                     \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Whelan, S. and Goldman, N. 2001. A general empirical model of protein      \n");
    MrBayesPrint ("         evolution derived from multiple protein families using a maximum-       \n");
    MrBayesPrint ("         likelihood approach. Molecular Biology and Evolution 18:691-699.        \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Rtrev --                                                                   \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Dimmic M.W., J.S. Rest, D.P. Mindell, and D. Goldstein. 2002. RArtREV:     \n");
    MrBayesPrint ("         An amino acid substitution matrix for inference of retrovirus and       \n");
    MrBayesPrint ("         reverse transcriptase phylogeny. Journal of Molecular Evolution         \n");
    MrBayesPrint ("         55: 65-73.                                                              \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Cprev --                                                                   \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Adachi, J., P. Waddell, W. Martin, and M. Hasegawa. 2000. Plastid          \n");
    MrBayesPrint ("         genome phylogeny and a model of amino acid substitution for proteins    \n");
    MrBayesPrint ("         encoded by chloroplast DNA. Journal of Molecular Evolution              \n");
    MrBayesPrint ("         50:348-358.                                                             \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Blosum --                                                                  \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Henikoff, S., and J. G. Henikoff. 1992. Amino acid substitution            \n");
    MrBayesPrint ("         matrices from protein blocks. Proc. Natl. Acad. Sci., U.S.A.            \n");
    MrBayesPrint ("         89:10915-10919. The matrix implemented in MrBayes is Blosum62.          \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Vt --                                                                      \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Muller, T., and M. Vingron. 2000. Modeling amino acid replacement.         \n");
    MrBayesPrint ("         Journal of Computational Biology 7:761-776.                             \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("   MrBayes implements a simple Jukes-Cantor-like model for restriction sites     \n");
    MrBayesPrint ("   and other binary data. A problem with some of these data is that there is a   \n");
    MrBayesPrint ("   coding bias, such that certain characters are missing from any observable     \n");
    MrBayesPrint ("   data matrix. It is impossible, for instance, to observe restriction sites that\n");
    MrBayesPrint ("   are absent in all the studied taxa. However, MrBayes corrects for this coding \n");
    MrBayesPrint ("   bias according to an idea described in                                        \n");
    MrBayesPrint ("                                                                                 \n");
	MrBayesPrint ("      Felsenstein, J. 1992. Phylogenies from restriction sites: A maximum-       \n");
	MrBayesPrint ("         likelihood approach. Evolution 46:159-173.                              \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("   The model used by MrBayes for 'standard' or morphological data is based on    \n");
    MrBayesPrint ("   the ideas originally presented by                                             \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Lewis, P. O. 2001. A likelihood approach to estimating phylogeny from      \n");
	MrBayesPrint ("         discrete morphological character data. Systematic Biology 50:913-925.   \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("                                                                                 \n");
	MrBayesPrint ("   For both DNA sequence and amino-acid data, the program allows rates to        \n");
    MrBayesPrint ("   change under a covarion-like model, first described by Tuffley and Steel      \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Tuffley, C., and M. Steel. 1998. Modeling the covarion hypothesis          \n");
    MrBayesPrint ("         of nucleotide substitution. Mathematical Biosciences 147:63-91.         \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("   and implemented by Huelsenbeck (2002)                                         \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Huelsenbeck, J. P. 2002. Testing a covariotide model of DNA sub-           \n");
    MrBayesPrint ("         stitution. Molecular Biology and Evolution 19(5):698-707.               \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("   Galtier (2001) implements a different variant of the covarion model in        \n");
    MrBayesPrint ("   a paper that is worth reading:                                                \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Galtier, N. 2001. Maximum-likelihood phylogenetic analysis under a         \n");
    MrBayesPrint ("         covarion-like model. Mol. Biol. Evol. 18:866-873.                       \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("   A number of models are available that allow rates to vary                     \n");
    MrBayesPrint ("   across the characters. The program implements the proportion                  \n");
    MrBayesPrint ("   of invariable sites model and two variants of gamma distributed               \n");
    MrBayesPrint ("   rate variation. Yang\'s (1993) paper is a good one to cite for                \n");
    MrBayesPrint ("   implementing a gamma-distributed rates model. In the 1994 paper he            \n");
    MrBayesPrint ("   provides a way to approximate the continuous gamma distribution:              \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Yang, Z. 1993. Maximum likelihood estimation of phylogeny from DNA         \n");
    MrBayesPrint ("         sequences when substitution rates differ over sites. Molecular          \n");
    MrBayesPrint ("         Biology and Evolution 10:1396-1401.                                     \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Yang, Z. 1994. Maximum likelihood phylogenetic estimation from DNA         \n");
    MrBayesPrint ("         sequences with variable rates over sites: Approximate methods.          \n");
    MrBayesPrint ("         Journal of Molecular Evolution 39:306-314.                              \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("   The program also implements Yang\'s autocorrelated gamma model. In            \n");
    MrBayesPrint ("   this model, the rate at one site depends to some extent on the rate at        \n");
    MrBayesPrint ("   an adjacent site. The appropriate citation for this model is:                 \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Yang, Z. 1995. A space-time process model for the evolution of             \n");
    MrBayesPrint ("         DNA sequences. Genetics 139:993-1005.                                   \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("   Ancestral state reconstruction. These two papers show how ancestral           \n");
    MrBayesPrint ("   states on a tree can be reconstructed. The Yang et al. paper                  \n");
    MrBayesPrint ("   implements an empirical Bayes approach. The Huelsenbeck and Bollback          \n");
    MrBayesPrint ("   paper implements a hierarchical Bayes approach. The method implemented        \n");
    MrBayesPrint ("   in MrBayes is hierarchical Bayes as it integrates over uncertainty in         \n");
    MrBayesPrint ("   model parameters. --                                                          \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Yang, Z., S. Kumar, and M. Nei. 1995. A new method of inference of         \n");
    MrBayesPrint ("         ancestral nucleotide and amino acid sequences. Genetics 141:1641        \n");
    MrBayesPrint ("         1650.                                                                   \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Huelsenbeck, J. P., and J. P. Bollback. 2001. Empirical and hier-          \n");
    MrBayesPrint ("         archical Bayesian estimation of ancestral states. Systematic            \n");
    MrBayesPrint ("         Biology 50:351-366.                                                     \n");
    MrBayesPrint ("                                                                                 \n");
	MrBayesPrint ("   You may also want to consult a more recent review of Bayesian reconstruction  \n");
	MrBayesPrint ("   of ancestral states and character evolution:                                  \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Ronquist, F. 2004. Bayesian inference of character evolution. Trends in    \n");
	MrBayesPrint ("         Ecology and Evolution 19: 475-481.                                      \n");
    MrBayesPrint ("                                                                                 \n");
	MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("   The program implements an incredibly parameter rich model, first described    \n");
    MrBayesPrint ("   by Tuffley and Steel (1997), that orders trees in the same way as the         \n");
    MrBayesPrint ("   so-called parsimony method of phylogenetic inference. The appropriate         \n");
    MrBayesPrint ("   citation is:                                                                  \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Tuffley, C., and M. Steel. 1997. Links between maximum likelihood          \n");
    MrBayesPrint ("         and maximum parsimony under a simple model of site substitution.        \n");
    MrBayesPrint ("         Bull. Math. Bio. 59:581-607.                                            \n");
	MrBayesPrint ("   ---------------------------------------------------------------------------   \n");

	return (NO_ERROR);
	
}





int DoConstraints (void)

{

	int			i, howMany;

	/* first add name to constraintNames */
	if (AddToString (tempSetName, constraintNames, &howMany) == ERROR)
		{
		MrBayesPrint ("%s   Problem adding constraint %s to list\n", spacer, tempSetName);
		return (ERROR);
		}
	if (howMany != numDefinedConstraints + 1)
		{
		MrBayesPrint ("%s   Problem adding constraint %s to list\n", spacer, tempSetName);
		if (RemoveLastFromString (constraintNames) == ERROR)
			return (ERROR);
		return (ERROR);
		}

	/* add set to tempSet */
	if (fromI >= 0 && toJ < 0)
		{
		if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
			{
			if (RemoveLastFromString (constraintNames) == ERROR)
				return (ERROR);
			return (ERROR);
			}
		}
	else if (fromI >= 0 && toJ >= 0 && everyK < 0)
		{
		if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
			{
			if (RemoveLastFromString (constraintNames) == ERROR)
				return (ERROR);
			return (ERROR);
			}
		}
	else if (fromI >= 0 && toJ >= 0 && everyK >= 0)
		{
		if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
			{
			if (RemoveLastFromString (constraintNames) == ERROR)
				return (ERROR);
			return (ERROR);
			}
		}
			
	/* check that this is not a stupid constraint */
	howMany = 0;
	for (i=0; i<numTaxa; i++)
		if (tempSet[i] != 0)
			howMany++;
	if (howMany >= numTaxa)
		{
		MrBayesPrint ("%s   This constraint includes all taxa and will not be defined\n", spacer);
		if (RemoveLastFromString (constraintNames) == ERROR)
			return (ERROR);
		return (ERROR);
		}
	else if (howMany == 0)
		{
		MrBayesPrint ("%s   This constraint does not include any taxa and will not be defined\n", spacer);
		if (RemoveLastFromString (constraintNames) == ERROR)
			return (ERROR);
		return (ERROR);
		}
	else if (howMany == 1 || howMany == numTaxa - 1)
		MrBayesPrint ("%s   This is a trivial constraint for unrooted trees\n", spacer);
		
	/* merge tempSet with constraints */
	for (i=0; i<numTaxa; i++)
		taxaInfo[i].constraints[numDefinedConstraints] = tempSet[i];
	
	/* increment number of char sets */
	numDefinedConstraints++;
	
	/* show taxset (for debugging) */
#	if 0
	for (i=0; i<numTaxa; i++)
		MrBayesPrint ("%4d  %4d\n", i+1, taxaInfo[i].constraints[numDefinedConstraints-1]);
#	endif

	return (NO_ERROR);
	
}





int DoConstraintsParm (char *parmName, char *tkn)

{

	int		i, howMany, tempInt;
	MrBFlt	tempD;
	
	if (defMatrix == NO)
		{
		MrBayesPrint ("%s   A matrix must be specified before constraints can be defined\n", spacer);
		return (ERROR);
		}

	if (expecting == Expecting(PARAMETER))
		{
		if (!strcmp(parmName, "Xxxxxxxxxx"))
			{
			/* check size of constraint name */
			if (strlen(tkn) > 99)
				{
				MrBayesPrint ("%s   Constraint name is too long\n", spacer);
				return (ERROR);
				}
				
			/* check to see if the name has already been used as a taxset */
			if (numDefinedConstraints > 0)
				{
				if (CheckString (tkn, constraintNames, &howMany) == ERROR)
					{
					/* if the constraint name has not been used, then we should have an ERROR returned */
					/* we _want_ to be here */

					}
				else
					{
					MrBayesPrint ("%s   Constraint name has been used previously\n", spacer);
					return (ERROR);
					}
				}
			else if (numDefinedConstraints > 30)
				{
				MrBayesPrint ("%s   You cannot define more than 30 constraints\n", spacer);
				return (ERROR);
				}
				
			/* add the name to the temporary constraint names string */
			strcpy (tempSetName, tkn);
			
			/* clear tempSet */
			for (i=0; i<numTaxa; i++)
				tempSet[i] = 0;
			
			fromI = toJ = everyK = -1;
			foundDash = foundSlash = NO;
			MrBayesPrint ("%s   Defining constraint called %s\n", spacer, tkn);
			foundExp = NO;
			foundFirst = YES;
			isNegative = NO;
			expecting = Expecting(ALPHA);
			expecting |= Expecting(NUMBER);
			expecting |= Expecting(DASH);
			}
		else
			return (ERROR);
		}

	else if (expecting == Expecting(EQUALSIGN))
		{
		expecting  = Expecting(ALPHA);
		expecting |= Expecting(NUMBER);
		}
	else if (expecting == Expecting(LEFTPAR))
		{
		isNegative = NO;
		expecting = Expecting(NUMBER);
		expecting |= Expecting(DASH);
		}
	else if (expecting == Expecting(RIGHTPAR))
		{
		isNegative = NO;
		foundExp = NO;
		expecting = Expecting(EQUALSIGN);
		}
	else if (expecting == Expecting(DASH))
		{
		if (foundExp == YES)
			isNegative = YES;
		else
			foundDash = YES;
		expecting = Expecting(NUMBER);
		}
	else if (expecting == Expecting(ALPHA))
		{
		if (foundFirst == YES)
			{
			/* We are filling in the probability for the constraint. Specifically, we expect exp(number). */
			if (IsSame ("Exp", tkn) == SAME || IsSame ("Exp", tkn) == CONSISTENT_WITH)
				{
				foundExp = YES;
				foundDash = NO;
				isNegative = NO;
				}
			else
				{
				MrBayesPrint ("%s   Do not understand %s\n", spacer, tkn);
				return (ERROR);
				}
			expecting  = Expecting(LEFTPAR);
			}
		else
			{
			/* We are defining a constraint in terms of a taxon set (called tkn, here) or we are referring to
			   the taxon name. We should be able to find tkn in the list of taxon set names or in the list
			   of taxon names. If we cannot, then we have a problem and return an error. */
			if (CheckString (tkn, taxaNames, &howMany) == ERROR)
				{
				if (numTaxaSets < 1)
					{
					MrBayesPrint ("%s   Could not find a taxset called %s\n", spacer, tkn);
					return (ERROR);
					}
				if (CheckString (tkn, taxaSetNames, &howMany) == ERROR)
					{
					MrBayesPrint ("%s   Could not find a taxset called %s\n", spacer, tkn);
					return (ERROR);
					}
				/* add taxa from taxset tkn to new tempSet */
				for (i=0; i<numTaxa; i++)
					{
					if (taxaInfo[i].taxaSet[howMany-1] == 1)
						tempSet[i] = 1;
					}
				}
			else
				{
				tempSet[howMany-1] = 1;
				}
			fromI = toJ = everyK = -1;

			expecting  = Expecting(ALPHA);
			expecting |= Expecting(NUMBER);
			expecting |= Expecting(SEMICOLON);
			}
		}
	else if (expecting == Expecting(NUMBER))
		{
		if (foundFirst == YES)
			{
			/* We are filling in the probability for the constraint. Specifically, we expect number. */
			sscanf (tkn, "%lf", &tempD);		
			if (foundExp == NO && tempD < 0.0)
				{
				MrBayesPrint ("%s   The probability of a clade cannot be less than zero\n", spacer, tkn);
				return (ERROR);
				}
			if (isNegative == YES || foundDash == YES)
				tempD *= -1.0;
			if (foundExp == YES)
				{
				relConstraintProbs[numDefinedConstraints] = tempD;
				expecting  = Expecting(RIGHTPAR);
				}
			else
				{
				if (tempD <= 0.0)
					relConstraintProbs[numDefinedConstraints] = NEG_INFINITY;
				else
					relConstraintProbs[numDefinedConstraints] = log(tempD);
				expecting  = Expecting(EQUALSIGN);
				}
			foundFirst = NO;
			foundDash = NO;
			}
		else
			{		
			if (strlen(tkn) == 1 && !strcmp(tkn, "."))
				{
				tempInt = numTaxa;
				}
			else
				{
				sscanf (tkn, "%d", &tempInt);
				if (tempInt <= 0 || tempInt > numTaxa)
					{
					MrBayesPrint ("%s   Taxon number %d is out of range (should be between %d and %d)\n", spacer, tempInt, 1, numTaxa);
					return (ERROR);
					}
				}
			tempInt--;
			if (foundDash == YES)
				{
				if (fromI >= 0)
					toJ = tempInt;
				else
					{
					MrBayesPrint ("%s   Improperly formatted constraint\n", spacer);
					return (ERROR);
					}
				foundDash = NO;
				}
			else if (foundSlash == YES)
				{
				tempInt++;
				if (tempInt <= 1)
					{
					MrBayesPrint ("%s   Improperly formatted constraint\n", spacer);
					return (ERROR);
					}
				if (fromI >= 0 && toJ >= 0 && fromI < toJ)
					everyK = tempInt;
				else
					{
					MrBayesPrint ("%s   Improperly formatted constraint\n", spacer);
					return (ERROR);
					}
				foundSlash = NO;
				}
			else
				{
				if (fromI >= 0 && toJ < 0)
					{
					if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
						return (ERROR);
					fromI = tempInt;
					}
				else if (fromI < 0 && toJ < 0)
					{
					fromI = tempInt;
					}
				else if (fromI >= 0 && toJ >= 0 && everyK < 0)
					{
					if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
						return (ERROR);
					fromI = tempInt;
					toJ = everyK = -1;
					}
				else if (fromI >= 0 && toJ >= 0 && everyK >= 0)
					{
					if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
						return (ERROR);
					fromI = tempInt;
					toJ = everyK = -1;
					}
				else
					{
					MrBayesPrint ("%s   Improperly formatted constraint\n", spacer);
						{
						return (ERROR);
						}
					}
				}

			expecting  = Expecting(ALPHA);
			expecting |= Expecting(NUMBER);
			expecting |= Expecting(SEMICOLON);
			expecting |= Expecting(DASH);
			expecting |= Expecting(BACKSLASH);
			}
		}
	else if (expecting == Expecting(BACKSLASH))
		{
		foundSlash = YES;
		expecting = Expecting(NUMBER);
		}
	else
		return (ERROR);

	return (NO_ERROR);

}





int DoCtype (void)

{

	int			i, foundIllegal, marks[5], numAppliedTo;

	/* add set to tempSet */
	if (fromI >= 0 && toJ < 0)
		{
		if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
			return (ERROR);
		}
	else if (fromI >= 0 && toJ >= 0 && everyK < 0)
		{
		if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
			return (ERROR);
		}
	else if (fromI >= 0 && toJ >= 0 && everyK >= 0)
		{
		if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
			return (ERROR);
		}
		
	/* merge tempSet with ctype */
	numAppliedTo = 0;
	for (i=0; i<5; i++)
		marks[i] = NO;
	for (i=0; i<numChar; i++)
		{
		if (tempSet[i] != 0)
			{
			foundIllegal = NO;
			if (charOrdering != UNORD)
				{
				if (charInfo[i].charType == DNA)
					{
					foundIllegal = YES;
					if (marks[0] == NO)
						MrBayesPrint ("%s   Ctype not applied to DNA states which must be unordered\n", spacer);
					marks[0] = YES;
					}
				else if (charInfo[i].charType == RNA)
					{
					foundIllegal = YES;
					if (marks[1] == NO)
						MrBayesPrint ("%s   Ctype not applied to RNA states which must be unordered\n", spacer);
					marks[1] = YES;
					}
				else if (charInfo[i].charType == PROTEIN)
					{
					foundIllegal = YES;
					if (marks[2] == NO)
						MrBayesPrint ("%s   Ctype not applied to amino acid states which must be unordered\n", spacer);
					marks[2] = YES;
					}
				else if (charInfo[i].charType == RESTRICTION)
					{
					foundIllegal = YES;
					if (marks[3] == NO)
						MrBayesPrint ("%s   Ctype not applied to restriction site states which must be unordered\n", spacer);
					marks[3] = YES;
					}
				else if (charInfo[i].charType == CONTINUOUS)
					{
					foundIllegal = YES;
					if (marks[4] == NO)
						MrBayesPrint ("%s   Ctype not applied to continuous characters\n", spacer);
					marks[4] = YES;
					}
				}
			if (foundIllegal == NO)
				{
				charInfo[i].ctype = charOrdering;
				numAppliedTo++;
				}
			}
		}
	if (numAppliedTo > 0)
		{
		MrBayesPrint ("%s   Ctype was applied to %d standard characters\n", spacer, numAppliedTo);
		}
	else
		{
		MrBayesPrint ("%s   No standard characters found to apply ctype to\n", spacer);
		}
	
#	if 0
	for (i=0; i<numChar; i++)
		MrBayesPrint ("%4d -- %d\n", i, ctype[i]);
#	endif

	return (NO_ERROR);
	
}





int DoCtypeParm (char *parmName, char *tkn)

{

	int		i, howMany, tempInt;
	
	if (defMatrix == NO)
		{
		MrBayesPrint ("%s   A matrix must be specified before typesets can be defined\n", spacer);
		return (ERROR);
		}

	if (expecting == Expecting(PARAMETER))
		{
		if (!strcmp(parmName, "Xxxxxxxxxx"))
			{
			if (IsSame ("Ordered", tkn) == SAME || IsSame ("Ordered", tkn) == CONSISTENT_WITH)
				charOrdering = ORD;
			else if (IsSame ("Unordered", tkn) == SAME || IsSame ("Unordered", tkn) == CONSISTENT_WITH)
				charOrdering = UNORD;
			else if (IsSame ("Dollo", tkn) == SAME || IsSame ("Dollo", tkn) == CONSISTENT_WITH)
				charOrdering = DOLLO;
			else if (IsSame ("Irreversible", tkn) == SAME || IsSame ("Irreversible", tkn) == CONSISTENT_WITH)
				charOrdering = IRREV;
			else
				{
				MrBayesPrint ("%s   Do not understand delimiter \"%s\"\n", spacer, tkn);
				return (ERROR);
				}
			
			/* clear tempSet */
			for (i=0; i<numChar; i++)
				tempSet[i] = 0;
			
			fromI = toJ = everyK = -1;
			foundDash = foundSlash = NO;
			MrBayesPrint ("%s   Setting characters to %s\n", spacer, tkn);
			expecting = Expecting(COLON);
			}
		else
			return (ERROR);
		}
	else if (expecting == Expecting(COLON))
		{
		expecting  = Expecting(ALPHA) | Expecting(NUMBER);
		}
	else if (expecting == Expecting(ALPHA))
		{
		/* first, check that we are not trying to put in another character ordering */
		if (IsSame ("Ordered", tkn) == SAME || IsSame ("Ordered", tkn) == CONSISTENT_WITH)
			{
			MrBayesPrint ("%s   You cannot specify more than one ordering with a single use of ctype\n", spacer, tkn);
			return (ERROR);
			}
		else if (IsSame ("Unordered", tkn) == SAME || IsSame ("Unordered", tkn) == CONSISTENT_WITH)
			{
			MrBayesPrint ("%s   You cannot specify more than one ordering with a single use of ctype\n", spacer, tkn);
			return (ERROR);
			}
		else if (IsSame ("Dollo", tkn) == SAME || IsSame ("Dollo", tkn) == CONSISTENT_WITH)
			{
			MrBayesPrint ("%s   You cannot specify more than one ordering with a single use of ctype\n", spacer, tkn);
			return (ERROR);
			}
		else if (IsSame ("Irreversible", tkn) == SAME || IsSame ("Irreversible", tkn) == CONSISTENT_WITH)
			{
			MrBayesPrint ("%s   You cannot specify more than one ordering with a single use of ctype\n", spacer, tkn);
			return (ERROR);
			}
		
		/* We are defining a type set in terms of another (called tkn, here). We should be able
		   to find tkn in the list of character set names. If we cannot, then we have a problem and
		   return an error. */
		if (IsSame ("All", tkn) == SAME)
			{
			for (i=0; i<numChar; i++)
				tempSet[i] = 1;
			fromI = toJ = everyK = -1;
			expecting = Expecting(SEMICOLON);
			}
		else
			{
			if (numCharSets < 1)
				{
				MrBayesPrint ("%s   Could not find a character set called %s\n", spacer, tkn);
				return (ERROR);
				}
			if (CheckString (tkn, charSetNames, &howMany) == ERROR)
				{
				MrBayesPrint ("%s   Could not find a character set called %s\n", spacer, tkn);
				return (ERROR);
				}
				
			/* add characters from charset tkn to new tempset */
			for (i=0; i<numChar; i++)
				{
				if (charInfo[i].charSet[howMany-1] == 1)
					tempSet[i] = 1;
				}
			fromI = toJ = everyK = -1;
			expecting  = Expecting(ALPHA);
			expecting |= Expecting(NUMBER);
			expecting |= Expecting(SEMICOLON);
			}
		}
	else if (expecting == Expecting(NUMBER))
		{
		if (strlen(tkn) == 1 && tkn[0] == '.')
			tempInt = numChar;
		else
			sscanf (tkn, "%d", &tempInt);
		if (tempInt <= 0 || tempInt > numChar)
			{
			MrBayesPrint ("%s   Character number %d is out of range (should be between %d and %d)\n", spacer, tempInt, 1, numChar);
			return (ERROR);
			}
		tempInt--;
		if (foundDash == YES)
			{
			if (fromI >= 0)
				toJ = tempInt;
			else
				{
				MrBayesPrint ("%s   Improperly formatted ctype\n", spacer);
				return (ERROR);
				}
			foundDash = NO;
			}
		else if (foundSlash == YES)
			{
			tempInt++;
			if (tempInt <= 1)
				{
				MrBayesPrint ("%s   Improperly formatted ctype\n", spacer);
				return (ERROR);
				}
			if (fromI >= 0 && toJ >= 0 && fromI < toJ)
				everyK = tempInt;
			else
				{
				MrBayesPrint ("%s   Improperly formatted ctype\n", spacer);
				return (ERROR);
				}
			foundSlash = NO;
			}
		else
			{
			if (fromI >= 0 && toJ < 0)
				{
				if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
					return (ERROR);
				fromI = tempInt;
				}
			else if (fromI < 0 && toJ < 0)
				{
				fromI = tempInt;
				}
			else if (fromI >= 0 && toJ >= 0 && everyK < 0)
				{
				if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
					return (ERROR);
				fromI = tempInt;
				toJ = everyK = -1;
				}
			else if (fromI >= 0 && toJ >= 0 && everyK >= 0)
				{
				if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
					return (ERROR);
				fromI = tempInt;
				toJ = everyK = -1;
				}
			else
				{
				MrBayesPrint ("%s   Improperly formatted ctype\n", spacer);
					{
					return (ERROR);
					}
				}
				
			}

		
		expecting  = Expecting(ALPHA);
		expecting |= Expecting(NUMBER);
		expecting |= Expecting(SEMICOLON);
		expecting |= Expecting(DASH);
		expecting |= Expecting(BACKSLASH);
		}
	else if (expecting == Expecting(DASH))
		{
		foundDash = YES;
		expecting = Expecting(NUMBER);
		}
	else if (expecting == Expecting(BACKSLASH))
		{
		foundSlash = YES;
		expecting = Expecting(NUMBER);
		}
	else
		return (ERROR);

	return (NO_ERROR);

}





int DoDelete (void)

{

	int			i, alreadyDone;

	MrBayesPrint ("%s   Excluding taxa\n", spacer);

	/* add set to tempSet */
	if (fromI >= 0 && toJ < 0)
		{
		if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
			return (ERROR);
		}
	else if (fromI >= 0 && toJ >= 0 && everyK < 0)
		{
		if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
			return (ERROR);
		}
	else if (fromI >= 0 && toJ >= 0 && everyK >= 0)
		{
		if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
			return (ERROR);
		}
		
	/* merge tempSet with taxset */
	alreadyDone = NO;
	for (i=0; i<numTaxa; i++)
		{
		if (tempSet[i] == 1)
			{
			if (taxaInfo[i].isDeleted == YES && alreadyDone == NO)
				{
				MrBayesPrint ("%s   Some taxa already excluded\n", spacer);
				alreadyDone = YES;
				}
			taxaInfo[i].isDeleted = YES;
			}
		}

	/* show tempSet (for debugging) */
#	if 0
	for (i=0; i<numTaxa; i++)
		MrBayesPrint ("%4d  %4d\n", i+1, tempSet[i]);
#	endif

	return (NO_ERROR);

}





int DoDeleteParm (char *parmName, char *tkn)

{

	int		i, howMany, tempInt;
		
	if (defMatrix == NO)
		{
		MrBayesPrint ("%s   A matrix must be specified before you can delete taxa\n", spacer);
		return (ERROR);
		}
		
	if (foundFirst == NO)
		{
		/* this is the first time in */
		fromI = toJ = everyK = -1;
		foundDash = NO;
		for (i=0; i<numTaxa; i++) /* clear tempSet */
			tempSet[i] = 0;
		foundFirst = YES;
		}

	if (expecting == Expecting(ALPHA))
		{
		if (IsSame ("All", tkn) == SAME || IsSame ("All", tkn) == CONSISTENT_WITH)
			{
			for (i=0; i<numTaxa; i++)
				tempSet[i] = 1;
			}
		else
			{
			if (CheckString (tkn, taxaNames, &howMany) == ERROR)
				{
				/* we are using a pre-defined taxa set */
				if (numTaxaSets < 1)
					{
					MrBayesPrint ("%s   Could not find a taxset called %s\n", spacer, tkn);
					return (ERROR);
					}
				if (CheckString (tkn, taxaSetNames, &howMany) == ERROR)
					{
					MrBayesPrint ("%s   Could not find a taxset called %s\n", spacer, tkn);
					return (ERROR);
					}
				/* add taxa from taxset tkn to new tempSet */
				for (i=0; i<numTaxa; i++)
					{
					if (taxaInfo[i].taxaSet[howMany-1] == 1)
						tempSet[i] = 1;
					}
				}
			else
				{
				/* we found the taxon name */
				if (fromI >= 0 && toJ < 0)
					{
					if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
						return (ERROR);
					}
				else if (fromI >= 0 && toJ >= 0)
					{
					if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
						return (ERROR);
					}
					
				tempSet[howMany-1] = 1;
				}
			}
		foundDash = NO;
		fromI = toJ = everyK = -1;
		
		expecting  = Expecting(ALPHA);
		expecting |= Expecting(NUMBER);
		expecting |= Expecting(SEMICOLON);
		}
	else if (expecting == Expecting(NUMBER))
		{
		if (strlen(tkn) == 1 && !strcmp(tkn, "."))
			tempInt = numTaxa;
		else
			{
			sscanf (tkn, "%d", &tempInt);
			if (tempInt <= 0 || tempInt > numTaxa)
				{
				MrBayesPrint ("%s   Taxon number %d is out of range (should be between %d and %d)\n", spacer, tempInt, 1, numTaxa);
				return (ERROR);
				}
			}
		tempInt--;
		if (foundDash == YES)
			{
			if (fromI >= 0)
				toJ = tempInt;
			else
				{
				MrBayesPrint ("%s   Improperly formatted delete set\n", spacer);
				return (ERROR);
				}
			foundDash = NO;
			}
		else
			{
			if (fromI >= 0 && toJ < 0)
				{
				if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
					return (ERROR);
				fromI = tempInt;
				}
			else if (fromI < 0 && toJ < 0)
				{
				fromI = tempInt;
				}
			else if (fromI >= 0 && toJ >= 0 && everyK < 0)
				{
				if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
					return (ERROR);
				fromI = tempInt;
				toJ = everyK = -1;
				}
			else if (fromI >= 0 && toJ >= 0 && everyK >= 0)
				{
				if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
					return (ERROR);
				fromI = tempInt;
				toJ = everyK = -1;
				}
			else
				{
				MrBayesPrint ("%s   Improperly formatted delete set\n", spacer);
					{
					return (ERROR);
					}
				}
			}
		expecting  = Expecting(ALPHA);
		expecting |= Expecting(NUMBER);
		expecting |= Expecting(SEMICOLON);
		expecting |= Expecting(DASH);
		}
	else if (expecting == Expecting(DASH))
		{
		foundDash = YES;
		expecting = Expecting(NUMBER);
		}
	else
		return (ERROR);

	return (NO_ERROR);
	MrBayesPrint ("%s", parmName);  /*just because I am tired of seeing the unused parameter error msg */
	
}





int DoDeroot (void)

{

	if (isUserTreeDefined == NO)
		{
		MrBayesPrint ("   No tree in memory\n", spacer);
		return (ERROR);
		}
	if (userTree->isRooted == YES)
		{
		if (DerootTree (userTree, outGroupNum) == ERROR)
			return (ERROR);
		MrBayesPrint ("%s   Derooting user tree\n", spacer);
		}
	else
		MrBayesPrint ("%s   Tree is already unrooted\n", spacer);
		
	return (NO_ERROR);
	
}





int DoDimensions (void)

{

	if (defTaxa == NO || defChars == NO)
		{
		MrBayesPrint ("%s   Expecting both Ntax and Nchar to be defined\n", spacer);
		return (ERROR);
		}
	if (inDataBlock == NO)
		{
		MrBayesPrint ("%s   Dimensions can only be defined in a data block\n", spacer);
		return (ERROR);
		}
	/* allocate matrix */
	if (AllocMatrix() == ERROR)
		return (ERROR);
	MrBayesPrint ("%s   Matrix has %d taxa and %d characters\n", spacer, numTaxa, numChar);
	return (NO_ERROR);
	
}





int DoDimensionsParm (char *parmName, char *tkn)

{

	if (expecting == Expecting(PARAMETER))
		{
		expecting = Expecting(EQUALSIGN);
		}
	else
		{
		/* set Ntax (numTaxa) *****************************************************************/
		if (!strcmp(parmName, "Ntax"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting(NUMBER);
			else if (expecting == Expecting(NUMBER))
				{
				sscanf (tkn, "%d", &numTaxa);
				defTaxa = YES;
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else
				return (ERROR);
			}
		/* set Nchar (numChar) ****************************************************************/
		else if (!strcmp(parmName, "Nchar"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting(NUMBER);
			else if (expecting == Expecting(NUMBER))
				{
				sscanf (tkn, "%d", &numChar);
				defChars = YES;
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else
				return (ERROR);
			}
		else
			return (ERROR);
		}

	return (NO_ERROR);

}





int DoDisclaimer (void)

{

	MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
	MrBayesPrint ("   Disclaimer                                                                    \n");
    MrBayesPrint ("                                                                                 \n");
	MrBayesPrint ("   Copyright 2003 by John P. Huelsenbeck and Fredrik Ronquist                    \n");
	MrBayesPrint ("                                                                                 \n");
	MrBayesPrint ("   This software package is provided \"as is\" and without a warranty of any     \n");
	MrBayesPrint ("   kind. In no event shall the authors be held responsible for any damage        \n");
	MrBayesPrint ("   resulting from the use of this software. The program--including source code,  \n");
	MrBayesPrint ("   example data sets, and executables--is distributed free of charge for         \n");
	MrBayesPrint ("   academic use only.                                                            \n");
	MrBayesPrint ("   ---------------------------------------------------------------------------   \n");

	return (NO_ERROR);
	
}





int DoEndBlock (void)

{

	if (inMrbayesBlock == YES)
		{
		MrBayesPrint ("   Exiting MrBayes block\n");
		inMrbayesBlock = NO;
		}
	else if (inDataBlock == YES)
		{
		MrBayesPrint ("   Exiting data block\n");
		inDataBlock = NO;
		}
	else if (inForeignBlock == YES)
		{
		MrBayesPrint ("   Exiting foreign block\n");
		inForeignBlock = NO;
		}
	else
		{
		MrBayesPrint ("   Unknown \"end\" statement\n");
		return (ERROR);
		}

	return (NO_ERROR);

}





int DoExecute (void)

{

	int			c, i, rc, lineTerm, longestLineLength, nErrors;
	char		*s;
	FILE		*fp;
#				if defined (MPI_ENABLED)
	int			sumErrors;
#				endif
		
	nErrors = 0;
	numOpenExeFiles++;
	s = NULL;
	
	if (numOpenExeFiles > 1)
		MrBayesPrint ("\n%s   Executing file \"%s\"\n", spacer, inputFileName);
	else
		MrBayesPrint ("%s   Executing file \"%s\"\n", spacer, inputFileName);

	/* open binary file */
	if ((fp = OpenBinaryFileR(inputFileName)) == NULL)
		nErrors++;

	/* set indentation to 0 */
	strcpy (spacer, "");
	
#	if defined (MPI_ENABLED)
	MPI_Allreduce (&nErrors, &sumErrors, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	if (sumErrors > 0)
		{
		MrBayesPrint ("%s   There was an error on at least one processor\n", spacer);
		goto errorExit;
		}
#	else
	if (nErrors > 0)
		goto errorExit;
#	endif

	/* find out what type of line termination is used */
	lineTerm = LineTermType (fp);
	if (lineTerm == LINETERM_MAC)
		MrBayesPrint ("%s   Macintosh line termination\n", spacer);
	else if (lineTerm == LINETERM_DOS)
		MrBayesPrint ("%s   DOS line termination\n", spacer);
	else if (lineTerm == LINETERM_UNIX)
		MrBayesPrint ("%s   UNIX line termination\n", spacer);
	else
		{
		MrBayesPrint ("%s   Unknown line termination\n", spacer);
		nErrors++;
		}
#	if defined (MPI_ENABLED)
	MPI_Allreduce (&nErrors, &sumErrors, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	if (sumErrors > 0)
		{
		MrBayesPrint ("%s   There was an error on at least one processor\n", spacer);
		goto errorExit;
		}
#	else
	if (nErrors > 0)
		goto errorExit;
#	endif
			
	/* find length of longest line */
	longestLineLength = LongestLine (fp);
	MrBayesPrint ("%s   Longest line length = %d\n", spacer, longestLineLength);
	longestLineLength += 50;
	
	/* check that longest line is not longer than CMD_STRING_LENGTH */
	if (longestLineLength >= CMD_STRING_LENGTH - 100)
		{
		MrBayesPrint ("%s   A maximum of %d characters is allowed on a single line\n", spacer, CMD_STRING_LENGTH - 100);
		MrBayesPrint ("%s   in a file. The longest line of the file %s\n", spacer, inputFileName);
		MrBayesPrint ("%s   contains at least one line with %d characters.\n", spacer, longestLineLength);
		nErrors++;
		}
#	if defined (MPI_ENABLED)
	MPI_Allreduce (&nErrors, &sumErrors, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	if (sumErrors > 0)
		{
		MrBayesPrint ("%s   There was an error on at least one processor\n", spacer);
		goto errorExit;
		}
#	else
	if (nErrors > 0)
		goto errorExit;
#	endif
	
	/* allocate a string long enough to hold a line */
	s = (char *)SafeMalloc((size_t) (longestLineLength * sizeof(char)));
	if (!s)
		{
		MrBayesPrint ("%s   Problem allocating string for reading file\n", spacer);
		nErrors++;
		}
#	if defined (MPI_ENABLED)
	MPI_Allreduce (&nErrors, &sumErrors, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	if (sumErrors > 0)
		{
		MrBayesPrint ("%s   There was an error on at least one processor\n", spacer);
		goto errorExit;
		}
#	else
	if (nErrors > 0)
		goto errorExit;
#	endif

	/* close binary file */
	SafeFclose (&fp);
	
	/* open text file */
	if ((fp = OpenTextFileR(inputFileName)) == NULL)
		{
		nErrors++;
		}
#	if defined (MPI_ENABLED)
	MPI_Allreduce (&nErrors, &sumErrors, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	if (sumErrors > 0)
		{
		MrBayesPrint ("%s   There was an error on at least one processor\n", spacer);
		goto errorExit;
		}
#	else
	if (nErrors > 0)
		goto errorExit;
#	endif
	
	/* parse file, reading each line in turn */
	MrBayesPrint ("%s   Parsing file\n", spacer);

	inMrbayesBlock = inDataBlock = inForeignBlock = inSumtBlock = NO;
	foundNewLine = NO;
	expecting = Expecting(COMMAND);

	do {
		/* read in a new line into s */
		i = 0;
		do {
			c = fgetc(fp);
			if (c == '\r' || c == '\n' || c == EOF)
				s[i++] = '\n';
			else
				s[i++] = c;
			} while (s[i-1] != '\n');
		s[i] = '\0';
		foundNewLine = YES;
		
		/* process string if not empty */
		if (strlen(s) > 1)
			{
			/* check that all characters in the string are valid */
			if (CheckStringValidity (s) == ERROR)
				{
				nErrors++;
				}
#			if defined (MPI_ENABLED)
			MPI_Allreduce (&nErrors, &sumErrors, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
			if (sumErrors > 0)
				{
				MrBayesPrint ("%s   There was an error on at least one processor\n", spacer);
				goto errorExit;
				}
#			else
			if (nErrors > 0)
				goto errorExit;
#			endif
				
			/* interpret commands on line */
			rc = ParseCommand (s);
			if (rc == ERROR)
				nErrors++;
#			if defined (MPI_ENABLED)
			MPI_Allreduce (&nErrors, &sumErrors, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
			if (sumErrors > 0)
				{
				MrBayesPrint ("%s   There was an error on at least one processor\n", spacer);
				goto errorExit;
				}
#			else
			if (nErrors > 0)
				goto errorExit;
#			endif
			if (rc == NO_ERROR_QUIT)
				nErrors++;
#			if defined (MPI_ENABLED)
			MPI_Allreduce (&nErrors, &sumErrors, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
			if (sumErrors > 0)
				goto quitExit;
#			else
			if (nErrors > 0)
				goto quitExit;
#			endif
			}
		} while (c != EOF); 
	
	MrBayesPrint ("%sReached end of file\n", spacer);

	if (inComment == YES)
		nErrors++;

#	if defined (MPI_ENABLED)
	MPI_Allreduce (&nErrors, &sumErrors, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	if (sumErrors > 0)
		{
		MrBayesPrint ("%s   There was an error on at least one processor\n", spacer);
		goto errorExit;
		}
#	else
	if (nErrors > 0)
		goto errorExit;
#	endif

	if (s)
		free (s);
	SafeFclose (&fp);
	numOpenExeFiles--;

	if (numOpenExeFiles > 0)
		{
		inMrbayesBlock = YES;
		MrBayesPrint ("\n   Returning execution to calling file ...\n");
		strcpy (spacer, "   ");
		}
	else
		strcpy (spacer, "");

	return (NO_ERROR);
	
	quitExit:
		if (s)
			free (s);
		SafeFclose (&fp);
		numOpenExeFiles--;
		if (numOpenExeFiles > 0)
			{
			inMrbayesBlock = YES;
			strcpy (spacer, "   ");
			}
		else
			strcpy (spacer, "");
		return (NO_ERROR_QUIT);
			
	errorExit:
		if (s)
			free (s);
		SafeFclose (&fp);
		numOpenExeFiles--;	/* we increase the value above even if no file is successfully opened */
		if (inComment == YES)
			{
			MrBayesPrint ("%s   ERROR: Reached end of file while in comment.\n", spacer);
			inComment = NO;
			numComments = 0;
			}
		if (numOpenExeFiles > 0)
			{
			inMrbayesBlock = YES;
			MrBayesPrint ("\n   Returning execution to calling file ...\n");
			strcpy (spacer, "   ");
			return (NO_ERROR);
			}
		else
			strcpy (spacer, "");

		strcpy (token, "Execute");
		i = 0;
		if (FindValidCommand (token, &i) == ERROR)
			MrBayesPrint ("%s   Could not find execute\n", spacer);

		return (ERROR);	
	
}





int DoExecuteParm (char *parmName, char *tkn)

{
	
	strcpy (inputFileName, tkn);
	
	expecting = Expecting (SEMICOLON);

	return (NO_ERROR);
	MrBayesPrint ("%s", parmName); /* just because I am tired of seeing the unused parameter error msg */

}





int DoExclude (void)

{

	int			i, alreadyDone;

	MrBayesPrint ("%s   Excluding character(s)\n", spacer);

	/* add set to tempSet */
	if (fromI >= 0 && toJ < 0)
		{
		if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
			return (ERROR);
		}
	else if (fromI >= 0 && toJ >= 0 && everyK < 0)
		{
		if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
			return (ERROR);
		}
	else if (fromI >= 0 && toJ >= 0 && everyK >= 0)
		{
		if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
			return (ERROR);
		}
		
	/* merge tempSet with charset */
	alreadyDone = NO;
	for (i=0; i<numChar; i++)
		{
		if (tempSet[i] == 1)
			{
			if (charInfo[i].isExcluded == YES && alreadyDone == NO)
				{
				MrBayesPrint ("%s   Some characters already excluded\n", spacer);
				alreadyDone = YES;
				}
			charInfo[i].isExcluded = YES;
			}
		}
		
	foundFirst = NO;

	return (NO_ERROR);
	
}





int DoExcludeParm (char *parmName, char *tkn)

{

	int		i, howMany, tempInt;
		
	if (defMatrix == NO)
		{
		MrBayesPrint ("%s   A matrix must be specified before you can exclude characters\n", spacer);
		return (ERROR);
		}
		
	if (foundFirst == NO)
		{
		/* this is the first time in */
		fromI = toJ = everyK = -1;
		foundDash = foundSlash = NO;
		for (i=0; i<numChar; i++) /* clear tempSet */
			tempSet[i] = 0;
		foundFirst = YES;
		}

	if (expecting == Expecting(ALPHA))
		{
		if (IsSame ("All", tkn) == SAME || IsSame ("All", tkn) == CONSISTENT_WITH)
			{
			for (i=0; i<numChar; i++)
				tempSet[i] = 1;
			}
		else if (IsSame ("Missambig", tkn) == SAME || IsSame ("Missambig", tkn) == CONSISTENT_WITH)
			{
			for (i=0; i<numChar; i++)
				{
				if (charInfo[i].isMissAmbig == YES)
					tempSet[i] = 1;
				}
			}
		else
			{
			/* we are using a pre-defined character set */
			if (numCharSets < 1)
				{
				MrBayesPrint ("%s   Could not find a character set called %s\n", spacer, tkn);
				return (ERROR);
				}
			if (CheckString (tkn, charSetNames, &howMany) == ERROR)
				{
				MrBayesPrint ("%s   Could not find a character set called %s\n", spacer, tkn);
				return (ERROR);
				}
			/* add characters from charset tkn to new tempSet */
			for (i=0; i<numChar; i++)
				{
				if (charInfo[i].charSet[howMany-1] == 1)
					tempSet[i] = 1;
				}
			fromI = toJ = everyK = -1;
			}

		expecting  = Expecting(ALPHA);
		expecting |= Expecting(NUMBER);
		expecting |= Expecting(SEMICOLON);
		}
	else if (expecting == Expecting(NUMBER))
		{
		if (strlen(tkn) == 1 && tkn[0] == '.')
			tempInt = numChar;
		else
			sscanf (tkn, "%d", &tempInt);
		if (tempInt <= 0 || tempInt > numChar)
			{
			MrBayesPrint ("%s   Character number %d is out of range (should be between %d and %d)\n", spacer, tempInt, 1, numChar);
			return (ERROR);
			}
		tempInt--;
		if (foundDash == YES)
			{
			if (fromI >= 0)
				toJ = tempInt;
			else
				{
				MrBayesPrint ("%s   Improperly formatted exclude set\n", spacer);
				return (ERROR);
				}
			foundDash = NO;
			}
		else if (foundSlash == YES)
			{
			tempInt++;
			if (tempInt <= 1)
				{
				MrBayesPrint ("%s   Improperly formatted exclude set\n", spacer);
				return (ERROR);
				}
			if (fromI >= 0 && toJ >= 0 && fromI < toJ)
				everyK = tempInt;
			else
				{
				MrBayesPrint ("%s   Improperly formatted exclude set\n", spacer);
				return (ERROR);
				}
			foundSlash = NO;
			}
		else
			{
			if (fromI >= 0 && toJ < 0)
				{
				if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
					return (ERROR);
				fromI = tempInt;
				}
			else if (fromI < 0 && toJ < 0)
				{
				fromI = tempInt;
				}
			else if (fromI >= 0 && toJ >= 0 && everyK < 0)
				{
				if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
					return (ERROR);
				fromI = tempInt;
				toJ = everyK = -1;
				}
			else if (fromI >= 0 && toJ >= 0 && everyK >= 0)
				{
				if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
					return (ERROR);
				fromI = tempInt;
				toJ = everyK = -1;
				}
			else
				{
				MrBayesPrint ("%s   Improperly formatted exclude set\n", spacer);
					{
					return (ERROR);
					}
				}
			}
		expecting  = Expecting(ALPHA);
		expecting |= Expecting(NUMBER);
		expecting |= Expecting(SEMICOLON);
		expecting |= Expecting(DASH);
		expecting |= Expecting(BACKSLASH);
		}
	else if (expecting == Expecting(DASH))
		{
		foundDash = YES;
		expecting = Expecting(NUMBER);
		}
	else if (expecting == Expecting(BACKSLASH))
		{
		foundSlash = YES;
		expecting = Expecting(NUMBER);
		}
	else
		return (ERROR);

	return (NO_ERROR);
	MrBayesPrint ("%s", parmName); /* just because I am tired of seeing the unused parameter error msg */

}





int DoFormat (void)

{

	int		howMany, nDivs;

	if (inDataBlock == NO)
		{
		MrBayesPrint ("%s   Formats can only be defined in a data block\n", spacer);
		return (ERROR);
		}
	if (CheckInitialPartitions() == ERROR)
		return (ERROR);
	/* enter name of default partition */
	if (AddToString ("Default", partitionNames, &howMany) == ERROR)
		{
		MrBayesPrint ("%s   Problem adding Default name to partition list\n", spacer);
		return (ERROR);
		}
	numDefinedPartitions = 1;
	nDivs = SetPartitionInfo (0);
	numCurrentDivisions = nDivs;
	if (nDivs == 1)
		MrBayesPrint ("%s   Setting default partition (does not divide up characters).\n", spacer);
	else
		MrBayesPrint ("%s   Setting default partition, dividing characters into %d parts.\n", spacer, nDivs);
	
	return (NO_ERROR);
	
}





int DoFormatParm (char *parmName, char *tkn)

{

	int			i, tempInt;
	char		tempStr[100];
	
	if (inDataBlock == NO)
		{
		MrBayesPrint ("%s   Formats can only be defined in a data block\n", spacer);
		return (ERROR);
		}
	if (defTaxa == NO || defChars == NO)
		{
		MrBayesPrint ("%s   The dimensions of the matrix must be defined before the format\n", spacer);
		return (ERROR);
		}
	
	if (expecting == Expecting(PARAMETER))
		{
		expecting = Expecting(EQUALSIGN);
		if (!strcmp(parmName, "Interleave"))
			{
			expecting = Expecting(EQUALSIGN) | Expecting(PARAMETER) | Expecting(SEMICOLON);
			isInterleaved = YES;
			}
		}
	else
		{
		/* set Datatype (dataType) ************************************************************/
		if (!strcmp(parmName, "Datatype"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting(ALPHA);
			else if (expecting == Expecting(ALPHA))
				{
				if (IsArgValid(tkn, tempStr) == NO_ERROR)
					{
					if (isMixed == NO)
						{
						if (!strcmp(tempStr, "Dna"))
							dataType = DNA;
						else if (!strcmp(tempStr, "Rna"))
							dataType = RNA;
						else if (!strcmp(tempStr, "Protein"))
							dataType = PROTEIN;
						else if (!strcmp(tempStr, "Restriction"))
							dataType = RESTRICTION;
						else if (!strcmp(tempStr, "Standard"))
							dataType = STANDARD;
						else if (!strcmp(tempStr, "Continuous"))
							dataType = CONTINUOUS;
						else if (!strcmp(tempStr, "Mixed"))
							{
							dataType = MIXED;
							isMixed = YES;
							for (i=0; i<numChar; i++)
								tempSet[i] = 0;
							fromI = toJ = everyK = -1;
							foundDash = foundSlash = NO;
							numPartitions = 0;
							MrBayesPrint ("%s   Data is Mixed\n", spacer);
							}
						if (dataType == MIXED)
							expecting = Expecting(LEFTPAR);
						else
							expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
						}
					else
						{
						if (!strcmp(tempStr, "Dna"))
							dataType = DNA;
						else if (!strcmp(tempStr, "Rna"))
							dataType = RNA;
						else if (!strcmp(tempStr, "Protein"))
							dataType = PROTEIN;
						else if (!strcmp(tempStr, "Restriction"))
							dataType = RESTRICTION;
						else if (!strcmp(tempStr, "Standard"))
							dataType = STANDARD;
						else if (!strcmp(tempStr, "Continuous"))
							dataType = CONTINUOUS;
						else if (!strcmp(tempStr, "Mixed"))
							{
							MrBayesPrint ("%s   Cannot have mixed datatype within a mixed datatype\n", spacer);
							return (ERROR);
							}
						expecting = Expecting(COLON);
						for (i=0; i<numChar; i++)
							tempSet[i] = 0;
						fromI = toJ = everyK = -1;
						foundDash = foundSlash = NO;
						}
					if (isMixed == NO)
						{
						numPartitions = 1;
						for (i=0; i<numChar; i++)
							{
							charInfo[i].charType = dataType;
							charInfo[i].partitionId[0] = numPartitions;
							}
						}
					}
				else
					{
					MrBayesPrint ("%s   Invalid data type argument\n", spacer);
					return (ERROR);
					}
				if (isMixed == NO)
					MrBayesPrint ("%s   Data is %s\n", spacer, tempStr);
				else if (strcmp(tempStr, "Mixed"))
					MrBayesPrint ("%s      Data for partition %d is %s\n", spacer, numPartitions+1, tempStr);
				}
			else if (expecting == Expecting(LEFTPAR))
				{
				expecting = Expecting(ALPHA);
				}
			else if (expecting == Expecting(COLON))
				{
				expecting = Expecting(NUMBER);
				}
			else if (expecting == Expecting(NUMBER))
				{
				if (strlen(tkn) == 1 && tkn[0] == '.')
					tempInt = numChar;
				else
					sscanf (tkn, "%d", &tempInt);
				if (tempInt <= 0 || tempInt > numChar)
					{
					MrBayesPrint ("%s   Character number %d is out of range (should be between %d and %d)\n", spacer, tempInt, 1, numChar);
					return (ERROR);
					}
				tempInt--;
				if (foundDash == YES)
					{
					if (fromI >= 0)
						toJ = tempInt;
					else
						{
						MrBayesPrint ("%s   Improperly formatted partition\n", spacer);
						return (ERROR);
						}
					foundDash = NO;
					}
				else if (foundSlash == YES)
					{
					tempInt++;
					if (tempInt <= 1)
						{
						MrBayesPrint ("%s   Improperly formatted partition\n", spacer);
						return (ERROR);
						}
					if (fromI >= 0 && toJ >= 0 && fromI < toJ)
						everyK = tempInt;
					else
						{
						MrBayesPrint ("%s   Improperly formatted partition\n", spacer);
						return (ERROR);
						}
					foundSlash = NO;
					}
				else
					{
					if (fromI >= 0 && toJ < 0)
						{
						if (AddToSet (fromI, toJ, everyK, numPartitions+1) == ERROR)
							return (ERROR);
						fromI = tempInt;
						}
					else if (fromI < 0 && toJ < 0)
						{
						fromI = tempInt;
						}
					else if (fromI >= 0 && toJ >= 0 && everyK < 0)
						{
						if (AddToSet (fromI, toJ, everyK, numPartitions+1) == ERROR)
							return (ERROR);
						fromI = tempInt;
						toJ = everyK = -1;
						}
					else if (fromI >= 0 && toJ >= 0 && everyK >= 0)
						{
						if (AddToSet (fromI, toJ, everyK, numPartitions+1) == ERROR)
							return (ERROR);
						fromI = tempInt;
						toJ = everyK = -1;
						}
					else
						{
						MrBayesPrint ("%s   Improperly formatted partition\n", spacer);
							{
							return (ERROR);
							}
						}
						
					}
				expecting  = Expecting(NUMBER);
				expecting |= Expecting(DASH);
				expecting |= Expecting(BACKSLASH);
				expecting |= Expecting(COMMA);
				expecting |= Expecting(RIGHTPAR);
				}
			else if (expecting == Expecting(DASH))
				{
				foundDash = YES;
				expecting = Expecting(NUMBER);
				}
			else if (expecting == Expecting(BACKSLASH))
				{
				foundSlash = YES;
				expecting = Expecting(NUMBER);
				}
			else if (expecting == Expecting(COMMA))
				{
				/* add set to tempSet */
				if (fromI >= 0 && toJ < 0)
					{
					if (AddToSet (fromI, toJ, everyK, numPartitions+1) == ERROR)
						{
						if (RemoveLastFromString (charSetNames) == ERROR)
							return (ERROR);
						return (ERROR);
						}
					}
				else if (fromI >= 0 && toJ >= 0 && everyK < 0)
					{
					if (AddToSet (fromI, toJ, everyK, numPartitions+1) == ERROR)
						{
						if (RemoveLastFromString (charSetNames) == ERROR)
							return (ERROR);
						return (ERROR);
						}
					}
				else if (fromI >= 0 && toJ >= 0 && everyK >= 0)
					{
					if (AddToSet (fromI, toJ, everyK, numPartitions+1) == ERROR)
						{
						if (RemoveLastFromString (charSetNames) == ERROR)
							return (ERROR);
						return (ERROR);
						}
					}
				for (i=0; i<numChar; i++)
					{
					if (tempSet[i] == numPartitions)
						charInfo[i].charType = dataType;
					}

				/* merge tempSet */
				for (i=0; i<numChar; i++)
					{
					if (tempSet[i] != 0)
						{
						if (charInfo[i].partitionId[0] == 0)
							{
							charInfo[i].charType = dataType;
							charInfo[i].partitionId[0] = numPartitions + 1;
							}
						else
							{
							MrBayesPrint ("%s   Improperly formatted partition (same character found in multiple partitions)\n", spacer);
							return (ERROR);
							}
						}
					}
				
				/* increment number of partitions */
				numPartitions++;				
				expecting = Expecting(ALPHA);
				}
			else if (expecting == Expecting(RIGHTPAR))
				{
				/* add set to tempSet */
				if (fromI >= 0 && toJ < 0)
					{
					if (AddToSet (fromI, toJ, everyK, numPartitions+1) == ERROR)
						{
						if (RemoveLastFromString (charSetNames) == ERROR)
							return (ERROR);
						return (ERROR);
						}
					}
				else if (fromI >= 0 && toJ >= 0 && everyK < 0)
					{
					if (AddToSet (fromI, toJ, everyK, numPartitions+1) == ERROR)
						{
						if (RemoveLastFromString (charSetNames) == ERROR)
							return (ERROR);
						return (ERROR);
						}
					}
				else if (fromI >= 0 && toJ >= 0 && everyK >= 0)
					{
					if (AddToSet (fromI, toJ, everyK, numPartitions+1) == ERROR)
						{
						if (RemoveLastFromString (charSetNames) == ERROR)
							return (ERROR);
						return (ERROR);
						}
					}
					
				/* merge tempSet */
				for (i=0; i<numChar; i++)
					{
					if (tempSet[i] != 0)
						{
						if (charInfo[i].partitionId[0] == 0)
							{
							charInfo[i].charType = dataType;
							charInfo[i].partitionId[0] = numPartitions + 1;
							}
						else
							{
							MrBayesPrint ("%s   Improperly formatted partition (same character found in multiple partitions)\n", spacer);
							return (ERROR);
							}
						}
					}
				
				/* increment number of partitions */
				numPartitions++;				
				if (isMixed == YES)
					dataType = MIXED;
					
				if (numPartitions > 1)
					MrBayesPrint ("%s   There are a total of %d default data divisions\n", spacer, numPartitions);
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else
				return (ERROR);
			}
		/* set Interleave (isInterleaved) *****************************************************/
		else if (!strcmp(parmName, "Interleave"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting(ALPHA);
			else if (expecting == Expecting(ALPHA))
				{
				if (IsArgValid(tkn, tempStr) == NO_ERROR)
					{
					if (!strcmp(tempStr, "Yes"))
						isInterleaved = YES;
					else
						isInterleaved = NO;
					}
				else
					{
					MrBayesPrint ("%s   Invalid argument for interleaved data\n", spacer);
					return (ERROR);
					}
				if (isInterleaved == YES)
					MrBayesPrint ("%s   Data matrix is interleaved\n", spacer);
				else
					MrBayesPrint ("%s   Data matrix is not interleaved\n", spacer);
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else
				return (ERROR);
			}
		/* set Gap (gapId) ********************************************************************/
		else if (!strcmp(parmName, "Gap"))
			{
			if (expecting == Expecting(EQUALSIGN))
				{
				expecting  = Expecting(ALPHA);
				expecting |= Expecting(QUESTIONMARK);
				expecting |= Expecting(DASH);
				expecting |= Expecting(NUMBER);
				expecting |= Expecting(ASTERISK);
				expecting |= Expecting(EXCLAMATIONMARK);
				expecting |= Expecting(PERCENT);
				expecting |= Expecting(WEIRD);
				}
			else if (((expecting & Expecting(ALPHA)) == Expecting(ALPHA)) || 
				 	 ((expecting & Expecting(QUESTIONMARK)) == Expecting(QUESTIONMARK)) || 
					 ((expecting & Expecting(DASH)) == Expecting(DASH)) || 
					 ((expecting & Expecting(NUMBER)) == Expecting(NUMBER)) || 
					 ((expecting & Expecting(ASTERISK)) == Expecting(ASTERISK)) || 
					 ((expecting & Expecting(EXCLAMATIONMARK)) == Expecting(EXCLAMATIONMARK)) || 
					 ((expecting & Expecting(PERCENT)) == Expecting(PERCENT)) || 
					 ((expecting & Expecting(WEIRD)) == Expecting(WEIRD)))
				{
				if (strlen(tkn) == 1)
					{
					if (tkn[0] == matchId || tkn[0] == missingId)
						{
						MrBayesPrint ("%s   Gap character matches matching or missing characters\n", spacer);
						return (ERROR);
						}
					gapId = tkn[0];
					}
				else
					{
					MrBayesPrint ("%s   Invalid gap argument %s\n", spacer, tkn);
					return (ERROR);
					}
				MrBayesPrint ("%s   Gaps coded as %s\n", spacer, tkn);
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else
				return (ERROR);
			}
		/* set Missing (missingId) ************************************************************/
		else if (!strcmp(parmName, "Missing"))
			{
			if (expecting == Expecting(EQUALSIGN))
				{
				expecting  = Expecting(ALPHA);
				expecting |= Expecting(QUESTIONMARK);
				expecting |= Expecting(DASH);
				expecting |= Expecting(NUMBER);
				expecting |= Expecting(ASTERISK);
				expecting |= Expecting(EXCLAMATIONMARK);
				expecting |= Expecting(PERCENT);
				expecting |= Expecting(WEIRD);
				}
			else if (((expecting & Expecting(ALPHA)) == Expecting(ALPHA)) || 
				 	 ((expecting & Expecting(QUESTIONMARK)) == Expecting(QUESTIONMARK)) || 
					 ((expecting & Expecting(DASH)) == Expecting(DASH)) || 
					 ((expecting & Expecting(NUMBER)) == Expecting(NUMBER)) || 
					 ((expecting & Expecting(ASTERISK)) == Expecting(ASTERISK)) || 
					 ((expecting & Expecting(EXCLAMATIONMARK)) == Expecting(EXCLAMATIONMARK)) || 
					 ((expecting & Expecting(PERCENT)) == Expecting(PERCENT)) || 
					 ((expecting & Expecting(WEIRD)) == Expecting(WEIRD)))
				{
				if (strlen(tkn) == 1)
					{
					if (tkn[0] == gapId || tkn[0] == matchId)
						{
						MrBayesPrint ("%s   Missing character matches matching or gap characters\n", spacer);
						return (ERROR);
						}
					missingId = tkn[0];
					}
				else
					{
					MrBayesPrint ("%s   Invalid missing argument %s\n", spacer, tkn);
					return (ERROR);
					}
				MrBayesPrint ("%s   Missing data coded as %s\n", spacer, tkn);
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else
				return (ERROR);
			}
		/* set Matchchar (matchId) ************************************************************/
		else if (!strcmp(parmName, "Matchchar"))
			{
			if (expecting == Expecting(EQUALSIGN))
				{
				expecting  = Expecting(ALPHA);
				expecting |= Expecting(QUESTIONMARK);
				expecting |= Expecting(DASH);
				expecting |= Expecting(NUMBER);
				expecting |= Expecting(ASTERISK);
				expecting |= Expecting(EXCLAMATIONMARK);
				expecting |= Expecting(PERCENT);
				expecting |= Expecting(WEIRD);
				}
			else if (((expecting & Expecting(ALPHA)) == Expecting(ALPHA)) || 
				 	 ((expecting & Expecting(QUESTIONMARK)) == Expecting(QUESTIONMARK)) || 
					 ((expecting & Expecting(DASH)) == Expecting(DASH)) || 
					 ((expecting & Expecting(NUMBER)) == Expecting(NUMBER)) || 
					 ((expecting & Expecting(ASTERISK)) == Expecting(ASTERISK)) || 
					 ((expecting & Expecting(EXCLAMATIONMARK)) == Expecting(EXCLAMATIONMARK)) || 
					 ((expecting & Expecting(PERCENT)) == Expecting(PERCENT)) || 
					 ((expecting & Expecting(WEIRD)) == Expecting(WEIRD)))
				{
				if (strlen(tkn) == 1)
					{
					if (tkn[0] == gapId || tkn[0] == missingId)
						{
						MrBayesPrint ("%s   Matching character matches gap or missing characters\n", spacer);
						return (ERROR);
						}
					matchId = tkn[0];
					}
				else
					{
					MrBayesPrint ("%s   Invalid matchchar argument %s\n", spacer, tkn);
					return (ERROR);
					}
				MrBayesPrint ("%s   Matching characters coded as %s\n", spacer, tkn);
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else
				return (ERROR);
			}
		else
			return (ERROR);
		}

	return (NO_ERROR);

}





int DoHelp (void)

{

	int			i, j, longestDescription;
	CmdType		*p;

	if (foundFirst == NO)
		{
		longestDescription = 0;
		for (i=1; i<NUMCOMMANDS; i++)
			{
			p = commands + i;
			if ((int)strlen(p->string) > longestDescription)
				longestDescription = (int) strlen(p->string);
			}
		
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		MrBayesPrint ("   Commands that are available from the command                                  \n");
		MrBayesPrint ("   line or from a MrBayes block include:                                         \n");
	    MrBayesPrint ("                                                                                 \n");
		for (i=1; i<NUMCOMMANDS; i++)
			{
			p = commands + i;
			if (p->cmdUse == IN_CMD && p->hiding == SHOW)
				{
				MrBayesPrint ("   %s", p->string);
				for (j=0; j<longestDescription - (int) strlen(p->string); j++)
				MrBayesPrint (" ");
				MrBayesPrint (" -- %s\n", p->cmdDescription);
				}
			}
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("   Commands that should be in a NEXUS file (data                                 \n");
		MrBayesPrint ("   block or trees block) include:                                                \n");
	    MrBayesPrint ("                                                                                 \n");
		for (i=1; i<NUMCOMMANDS; i++)
			{
			p = commands + i;
			if (p->cmdUse == IN_FILE && p->hiding == SHOW)
				{
				MrBayesPrint ("   %s", p->string);
				for (j=0; j<longestDescription - (int) strlen(p->string); j++)
				MrBayesPrint (" ");
				MrBayesPrint (" -- %s\n", p->cmdDescription);
				}
			}
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   Note that this program supports the use of the shortest unambiguous           \n"); 
	    MrBayesPrint ("   spelling of the above commands (e.g., \"exe\" instead of \"execute\").        \n");
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		}
	foundFirst = NO;

	return (NO_ERROR);
	
}





int DoHelpParm (char *parmName, char *tkn)

{
		
	int			i, j, tkLen, targetLen, numDiff, numMatches;
	CmdType		*p, *q=NULL;

	if (expecting == Expecting(ALPHA))
		{
		p = commands + 0;
		tkLen = (int) strlen(tkn);
		numMatches = 0;
		for (i=0; i<NUMCOMMANDS; i++)
			{
			targetLen = (int) strlen(p->string);
			if (tkLen <= targetLen)
				{
				for (j=0, numDiff=0; j<tkLen; j++)
					{
					if (ChangeCase(tkn[j]) != ChangeCase(p->string[j]))
						numDiff++;
					}
				if (numDiff == 0)
					{
					numMatches++;
					q = p;
					if (tkLen == targetLen)
						break;
					}
				}		
			p++;
			}
		if (numMatches == 0)
			{
			MrBayesPrint ("%s   Could not find command \"%s\"\n", spacer, tkn);
			return (ERROR);
			}
		else if (numMatches == 1)
			{
			if (GetUserHelp (q->string) == ERROR)
				{
				MrBayesPrint ("%s   Problem getting help for command \"%s\"\n", spacer, q->string);
				}
			}
		else 
			{
			MrBayesPrint ("%s   Ambiguous command \"%s\"\n", spacer, tkn);
			return (ERROR);
			}
			
		expecting = Expecting(SEMICOLON);
		foundFirst = YES;
		}
	else
		return (ERROR);

	return (NO_ERROR);
	MrBayesPrint ("%s", parmName); /* just because I am tired of seeing the unused parameter error msg */

}





int DoInclude (void)

{

	int			i, alreadyDone;

	MrBayesPrint ("%s   Including character(s)\n", spacer);

	/* add set to tempSet */
	if (fromI >= 0 && toJ < 0)
		{
		if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
			return (ERROR);
		}
	else if (fromI >= 0 && toJ >= 0 && everyK < 0)
		{
		if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
			return (ERROR);
		}
	else if (fromI >= 0 && toJ >= 0 && everyK >= 0)
		{
		if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
			return (ERROR);
		}
		
	/* merge tempSet with excludedChars */
	alreadyDone = NO;
	for (i=0; i<numChar; i++)
		{
		if (tempSet[i] == 1)
			{
			if (charInfo[i].isExcluded == NO && alreadyDone == NO)	
				{
				MrBayesPrint ("%s   Some characters already included\n", spacer);
				alreadyDone = YES;
				}
			charInfo[i].isExcluded = NO;
			}
		}

	return (NO_ERROR);
	
}





int DoIncludeParm (char *parmName, char *tkn)

{

	int		i, howMany, tempInt;
		
	if (defMatrix == NO)
		{
		MrBayesPrint ("%s   A matrix must be specified before you can include characters\n", spacer);
		return (ERROR);
		}
		
	if (foundFirst == NO)
		{
		/* this is the first time in */
		fromI = toJ = everyK = -1;
		foundDash = foundSlash = NO;
		for (i=0; i<numChar; i++) /* clear tempSet */
			tempSet[i] = 0;
		foundFirst = YES;
		}

	if (expecting == Expecting(ALPHA))
		{
		if (IsSame ("All", tkn) == SAME || IsSame ("All", tkn) == CONSISTENT_WITH)
			{
			for (i=0; i<numChar; i++)
				tempSet[i] = 1;
			}
		else if (IsSame ("Missambig", tkn) == SAME || IsSame ("Missambig", tkn) == CONSISTENT_WITH)
			{
			for (i=0; i<numChar; i++)
				{
				if (charInfo[i].isMissAmbig == YES)
					tempSet[i] = 1;
				}
			}
		else
			{
			/* we are using a pre-defined character set */
			if (numCharSets < 1)
				{
				MrBayesPrint ("%s   Could not find a character set called %s\n", spacer, tkn);
				return (ERROR);
				}
			if (CheckString (tkn, charSetNames, &howMany) == ERROR)
				{
				MrBayesPrint ("%s   Could not find a character set called %s\n", spacer, tkn);
				return (ERROR);
				}
			/* add characters from charset tkn to new tempSet */
			for (i=0; i<numChar; i++)
				{
				if (charInfo[i].charSet[howMany-1] == 1)
					tempSet[i] = 1;
				}
			fromI = toJ = everyK = -1;
			}

		expecting  = Expecting(ALPHA);
		expecting |= Expecting(NUMBER);
		expecting |= Expecting(SEMICOLON);
		}
	else if (expecting == Expecting(NUMBER))
		{
		if (strlen(tkn) == 1 && tkn[0] == '.')
			tempInt = numChar;
		else
			sscanf (tkn, "%d", &tempInt);
		if (tempInt <= 0 || tempInt > numChar)
			{
			MrBayesPrint ("%s   Character number %d is out of range (should be between %d and %d)\n", spacer, tempInt, 1, numChar);
			return (ERROR);
			}
		tempInt--;
		if (foundDash == YES)
			{
			if (fromI >= 0)
				toJ = tempInt;
			else
				{
				MrBayesPrint ("%s   Improperly formatted include set\n", spacer);
				return (ERROR);
				}
			foundDash = NO;
			}
		else if (foundSlash == YES)
			{
			tempInt++;
			if (tempInt <= 1)
				{
				MrBayesPrint ("%s   Improperly formatted include set\n", spacer);
				return (ERROR);
				}
			if (fromI >= 0 && toJ >= 0 && fromI < toJ)
				everyK = tempInt;
			else
				{
				MrBayesPrint ("%s   Improperly formatted include set\n", spacer);
				return (ERROR);
				}
			foundSlash = NO;
			}
		else
			{
			if (fromI >= 0 && toJ < 0)
				{
				if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
					return (ERROR);
				fromI = tempInt;
				}
			else if (fromI < 0 && toJ < 0)
				{
				fromI = tempInt;
				}
			else if (fromI >= 0 && toJ >= 0 && everyK < 0)
				{
				if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
					return (ERROR);
				fromI = tempInt;
				toJ = everyK = -1;
				}
			else if (fromI >= 0 && toJ >= 0 && everyK >= 0)
				{
				if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
					return (ERROR);
				fromI = tempInt;
				toJ = everyK = -1;
				}
			else
				{
				MrBayesPrint ("%s   Improperly formatted include set\n", spacer);
					{
					return (ERROR);
					}
				}
			}
		expecting  = Expecting(ALPHA);
		expecting |= Expecting(NUMBER);
		expecting |= Expecting(SEMICOLON);
		expecting |= Expecting(DASH);
		expecting |= Expecting(BACKSLASH);
		}
	else if (expecting == Expecting(DASH))
		{
		foundDash = YES;
		expecting = Expecting(NUMBER);
		}
	else if (expecting == Expecting(BACKSLASH))
		{
		foundSlash = YES;
		expecting = Expecting(NUMBER);
		}
	else
		return (ERROR);

	return (NO_ERROR);
	MrBayesPrint ("%s", parmName); /* just because I am tired of seeing the unused parameter error msg */

}





int DoLog (void)

{

	if (logToFile == YES)
		{
		SafeFclose (&logFileFp);
		if (replaceLogFile == YES)
			{
			if ((logFileFp = fopen (logFileName, "w")) == NULL)  
				{
				MrBayesPrint ("%s   Could not open log file \"%s\"\n", spacer, logFileName);
				logToFile = NO;
				return (ERROR);
				}
			}
		else
			{
			if ((logFileFp = fopen (logFileName, "a")) == NULL)  
				{
				MrBayesPrint ("%s   Could not open log file \"%s\"\n", spacer, logFileName);
				logToFile = NO;
				return (ERROR);
				}
			}
		MrBayesPrint ("%s   Logging screen output to file \"%s\"\n", spacer, logFileName);
		}
	else
		{
		SafeFclose (&logFileFp);
		MrBayesPrint ("%s   Terminating log output\n", spacer);
		}

	return (NO_ERROR);
	
}





int DoLogParm (char *parmName, char *tkn)

{
	
	if (expecting == Expecting(PARAMETER))
		{
		if (!strcmp(parmName, "Start"))
			{
			if (logToFile == YES)
				MrBayesPrint ("%s   Logging to file is already on\n", spacer, logFileName);
			else
				logToFile = YES;
			expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
			}
		else if (!strcmp(parmName, "Stop"))
			{
			if (logToFile == NO)
				MrBayesPrint ("%s   Logging to file is already off\n", spacer, logFileName);
			else
				logToFile = NO;
			expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
			}
		else if (!strcmp(parmName, "Replace"))
			{
			replaceLogFile = YES;
			expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
			}
		else if (!strcmp(parmName, "Append"))
			{
			replaceLogFile = NO;
			expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
			}
		else
			expecting = Expecting(EQUALSIGN);
		}
	else
		{
		if (!strcmp(parmName, "Filename"))
			{
			if (expecting == Expecting(EQUALSIGN))
				{
				expecting = Expecting(ALPHA);
				readWord = YES;
				}
			else if (expecting == Expecting(ALPHA))
				{
				strcpy (logFileName, tkn);
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else
				return (ERROR);
			}
		else
			{
			MrBayesPrint ("%s   Unknown parameter in Log\n", spacer);
			return (ERROR);
			}
		}

	return (NO_ERROR);
	
}





int DoManual (void)

{

	int		i, j, logSetting;
	char	title[100];
	FILE	*fp, *logfp;
	CmdType	*p;
	
	/* try to open file, return error if present */
	if ((fp = fopen(manFileName,"r")) != NULL)
		{
		MrBayesPrint ("%s   File \"%s\" already exists \n", spacer, manFileName);
		SafeFclose(&fp);
		return (ERROR);
		}

	/* try to open file for writing, return error if not possible */
	if ((fp = fopen(manFileName,"w")) == NULL)
		{
		MrBayesPrint ("%s   Could not open file for writing\n", spacer, manFileName);
		return (ERROR);
		}

	/* print message */
	MrBayesPrint ("%s   Producing command reference file \"%s\"\n", spacer, manFileName);

	/* temporarily disable normal logging and switch echoing off */
	logSetting = logToFile;
	logfp = logFileFp;
	echoMB = NO;
	logToFile = YES;
	logFileFp = fp;
	
	/* produce command reference file */
	/* header */
	strcpy (title, "Command Reference for MrBayes ver. ");
	strcat (title, VERSION_NUMBER);

	i = (70 - (int) strlen (title)) / 2;
	j = 70 - i - (int) strlen(title);

	MrBayesPrint ("                                                                                 \n");
	MrBayesPrint ("                                                                                 \n");
	MrBayesPrint ("                                                                                 \n");
	MrBayesPrint ("                                                                                 \n");
	MrBayesPrint ("      %*c%s%*c      \n", i, ' ', title, j, ' ');
	MrBayesPrint ("                                                                                 \n");
	MrBayesPrint ("                   (c) John P. Huelsenbeck and Fredrik Ronquist                  \n");
	MrBayesPrint ("                                                                                 \n");
	MrBayesPrint ("                                                                                 \n");

	/* summary */
	MrBayesPrint ("                                                                                 \n");
	MrBayesPrint ("   ***************************************************************************   \n");
	MrBayesPrint ("   *                                                                         *   \n");
	MrBayesPrint ("   *  1. Command summary                                                     *   \n");
	MrBayesPrint ("   *                                                                         *   \n");
	MrBayesPrint ("   ***************************************************************************   \n");
	MrBayesPrint ("                                                                                 \n");
	foundFirst = NO;
	if (DoHelp() == ERROR)
		{
		MrBayesPrint ("%s   Could not produce command reference summary\n", spacer);
		goto errorExit;
		}
	
	/* list of MrBayes commands */
	MrBayesPrint ("                                                                                 \n");
	MrBayesPrint ("   ***************************************************************************   \n");
	MrBayesPrint ("   *                                                                         *   \n");
	MrBayesPrint ("   *  2. MrBayes commands                                                    *   \n");
	MrBayesPrint ("   *                                                                         *   \n");
	MrBayesPrint ("   ***************************************************************************   \n");
	MrBayesPrint ("                                                                                 \n");
	for (i=1; i<NUMCOMMANDS; i++)
		{
		p = commands + i;
		if (p->cmdUse == IN_CMD && p->hiding == SHOW)
			{
			if (GetUserHelp(p->string)==ERROR)
				goto errorExit;
			}
		}

	/* list of data or tree block commands */
	MrBayesPrint ("                                                                                 \n");
	MrBayesPrint ("   ***************************************************************************   \n");
	MrBayesPrint ("   *                                                                         *   \n");
	MrBayesPrint ("   *  3. 'Data' or 'tree' block commands (in #NEXUS file)                    *   \n");
	MrBayesPrint ("   *                                                                         *   \n");
	MrBayesPrint ("   ***************************************************************************   \n");
	MrBayesPrint ("                                                                                 \n");
	for (i=1; i<NUMCOMMANDS; i++)
		{
		p = commands + i;
		if (p->cmdUse == IN_FILE && p->hiding == SHOW)
			{
			if (GetUserHelp(p->string) == ERROR)
				goto errorExit;
			}
		}

	/* return logging to previous setings and switch echoing on */
	SafeFclose (&fp);
	logToFile = logSetting;
	logFileFp = logfp;
	echoMB = YES;

	MrBayesPrint ("%s   Successfully produced command reference file \"%s\"\n", spacer, manFileName);

	return (NO_ERROR);

	errorExit:
		SafeFclose (&fp);
		logToFile = logSetting;
		logFileFp = logfp;
		echoMB = YES;

		return (ERROR);
		
}





int DoManualParm (char *parmName, char *tkn)

{
	
	if (expecting == Expecting(PARAMETER))
		{
		expecting = Expecting(EQUALSIGN);
		}
	else
		{
		if (!strcmp(parmName, "Filename"))
			{
			if (expecting == Expecting(EQUALSIGN))
				{
				expecting = Expecting(ALPHA);
				readWord = YES;
				}
			else if (expecting == Expecting(ALPHA))
				{
				strcpy (manFileName, tkn);
				expecting = Expecting(SEMICOLON);
				}
			else
				return (ERROR);
			}
		else
			{
			MrBayesPrint ("%s   Unknown parameter in Manual\n", spacer);
			return (ERROR);
			}
		}

	return (NO_ERROR);
	
}





int DoMatrix (void)

{

	int			i, j, hasMissingAmbig;
	
	if (taxonCount != numTaxa)
		{
		MrBayesPrint ("%s   Problem with number of taxa read in (%d taxa read in, while expecting %d)\n", spacer, taxonCount, numTaxa);
		FreeMatrix();
		return (ERROR);
		}
	for (i=0; i<numTaxa; i++)
		{
		if (taxaInfo[i].charCount != numChar)
			{
			MrBayesPrint ("%s   Problem with number of characters read in (%d expected for taxon %d, %d read in)\n", spacer, numChar, i, taxaInfo[i].charCount);
			FreeMatrix();
			return (ERROR);
			}
		}
		
	/* find out which characters have missing or ambiguous states (one time only, so no special function) */
	for (i=0; i<numChar; i++)
		{
		hasMissingAmbig = NO;
		for (j=0; j<numTaxa; j++)
			{
			if (IsMissing (matrix[pos(j,i,numChar)], charInfo[i].charType) == YES)
				hasMissingAmbig = YES;
			if (IsAmbig (matrix[pos(j,i,numChar)], charInfo[i].charType) == YES)
				hasMissingAmbig = YES;
			}
		if (hasMissingAmbig == YES)
			charInfo[i].isMissAmbig = YES;
		}
		
	/* set the type of patterns sampled */
	for (i=0; i<numChar; i++)
		{
		if (charInfo[i].charType == RESTRICTION)
			strcpy(modelParams[charInfo[i].partitionId[0]-1].coding, "Noabsencesites"); 
		else if (charInfo[i].charType == STANDARD)
			strcpy(modelParams[charInfo[i].partitionId[0]-1].coding, "Variable"); 
		}

	strcpy (sumtParams.sumtFileName, inputFileName);
	strcpy (sumpParams.sumpFileName, inputFileName);
	strcpy (sumpParams.sumpOutfile, sumpParams.sumpFileName);
	strcat (sumpParams.sumpOutfile, ".stat");

	if (chainParams.numRuns == 1)
		sprintf (comptreeParams.comptFileName1, "%s.t", inputFileName);
	else /* if (chainParams.numRuns > 1) */
		sprintf (comptreeParams.comptFileName1, "%s.run1.t", inputFileName);
	strcpy (comptreeParams.comptFileName2, comptreeParams.comptFileName1);

	if (chainParams.numRuns == 1)
		sprintf (plotParams.plotFileName, "%s.p", inputFileName);
	else /* if (chainParams.numRuns > 1) */
		sprintf (plotParams.plotFileName, "%s.run1.p", inputFileName);

	strcpy (chainParams.chainFileName, inputFileName);

	if (chainParams.numRuns > 1)
		MrBayesPrint ("%s   Setting output file names to \"%s.run<i>.<p/t>\"\n", spacer, chainParams.chainFileName);
	else
		MrBayesPrint ("%s   Setting output file names to \"%s.<p/t>\"\n", spacer, chainParams.chainFileName);

	MrBayesPrint ("%s   Successfully read matrix\n", spacer);
	if (matrixHasPoly == YES)
		MrBayesPrint ("%s   Matrix  contains polymorphisms, interpreted as ambiguity\n", spacer);
	defMatrix = YES;
	
#	if 0
	for (i=0; i<numChar; i++)
		{
		int		j;
		MrBayesPrint ("%4d -- ", i+1);
		for (j=0; j<numTaxa; j++)
			MrBayesPrint ("%2d ", matrix[pos(j,i,numChar)]);
		MrBayesPrint ("\n");
		}
#	endif
	return (NO_ERROR);

}





int DoMatrixParm (char *parmName, char *tkn)

{

	int				i, j, charCode, howMany;
	MrBFlt			charValue;

	expecting  = Expecting(ALPHA);
	expecting |= Expecting(QUESTIONMARK);
	expecting |= Expecting(DASH);
	expecting |= Expecting(NUMBER);
	expecting |= Expecting(ASTERISK);
	expecting |= Expecting(EXCLAMATIONMARK);
	expecting |= Expecting(PERCENT);
	expecting |= Expecting(WEIRD);
	expecting |= Expecting(SEMICOLON);
	expecting |= Expecting(LEFTPAR);
	expecting |= Expecting(RIGHTPAR);
	expecting |= Expecting(LEFTCURL);
	expecting |= Expecting(RIGHTCURL);

	if (defTaxa == NO || defChars == NO)
		{
		MrBayesPrint ("%s   Number of taxa and characters needs to be defined before matrix is read\n", spacer);
		goto errorExit;
		}
	if (inDataBlock == NO)
		{
		MrBayesPrint ("%s   Must be in data block to read in character matrix\n", spacer);
		goto errorExit;
		}

	if (isFirstMatrixRead == YES)
		{
		foundNewLine = YES;
		isFirstInterleavedBlock = YES;
		taxonCount = 0;
		isNegative = NO;
		}
	isFirstMatrixRead = NO;
	
	/* allow line breaks in non-interleaved matrices */
	if (isInterleaved == NO)
		{
		if (foundNewLine == YES && taxonCount > 0)
			{
			if (taxaInfo[taxonCount-1].charCount < numChar)
				foundNewLine = NO;
			}
		}

	if (taxonCount >= numTaxa && foundNewLine == YES)
		{
		if (isInterleaved == YES)
			{
			taxonCount = 0;
			isFirstInterleavedBlock = NO;
			}
		else
			{
			MrBayesPrint ("%s   Too many taxa in matrix\n", spacer);
			goto errorExit;
			}
		}
	
	if (taxaInfo[0].charCount > 4010)
		i = 1;

	if (foundNewLine == YES)
		{
		/* Should be a taxon. */
		if (isFirstInterleavedBlock == YES)
			{
			/* If this is the first interleaved block, then we need to add the 
			   taxon to the string of taxon names. */
			if (AddToString (tkn, taxaNames, &howMany) == ERROR)
				{
				MrBayesPrint ("%s   Problem adding taxon %s to list\n", spacer, tkn);
				goto errorExit;
				}
			if (howMany - 1 != taxonCount)
				{
				MrBayesPrint ("%s   Problem adding taxon %s to list\n", spacer, tkn);
				goto errorExit;
				}
			if (numTaxa < 10)
				MrBayesPrint ("%s   Taxon %d -> %s\n", spacer, taxonCount+1, tkn);
			else if (numTaxa < 100 && numTaxa >= 10)
				MrBayesPrint ("%s   Taxon %2d -> %s\n", spacer, taxonCount+1, tkn);
			else if (numTaxa < 1000 && numTaxa >= 100)
				MrBayesPrint ("%s   Taxon %3d -> %s\n", spacer, taxonCount+1, tkn);
			else
				MrBayesPrint ("%s   Taxon %4d -> %s\n", spacer, taxonCount+1, tkn);
			}
		else
			{
			/* If this is not the first interleaved block, then we need to
			   check to see if taxon name is present and in correct place. */
			if (CheckString (tkn, taxaNames, &howMany) == ERROR)
				{
				MrBayesPrint ("%s   Could not find taxon %s in list of taxa\n", spacer, tkn);
				goto errorExit;
				}
			if (howMany - 1 != taxonCount)
				{
				MrBayesPrint ("%s   Could not find taxon %s in correct position in list of taxa\n", spacer, tkn);
				goto errorExit;
				}
			}
		foundNewLine = NO;
		isNegative = NO;
		taxonCount++;
		}
	else
		{
		/* Should be a character (either continuous or otherwise). */
		if (charInfo[taxaInfo[taxonCount-1].charCount].charType == CONTINUOUS)
			{
			/* If we have a CONTINUOUS character, then the entire token should either be
			   a number or a dash (for a negative sign). */
			if (!strcmp(tkn, "-"))
				{
				/* Dealing with a negative number. We will multiply the next tkn, which
				   had better be a number, by -1. */
				isNegative = YES;
				}
			else
				{
				/* We have a number, we hope. */
				if (tkn[0] == matchId)
					{
					/* If the token is a matchchar, then things are simple. */
					if (taxonCount == 1)
						{
						MrBayesPrint ("%s   Matching characters cannot be in first taxon\n", spacer);
						goto errorExit;
						}
					charCode = matrix[pos(0,taxaInfo[taxonCount-1].charCount,numChar)];
					matrix[pos(taxonCount-1,taxaInfo[taxonCount-1].charCount,numChar)] = charCode;
					}
				else
					{
					/* Otherwise, we have a number. Check that it is a valid number first... */
					if (!IsIn(tkn[0],"0123456789."))
						{
						MrBayesPrint ("%s   Expecting a number for the continuous character\n", spacer);
						goto errorExit;
						}
					/* ... and then put the character into the matrix. Note that matrix
					   is defined as an integer, but we may have floating precision continuous
					   characters. To get around this, we multiply the value of the character
					   by 1000 before putting it into matrix. We will divide by 1000 later on
					   when/if we use the characters. */
					sscanf (tkn, "%lf", &charValue);
					charValue *= 1000.0;
					if (isNegative == YES)
						{
						charValue *= -1.0;
						isNegative = NO;
						}
					/*MrBayesPrint ("%d \n", (int)charValue);*/
					matrix[pos(taxonCount-1,taxaInfo[taxonCount-1].charCount++,numChar)] = (int)charValue;
					}
				}
			}
		else
			{
			/* Otherwise, we are dealing with a run-of-the-mill character, and we
		       cannot expect the entire token to contain only a single character. We
		       must, therefore, go through the token character-by-character. */
			i = 0;
			while (tkn[i] != '\0')
				{
				/*MrBayesPrint ("%c", tkn[i]);*/
				if (tkn[i] == matchId)
					{
					if (taxonCount == 1)
						{
						MrBayesPrint ("%s   Matching characters cannot be in first taxon\n", spacer);
						goto errorExit;
						}
					charCode = matrix[pos(0,taxaInfo[taxonCount-1].charCount,numChar)];
					matrix[pos(taxonCount-1,taxaInfo[taxonCount-1].charCount++,numChar)] = charCode;
					}
				else
					{
					if ((tkn[i] == ')' && isInAmbig == YES) || (tkn[i] == '}' && isInPoly == YES))
						{
						isInAmbig = isInPoly = NO;
						charCode = theAmbigChar;
						j = CharacterNumber (charCode, charInfo[taxaInfo[taxonCount-1].charCount].charType);
						if (j > charInfo[taxaInfo[taxonCount-1].charCount].numStates)
							charInfo[taxaInfo[taxonCount-1].charCount].numStates = j;
						matrix[pos(taxonCount-1,taxaInfo[taxonCount-1].charCount++,numChar)] = charCode;
						theAmbigChar = 0;
						}
					else if ((tkn[i] == '(' && isInAmbig == YES) || (tkn[i] == '{' && isInPoly == YES))
						{
						if (isInAmbig == YES)
							MrBayesPrint ("%s   Found an inappropriate \"(\"\n", spacer);
						else
							MrBayesPrint ("%s   Found an inappropriate \"{\"\n", spacer);
						goto errorExit;
						}
					else if (isInAmbig == YES || isInPoly == YES)
						{
						if (tkn[i] == ',')
							expecting |= Expecting (COMMA);
						else 
							{
							if (CharacterCode(tkn[i], &charCode, charInfo[taxaInfo[taxonCount-1].charCount].charType) == ERROR)
								goto errorExit;
							if (charCode == MISSING || charCode == GAP)
								goto errorExit;
							theAmbigChar |= charCode;
							expecting ^= Expecting (COMMA);
							}
						}
					else if (tkn[i] == '{' && isInPoly == NO && isInAmbig == NO)
						{
						isInPoly = YES;
						matrixHasPoly = YES;
						theAmbigChar = 0;
						}	
					else if (tkn[i] == '(' && isInPoly == NO && isInAmbig == NO)
						{
						isInAmbig = YES;
						theAmbigChar = 0;
						}
					else if (tkn[i] == '(' && isInPoly == NO && isInAmbig == NO)
						{
						isInAmbig = YES;
						theAmbigChar = 0;
						}
					else
						{
						if (CharacterCode(tkn[i], &charCode, charInfo[taxaInfo[taxonCount-1].charCount].charType) == ERROR)
							{
							MrBayesPrint ("%s   Error while reading character position %d (charCode %d)\n", spacer, taxaInfo[taxonCount-1].charCount, charCode);
							goto errorExit;
							}
						if (charCode != MISSING && charCode != GAP)
							{
							j = CharacterNumber (charCode, charInfo[taxaInfo[taxonCount-1].charCount].charType);
							if (j > charInfo[taxaInfo[taxonCount-1].charCount].numStates)
								charInfo[taxaInfo[taxonCount-1].charCount].numStates = j;
							}
						matrix[pos(taxonCount-1,taxaInfo[taxonCount-1].charCount++,numChar)] = charCode;
						}
					}
				i++; 
				}
			}
		}

	return (NO_ERROR);
	MrBayesPrint ("%s", parmName); /* just because I am tired of seeing the unused parameter error msg */
	errorExit:
		FreeMatrix();
		return (ERROR);

}





int DoNexusParm (char *parmName, char *tkn)

{

	if (!strcmp(parmName, "NEXUS"))
		{
		MrBayesPrint ("%s   Expecting NEXUS formatted file\n", spacer);
		expecting = Expecting(COMMAND);
		}
	else
		{
		MrBayesPrint ("%s   Found %s\n", spacer, tkn);
		return (ERROR);
		}
	
	return (NO_ERROR);
	
}





int DoOutgroup (void)

{

	char		tempName[100];

	if (GetNameFromString (taxaNames, tempName, outGroupNum + 1) == ERROR)
		{
		MrBayesPrint ("%s   Error getting taxon names \n", spacer);
		return (ERROR);
		}
	MrBayesPrint ("%s   Setting outgroup to taxon \"%s\"\n", spacer, tempName);
	
	if (isUserTreeDefined == YES && userTree->isRooted == NO)
		{
		MoveCalculationRoot (userTree, outGroupNum);
		}

	return (NO_ERROR);
	
}





int DoOutgroupParm (char *parmName, char *tkn)


{

	int		howMany, tempInt;

	if (expecting == Expecting(ALPHA))
		{
		if (CheckString (tkn, taxaNames, &howMany) == ERROR)
			{
			MrBayesPrint ("%s   Could not find taxon %s in list of taxa\n", spacer, tkn);
			return (ERROR);
			}
		outGroupNum = howMany - 1;
		
		expecting = Expecting(SEMICOLON);
		}
	else if (expecting == Expecting(NUMBER))
		{
		if (CheckString (tkn, taxaNames, &howMany) == ERROR)
			{
			/* OK, as we expect, the taxon is not a digit. So, now we assume that
			   the user is assigning the outgroup by its number */
			sscanf (tkn, "%d", &tempInt);
			if (tempInt < 1 || tempInt > numTaxa)
				{
				MrBayesPrint ("%s   Taxon number %d is out of range\n", spacer, tempInt);
				return (ERROR);
				}
			outGroupNum = tempInt - 1;
			}
		else
			{
			outGroupNum = howMany - 1;
			}
		
		expecting = Expecting(SEMICOLON);
		}
	else
		return (ERROR);

	return (NO_ERROR);
	MrBayesPrint ("%s", parmName); /* just because I am tired of seeing the unused parameter error msg */
	
}





int DoPairs (void)

{

	MrBayesPrint ("\n");
	MrBayesPrint ("%s   Successfully defined character pairings\n", spacer);

	defPairs = YES;
	foundFirst = NO;
	
	return (NO_ERROR);
	
}





int DoPairsParm (char *parmName, char *tkn)

{

	int		i, tempInt;
		
	if (defMatrix == NO)
		{
		MrBayesPrint ("%s   A matrix must be specified before you can define pairs of characters\n", spacer);
		return (ERROR);
		}
	
	if (defPairs == YES)
		{
		MrBayesPrint ("%s   Character pairs have been previously defined \n", spacer);
		MrBayesPrint ("%s   Now overwriting old pairings\n", spacer);
		for (i=0; i<numChar; i++)
			charInfo[i].pairsId = 0;
		defPairs = NO;
		}
		
	if (foundFirst == NO)
		{
		/* this is the first time in */
		pairId = 1;
		firstPair = YES;
		foundFirst = YES;
		MrBayesPrint ("%s   Defining character pairings:\n\n", spacer);
		MrBayesPrint ("%s      Pair --  First Second \n", spacer);
		}

	if (expecting == Expecting(NUMBER))
		{
		sscanf (tkn, "%d", &tempInt);
		if (tempInt <= 0 || tempInt > numChar)
			{
			MrBayesPrint ("\n");
			MrBayesPrint ("%s   Character number %d is out of range (should be between %d and %d)\n", spacer, tempInt, 1, numChar);
			for (i=0; i<numChar; i++)
				charInfo[i].pairsId = 0;
			return (ERROR);
			}
		tempInt--;
		
		if (charInfo[tempInt].pairsId != 0)
			{
			MrBayesPrint ("\n");
			MrBayesPrint ("%s   Character number %d has already been included in a pairing\n", spacer, tempInt+1);
			for (i=0; i<numChar; i++)
				charInfo[i].pairsId = 0;
			return (ERROR);
			}
		if (charInfo[tempInt].charType != DNA && charInfo[tempInt].charType != RNA)
			{
			MrBayesPrint ("\n");
			MrBayesPrint ("%s  Pairings may only include nucleotide data\n", spacer);
			if (charInfo[tempInt].charType == PROTEIN)
				MrBayesPrint ("%s  Character %d is an amino acid character\n", spacer, tempInt+1);
			else if (charInfo[tempInt].charType == RESTRICTION)
				MrBayesPrint ("%s  Character %d is a restriction site character\n", spacer, tempInt+1);
			else if (charInfo[tempInt].charType == STANDARD)
				MrBayesPrint ("%s  Character %d is a \"standard\" character\n", spacer, tempInt+1);
			else if (charInfo[tempInt].charType == CONTINUOUS)
				MrBayesPrint ("%s  Character %d is a continuously varying character\n", spacer, tempInt+1);
			for (i=0; i<numChar; i++)
				charInfo[i].pairsId = 0;
			return (ERROR);
			}
			
		charInfo[tempInt].pairsId = pairId;
		
		if (firstPair == YES)
			{
			MrBayesPrint ("%s      %4d --  %5d  ", spacer, pairId, tempInt+1);
			expecting  = Expecting(COLON);
			firstPair = NO;
			}
		else
			{
			MrBayesPrint ("%5d\n", tempInt+1);
			expecting  = (Expecting(COMMA) | Expecting(SEMICOLON));
			firstPair = YES;
			}
		}
	else if (expecting == Expecting(COMMA))
		{
		pairId++;
		expecting = Expecting(NUMBER);
		}
	else if (expecting == Expecting(COLON))
		{
		expecting = Expecting(NUMBER);
		}
	else
		{
		for (i=0; i<numChar; i++)
			charInfo[i].pairsId = 0;
		return (ERROR);
		}

	return (NO_ERROR);
	MrBayesPrint ("%s", parmName); /* just because I am tired of seeing the unused parameter error msg */

}





int DoPartition (void)

{

	int		i, j, howMany, numDivs, partTypes[MAX_NUM_DIVS], checkEachPart[MAX_NUM_DIVS];
	
	/* add set to tempSet */
	if (fromI >= 0)
		if (AddToSet (fromI, toJ, everyK, whichPartition+1) == ERROR)
			return (ERROR);

	/* add temporary information to partition */
	for (i=0; i<numChar; i++)
		{
		if (tempSet[i] == 0)
			{
			MrBayesPrint ("%s   Character %d not included in partition\n", spacer, i+1);
			return (ERROR);
			}
		/*MrBayesPrint ("%4d %4d \n", i, tempSet[i]);*/
		}
		
	/* check how many partitions were found against how many were expected */
	if (whichPartition != numPartitions - 1)
		{
		MrBayesPrint ("%s   Did not find correct number of partitions (expecting %d, found %d)\n", spacer, numPartitions, whichPartition + 1);
		return (ERROR);
		}
		
	/* add name to list of valid partitions */
	if (AddToString (tempSetName, partitionNames, &howMany) == ERROR)
		{
		MrBayesPrint ("%s   Problem adding charset %s to list\n", spacer, tempSetName);
		return (ERROR);
		}
		
	/* increment number of defined partitions */
	numDefinedPartitions++;

	/* enter partition information */
	for (i=0; i<numChar; i++)
		charInfo[i].partitionId[numDefinedPartitions-1] = tempSet[i];
		
	/* make certain that the partition labels go from 1 - numPartitions, inclusive */
	for (i=0; i<numPartitions; i++)
		checkEachPart[i] = NO;
	for (i=0; i<numChar; i++)
		checkEachPart[ charInfo[i].partitionId[numDefinedPartitions-1] - 1 ] = YES;
	for (i=0; i<numPartitions; i++)
		{
		if (checkEachPart[i] == NO)
			{
			MrBayesPrint ("%s   Could not find a division for character division %d\n", spacer, i+1);
			return (ERROR);
			}
		}
		
	/* check if partition overruns data types */
	numDivs = GetNumPartDivisions (numDefinedPartitions);
	for (i=0; i<MAX_NUM_DIVS; i++)
		partTypes[i] = -1;
	for (i=0; i<numChar; i++)
		{
		if (partTypes[ charInfo[i].partitionId[numDefinedPartitions-1]-1 ] == -1)
			partTypes[ charInfo[i].partitionId[numDefinedPartitions-1]-1 ] = charInfo[i].charType;
		else
			{
			if (partTypes[ charInfo[i].partitionId[numDefinedPartitions-1]-1 ] != charInfo[i].charType)
				{
				MrBayesPrint ("%s   There are two different data types for partition division %d\n", spacer, charInfo[i].partitionId[numDefinedPartitions-1]);
				if (RemoveLastFromString (partitionNames) == ERROR)
					return (ERROR);
				numDefinedPartitions--;
				for (j=0; j<numChar; j++)
					charInfo[j].partitionId[numDefinedPartitions] = 0;
				return (ERROR);
				}
			}
		}
	
	return (NO_ERROR);
	
}





int DoPartitionParm (char *parmName, char *tkn)

{

	int		i, howMany, tempInt;
	
	if (defMatrix == NO)
		{
		MrBayesPrint ("%s   A matrix must be specified before partitions can be defined\n", spacer);
		return (ERROR);
		}

	if (expecting == Expecting(PARAMETER))
		{
		/* set Partition ( ) ******************************************************************/
		if (!strcmp(parmName, "Xxxxxxxxxx"))
			{
			/* check size of partition name */
			if (strlen(tkn) > 99)
				{
				MrBayesPrint ("%s   Partition name is too long\n", spacer);
				return (ERROR);
				}
				
			/* check to see if the name has already been used as a partition */
			if (numDefinedPartitions > 1)
				{
				if (CheckString (tkn, partitionNames, &howMany) == ERROR)
					{
					/* if the partition name has not been used, then we should have an ERROR returned */
					/* we _want_ to be here */

					}
				else
					{
					MrBayesPrint ("%s   Partition name '%s' has been used previously\n", spacer, tkn);
					return (ERROR);
					}
				}
			else if (numDefinedPartitions > MAX_NUM_PARTITIONS)
				{
				MrBayesPrint ("%s   You cannot define more than %d partitions\n", spacer, MAX_NUM_PARTITIONS);
				return (ERROR);
				}
				
			/* add the name temporarily to tempSetName */
			strcpy (tempSetName, tkn);
			
			/* clear tempSet */
			for (i=0; i<numChar; i++)
				tempSet[i] = 0;
			
			fromI = toJ = everyK = -1;
			foundDash = foundSlash = NO;
			whichPartition = 0;
			foundFirst = NO;
			numPartitions = 0;
			MrBayesPrint ("%s   Defining partition called %s\n", spacer, tkn);
			expecting = Expecting(EQUALSIGN);
			}
		else
			return (ERROR);
		}
	else if (expecting == Expecting(EQUALSIGN))
		{
		expecting = Expecting(NUMBER);
		}
	else if (expecting == Expecting(ALPHA))
		{
		/* We are defining a partition in terms of a character set (called tkn, here). We should be able
		   to find tkn in the list of character set names. If we cannot, then we have a problem and
		   return an error. */
		if (numCharSets < 1)
			{
			MrBayesPrint ("%s   Could not find a character set called %s\n", spacer, tkn);
			return (ERROR);
			}
		if (CheckString (tkn, charSetNames, &howMany) == ERROR)
			{
			MrBayesPrint ("%s   Could not find a character set called %s\n", spacer, tkn);
			return (ERROR);
			}
		/* add characters from charset tkn to new tempSet */
		for (i=0; i<numChar; i++)
			{
			if (charInfo[i].charSet[howMany-1] == 1)
				tempSet[i] = whichPartition + 1;
			}
		fromI = toJ = everyK = -1;

		expecting  = Expecting(ALPHA);
		expecting |= Expecting(NUMBER);
		expecting |= Expecting(SEMICOLON);
		expecting |= Expecting(COMMA);
		}
	else if (expecting == Expecting(NUMBER))
		{
		if (foundFirst == NO)
			{
			sscanf (tkn, "%d", &tempInt);
			numPartitions = tempInt;
			expecting  = Expecting(COLON);
			foundFirst = YES;
			}
		else
			{
			if (strlen(tkn) == 1 && tkn[0] == '.')
				tempInt = numChar;
			else
				sscanf (tkn, "%d", &tempInt);
			if (tempInt <= 0 || tempInt > numChar)
				{
				MrBayesPrint ("%s   Character number %d is out of range (should be between %d and %d)\n", spacer, tempInt, 1, numChar);
				return (ERROR);
				}
			tempInt--;
			if (foundDash == YES)
				{
				if (fromI >= 0)
					toJ = tempInt;
				else
					{
					MrBayesPrint ("%s   Improperly formatted partition\n", spacer);
					return (ERROR);
					}
				foundDash = NO;
				}
			else if (foundSlash == YES)
				{
				tempInt++;
				if (tempInt <= 1)
					{
					MrBayesPrint ("%s   Improperly formatted charset\n", spacer);
					return (ERROR);
					}
				if (fromI >= 0 && toJ >= 0 && fromI < toJ)
					everyK = tempInt;
				else
					{
					MrBayesPrint ("%s   Improperly formatted charset\n", spacer);
					return (ERROR);
					}
				foundSlash = NO;
				}
			else
				{
				if (fromI >= 0 && toJ < 0)
					{
					if (AddToSet (fromI, toJ, everyK, whichPartition+1) == ERROR)
						return (ERROR);
					fromI = tempInt;
					}
				else if (fromI < 0 && toJ < 0)
					{
					fromI = tempInt;
					}
				else if (fromI >= 0 && toJ >= 0 && everyK < 0)
					{
					if (AddToSet (fromI, toJ, everyK, whichPartition+1) == ERROR)
						return (ERROR);
					fromI = tempInt;
					toJ = everyK = -1;
					}
				else if (fromI >= 0 && toJ >= 0 && everyK >= 0)
					{
					if (AddToSet (fromI, toJ, everyK, whichPartition+1) == ERROR)
						return (ERROR);
					fromI = tempInt;
					toJ = everyK = -1;
					}
				else
					{
					MrBayesPrint ("%s   Improperly formatted charset\n", spacer);
						{
						return (ERROR);
						}
					}
				}
			
			expecting  = Expecting(ALPHA);
			expecting |= Expecting(NUMBER);
			expecting |= Expecting(SEMICOLON);
			expecting |= Expecting(DASH);
			expecting |= Expecting(BACKSLASH);
			expecting |= Expecting(COMMA);
			}
		}
	else if (expecting == Expecting(COMMA))
		{
		/* add set to tempSet */
		if (fromI >= 0)
			if (AddToSet (fromI, toJ, everyK, whichPartition+1) == ERROR)
				return (ERROR);

		fromI = toJ = everyK = -1;
		foundDash = foundSlash = NO;
		whichPartition++;
		if (whichPartition > numPartitions)
			{
			MrBayesPrint ("%s   Too many partitions of the data (expecting %d)\n", spacer, numPartitions);
			return (ERROR);
			}
		expecting  = Expecting(NUMBER);
		expecting |= Expecting(ALPHA);
		}
	else if (expecting == Expecting(COLON))
		{
		expecting  = Expecting(NUMBER);
		expecting |= Expecting(ALPHA);
		}
	else if (expecting == Expecting(DASH))
		{
		foundDash = YES;
		expecting = Expecting(NUMBER);
		}
	else if (expecting == Expecting(BACKSLASH))
		{
		foundSlash = YES;
		expecting = Expecting(NUMBER);
		}
	else
		return (ERROR);

	return (NO_ERROR);

}





int DoProps (void)

{

	int			i, j, j1, j2, k, index, whichParamToChange, invalidProp;
	MrBFlt		tempD;
	char		tempString[100];
	MoveType	*mt = NULL;
	
	for (j2=0;j2<10;j2++)
		{
		/* display options */
		MrBayesPrint ("%s   Available proposal mechanisms:\n\n", spacer);
		for (i=index=0; i<numMoveTypes; i++)
			{
			mt = &moveTypes[i];

			if (userLevel < mt->level)
				continue;
			
			MrBayesPrint ("%s    %2d -- Change %s\n", spacer, ++index, mt->name);
			MrBayesPrint ("%s             proposal rate: %1.3lf\n", spacer, mt->relProposalProb);
			for (j=0; j<mt->numTuningParams; j++)
				MrBayesPrint ("%s             %s: %1.3lf\n", spacer, mt->nameTuning[j], mt->tuningParam[j]);
			}
		MrBayesPrint ("\n");

		/* pick an option */
		for (j1=0; j1<10; j1++) 
			{
			MrBayesPrint ("%s   Select a parameter to change (1 - %d; 0 to exit; %d to zero all proposal rates): ", spacer, index, index + 1);
			fgets (tempString, 100, stdin);
			sscanf (tempString, "%d", &whichParamToChange);
			if (whichParamToChange >= 0 && whichParamToChange <= index + 1)
				break;
			}

		if (whichParamToChange == 0)
			{
			/* exit */
			MrBayesPrint ("%s   Exit props\n", spacer);
			return (NO_ERROR);
			}
		else if (whichParamToChange == index + 1)
			{
			/* zero all proposal rates */
			MrBayesPrint ("%s   Setting all proposal rates to zero (this is very dangerous)\n", spacer);
			for (i=0; i<numMoveTypes; i++)
				{
				mt = &moveTypes[i];
				if (userLevel >= mt->level)
					mt->relProposalProb = 0.0;
				}
			}
		else if (whichParamToChange > 0 && whichParamToChange <= index)
			{
			/* find the proposal to change */
			for (i=index=0; i<numMoveTypes; i++)
				{
				mt = &moveTypes[i];
				if (userLevel < mt->level)
					continue;
				if (++index == whichParamToChange)
					break;
				}

			/* change the parameters for the proposal mechanism */
		    MrBayesPrint ("%s   Proposal %d: Change %s\n", spacer, whichParamToChange, mt->name);
			for (j1=0; j1<10; j1++)
				{
				invalidProp = NO;
				MrBayesPrint ("%s      New proposal rate (<return> to keep old = %1.3lf):  ", spacer, mt->relProposalProb);
				fgets (tempString, 100, stdin);
				k = sscanf (tempString, "%lf", &tempD);
				if (k == 1 && tempD > 0.0 && tempD < 1000.0)
					mt->relProposalProb = tempD;
				else if (k != EOF)
					invalidProp = YES;
				if (invalidProp == YES)
					MrBayesPrint ("%s      Incorrect value...try again\n", spacer);
				else
					break;
				}
			for (j=0; j<mt->numTuningParams; j++)
				{
				for (j1=0; j1<10; j1++)
					{
					invalidProp = NO;
					MrBayesPrint ("%s      New %s (<return> to keep old = %1.3lf): ", spacer, mt->nameTuning[j], mt->tuningParam[j]);
					fgets (tempString, 100, stdin);
					k = sscanf (tempString, "%lf", &tempD);
					if (k == 1 && tempD >= mt->minimum[j] && tempD <= mt->maximum[j])
						mt->tuningParam[j] = tempD;
					else if (k != EOF)
						invalidProp = YES;
					if (invalidProp == YES)
						MrBayesPrint ("%s      Incorrect value...try again\n", spacer);
					else
						break;
					}
				}
			}
		else
			MrBayesPrint ("%s      Too many invalid entries\n", spacer);
		}
		
	return (NO_ERROR);
	
}





int DoRestore (void)

{

	int			i, alreadyDone;

	MrBayesPrint ("%s   Restore taxa\n", spacer);

	/* add set to tempSet */
	if (fromI >= 0 && toJ < 0)
		{
		if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
			return (ERROR);
		}
	else if (fromI >= 0 && toJ >= 0 && everyK < 0)
		{
		if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
			return (ERROR);
		}
	else if (fromI >= 0 && toJ >= 0 && everyK >= 0)
		{
		if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
			return (ERROR);
		}
		
	/* merge tempSet with excludedTaxa */
	alreadyDone = NO;
	for (i=0; i<numTaxa; i++)
		{
		if (tempSet[i] == 1)
			{
			if (taxaInfo[i].isDeleted == NO && alreadyDone == NO)
				{
				MrBayesPrint ("%s   Some taxa already included\n", spacer);
				alreadyDone = YES;
				}
			taxaInfo[i].isDeleted = NO;
			}
		}

	/* show tempSet (for debugging) */
#	if 0
	for (i=0; i<numTaxa; i++)
		MrBayesPrint ("%4d  %4d\n", i+1, tempSet[i]);
#	endif
		
	return (NO_ERROR);
	
}





int DoRestoreParm (char *parmName, char *tkn)

{

	int		i, howMany, tempInt;
		
	if (defMatrix == NO)
		{
		MrBayesPrint ("%s   A matrix must be specified before you can restore taxa\n", spacer);
		return (ERROR);
		}
		
	if (foundFirst == NO)
		{
		/* this is the first time in */
		fromI = toJ = everyK = -1;
		foundDash = NO;
		for (i=0; i<numTaxa; i++) /* clear tempSet */
			tempSet[i] = 0;
		foundFirst = YES;
		}

	if (expecting == Expecting(ALPHA))
		{
		if (IsSame ("All", tkn) == SAME || IsSame ("All", tkn) == CONSISTENT_WITH)
			{
			for (i=0; i<numTaxa; i++)
				tempSet[i] = 1;
			}
		else
			{
			if (CheckString (tkn, taxaNames, &howMany) == ERROR)
				{
				/* we are using a pre-defined taxa set */
				if (numTaxaSets < 1)
					{
					MrBayesPrint ("%s   Could not find a taxset called %s\n", spacer, tkn);
					return (ERROR);
					}
				if (CheckString (tkn, taxaSetNames, &howMany) == ERROR)
					{
					MrBayesPrint ("%s   Could not find a taxset called %s\n", spacer, tkn);
					return (ERROR);
					}
				/* add taxa from taxaSet tkn to new tempSet */
				for (i=0; i<numTaxa; i++)
					{
					if (taxaInfo[i].taxaSet[howMany-1] == 1)
						tempSet[i] = 1;
					}
				}
			else
				{
				/* we found the taxon name */
				if (fromI >= 0 && toJ < 0)
					{
					if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
						return (ERROR);
					}
				else if (fromI >= 0 && toJ >= 0)
					{
					if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
						return (ERROR);
					}
				tempSet[howMany-1] = 1;
				}
			fromI = toJ = everyK = -1;
			}

		expecting  = Expecting(ALPHA);
		expecting |= Expecting(NUMBER);
		expecting |= Expecting(SEMICOLON);
		}
	else if (expecting == Expecting(NUMBER))
		{
		if (strlen(tkn) == 1 && !strcmp(tkn, "."))
			{
			tempInt = numTaxa;
			}
		else
			{
			sscanf (tkn, "%d", &tempInt);
			if (tempInt <= 0 || tempInt > numTaxa)
				{
				MrBayesPrint ("%s   Taxon number %d is out of range (should be between %d and %d)\n", spacer, tempInt, 1, numTaxa);
				return (ERROR);
				}
			}
		tempInt--;
		if (foundDash == YES)
			{
			if (fromI >= 0)
				toJ = tempInt;
			else
				{
				MrBayesPrint ("%s   Improperly formatted restore set\n", spacer);
				return (ERROR);
				}
			foundDash = NO;
			}
		else
			{
			if (fromI >= 0 && toJ < 0)
				{
				if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
					return (ERROR);
				fromI = tempInt;
				}
			else if (fromI < 0 && toJ < 0)
				{
				fromI = tempInt;
				}
			else if (fromI >= 0 && toJ >= 0 && everyK < 0)
				{
				if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
					return (ERROR);
				fromI = tempInt;
				toJ = everyK = -1;
				}
			else if (fromI >= 0 && toJ >= 0 && everyK >= 0)
				{
				if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
					return (ERROR);
				fromI = tempInt;
				toJ = everyK = -1;
				}
			else
				{
				MrBayesPrint ("%s   Improperly formatted restore set\n", spacer);
					{
					return (ERROR);
					}
				}
			}
		expecting  = Expecting(ALPHA);
		expecting |= Expecting(NUMBER);
		expecting |= Expecting(SEMICOLON);
		expecting |= Expecting(DASH);
		}
	else if (expecting == Expecting(DASH))
		{
		foundDash = YES;
		expecting = Expecting(NUMBER);
		}
	else
		return (ERROR);

	return (NO_ERROR);
	MrBayesPrint ("%s", parmName); /* just because I am tired of seeing the unused parameter error msg */

}





int DoRoot (void)

{
	char	tempName[100];
	
	if (isUserTreeDefined == NO)
		{
		MrBayesPrint ("   No tree in memory\n", spacer);
		return (ERROR);
		}
	
	if (userTree->isRooted == NO)
		{
		if (RootTree (userTree) == ERROR)
			return (ERROR);
		GetNameFromString (taxaNames, tempName, outGroupNum);
		MrBayesPrint ("%s   Rooting user tree midway between the outgroup (%s) and the ingroup\n", spacer, tempName);
		}
	else
		MrBayesPrint ("%s   Tree is already rooted\n", spacer);

	return (NO_ERROR);
	
}





int DoSet (void)

{

	return (NO_ERROR);
	
}





int DoSetParm (char *parmName, char *tkn)

{

	int			howMany, nDivs;
	char		tempStr[100];

	if (expecting == Expecting(PARAMETER))
		{
		expecting = Expecting(EQUALSIGN);
		}
	else
		{
		/* set Autoclose (autoClose) **********************************************************/
		if (!strcmp(parmName, "Autoclose"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting(ALPHA);
			else if (expecting == Expecting(ALPHA))
				{
				if (IsArgValid(tkn, tempStr) == NO_ERROR)
					{
					if (!strcmp(tempStr, "Yes"))
						autoClose = YES;
					else
						autoClose = NO;
					}
				else
					{
					MrBayesPrint ("%s   Invalid argument for autoclose\n", spacer);
					return (ERROR);
					}
				if (autoClose == YES)
					MrBayesPrint ("%s   Setting autoclose to yes\n", spacer);
				else
					MrBayesPrint ("%s   Setting autoclose to no\n", spacer);
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else
				return (ERROR);
			}
		/* set Nowarnings (noWarn) **********************************************************/
		else if (!strcmp(parmName, "Nowarnings"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting(ALPHA);
			else if (expecting == Expecting(ALPHA))
				{
				if (IsArgValid(tkn, tempStr) == NO_ERROR)
					{
					if (!strcmp(tempStr, "Yes"))
						noWarn = YES;
					else
						noWarn = NO;
					}
				else
					{
					MrBayesPrint ("%s   Invalid argument for nowarnings\n", spacer);
					return (ERROR);
					}
				if (noWarn == YES)
					MrBayesPrint ("%s   Setting nowarnings to yes\n", spacer);
				else
					MrBayesPrint ("%s   Setting nowarnings to no\n", spacer);
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else
				return (ERROR);
			}
		/* set Quitonerror (quitOnError) **************************************************/
		else if (!strcmp(parmName, "Mode"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting(ALPHA);
			else if (expecting == Expecting(ALPHA))
				{
				if (IsArgValid(tkn, tempStr) == NO_ERROR)
					{
					if (!strcmp(tempStr, "Yes"))
						quitOnError = YES;
					else
						quitOnError = NO;
					}
				else
					{
					MrBayesPrint ("%s   Invalid argument for quitonerror\n", spacer);
					return (ERROR);
					}
				if (quitOnError == YES)
					MrBayesPrint ("%s   Setting quitonerror to yes\n", spacer);
				else
					MrBayesPrint ("%s   Setting quitonerror to no\n", spacer);
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else
				return (ERROR);
			}
		/* set AutoOverwrite (autoOverwrite) **************************************************/
		else if (!strcmp(parmName, "Autooverwrite"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting(ALPHA);
			else if (expecting == Expecting(ALPHA))
				{
				if (IsArgValid(tkn, tempStr) == NO_ERROR)
					{
					if (!strcmp(tempStr, "Yes"))
						{
						autoOverwrite = YES;
						MrBayesPrint ("%s   Setting autooverwrite to yes\n", spacer);
						}
					else
						{
						autoOverwrite = NO;
						MrBayesPrint ("%s   Setting autooverwrite to no\n", spacer);
						}
					}
				else
					{
					MrBayesPrint ("%s   Invalid argument for autooverwrite\n", spacer);
					return (ERROR);
					}					
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else
				return (ERROR);
			}			
		/* set Partition (partitionNum) *******************************************************/
		else if (!strcmp(parmName, "Partition"))
			{
			if (defMatrix == NO)
				{
				MrBayesPrint ("%s   A character matrix must be defined first\n", spacer);
				return (ERROR);
				}
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting(ALPHA) | Expecting(NUMBER);
			else if (expecting == Expecting(ALPHA))
				{
				/* first check to see if name is there */
				if (CheckString (tkn, partitionNames, &howMany) == ERROR)
					{
					MrBayesPrint ("%s   Could not find \"%s\" as a defined partition\n", spacer, tkn);
					return (ERROR);
					}
				SetModelDefaults ();
				partitionNum = howMany;
				nDivs = SetPartitionInfo (howMany-1);
				numCurrentDivisions = nDivs;
				if (nDivs == 1)
					MrBayesPrint ("%s   Setting %s as the partition (does not divide up characters).\n", spacer, tkn); 
				else
					MrBayesPrint ("%s   Setting %s as the partition, dividing characters into %d parts.\n", spacer, tkn, nDivs); 
				MrBayesPrint ("%s   Resetting model values to defaults (NB! Any existing model settings will be deleted!)\n", spacer);
				MrBayesPrint ("%s   Reinitializing link table (linking all parameters)\n", spacer);
				if (InitializeLinks () == ERROR)
					{
					MrBayesPrint ("%s   Problem initializing link table\n", spacer);
					return (ERROR);
					}
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else if (expecting == Expecting(NUMBER))
				{
				sscanf (tkn, "%d", &howMany);
				if (howMany > numDefinedPartitions) 
					{
					MrBayesPrint ("%s   Partition number %d is not a valid parition. Only %d partitions\n", spacer, howMany, numDefinedPartitions);
					MrBayesPrint ("%s   have been defined.\n", spacer);
					return (ERROR);
					}
				if (howMany < 1)
					{
					MrBayesPrint ("%s   Partition number %d is not a valid parition. Must be between 1 and %d.\n", spacer, howMany, numDefinedPartitions);
					return (ERROR);
					}
				if (GetNameFromString (partitionNames, tempStr, howMany) == ERROR)
					{
					MrBayesPrint ("%s   Could not find name for partition number %d.\n", spacer, howMany);
					return (ERROR);
					}
				SetModelDefaults ();
				partitionNum = howMany;
				nDivs = SetPartitionInfo (howMany-1);
				numCurrentDivisions = nDivs;
				if (nDivs == 1)
					MrBayesPrint ("%s   Setting %s as the partition (does not divide up characters).\n", spacer, tempStr); 
				else
					MrBayesPrint ("%s   Setting %s as the partition, dividing characters into %d parts.\n", spacer, tempStr, nDivs); 
				MrBayesPrint ("%s   Resetting model values to defaults\n", spacer);
				MrBayesPrint ("%s   Reinitializing link table (linking all parameters)\n", spacer);
				if (InitializeLinks () == ERROR)
					{
					MrBayesPrint ("%s   Problem initializing link table\n", spacer);
					return (ERROR);
					}
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else
				return (ERROR);
			}
			
		else
			return (ERROR);
			
			
		}

	return (NO_ERROR);
	
}





int DoShowMatrix (void)

{

	int			i, j, nameLen, start, finish, ct, longestName;
	char		tempStr[100], stride;
	
	if (defMatrix == NO)
		{
		MrBayesPrint ("%s   A character matrix must be defined first\n", spacer);
		return (ERROR);
		}
			
	longestName = 0;
	for (i=0; i<numTaxa; i++)
		{
		if (GetNameFromString (taxaNames, tempStr, i+1) == ERROR)
			{
			MrBayesPrint ("%s   Could not find taxon %d\n", spacer, i+1);
			return (ERROR);
			}
		nameLen = (int) strlen(tempStr);
		if (nameLen > longestName)
			longestName = nameLen;
		}
			
	stride = 50;
	start = finish = 0;
	do
		{
		finish += stride;
		if (finish > numChar)
			finish = numChar;

		MrBayesPrint ("%s   ", spacer);
		for (j=0; j<longestName; j++)
			MrBayesPrint (" ");
		MrBayesPrint ("  ");
		MrBayesPrint ("%d\n", start+1);


		for (i=0; i<numTaxa; i++)
			{
			if (GetNameFromString (taxaNames, tempStr, i+1) == ERROR)
				{
				MrBayesPrint ("%s   Could not find taxon %d\n", spacer, i+1);
				return (ERROR);
				}
			nameLen = (int) strlen(tempStr);
			
			MrBayesPrint ("%s   ", spacer);
			if (nameLen >= longestName)
				{
				for (j=0; j<longestName; j++)
					MrBayesPrint ("%c", tempStr[j]);
				}
			else
				{
				MrBayesPrint ("%s", tempStr);
				for (j=0; j<longestName-nameLen; j++)
					MrBayesPrint (" ");
				}
			MrBayesPrint ("  ");

			for (j=start; j<finish; j++)
				{
				ct = charInfo[j].charType;
				if (ct == DNA || ct == RNA)
					MrBayesPrint ("%c", WhichNuc(matrix[pos(i,j,numChar)]));
				else if (ct == PROTEIN)
					MrBayesPrint ("%c", WhichAA(matrix[pos(i,j,numChar)]));
				else if (ct == STANDARD)
					MrBayesPrint ("%c", WhichStand(matrix[pos(i,j,numChar)]));
				else if (ct == RESTRICTION)
					MrBayesPrint ("%c", WhichRes(matrix[pos(i,j,numChar)]));
				else if (ct == CONTINUOUS)
					{
					if (WhichCont(matrix[pos(i,j,numChar)]) < 0.0)
						MrBayesPrint (" %2.2lf", WhichCont(matrix[pos(i,j,numChar)]));
					else
						MrBayesPrint ("  %2.2lf", WhichCont(matrix[pos(i,j,numChar)]));
					}
				else
					{
					MrBayesPrint ("%s   Unknown data type\n", spacer);
					return (ERROR);
					}
				
				}
			MrBayesPrint ("\n");
			}
		MrBayesPrint ("\n");
		start = finish;
		} while (finish != numChar);

	return (NO_ERROR);

}





int DoShowtree (void)

{


	if (isUserTreeDefined == YES)
		{
		if (userTree->isRooted == YES)
			MrBayesPrint ("\n   Tree is rooted:\n\n");
		else
			MrBayesPrint ("\n   Tree is unrooted:\n\n");
		if (ShowTree (userTree->root, userTree->isRooted, numTaxa) == ERROR)
			{
			return (ERROR);
			}
		else
			MrBayesPrint ("\n");
		}
	else
		{
		MrBayesPrint ("%s   A user tree has not been defined\n", spacer);
		}

	return (NO_ERROR);
	
}





int DoTaxaset (void)

{

	int			i, howMany;

	/* first add name to taxaSetName */
	if (AddToString (tempSetName, taxaSetNames, &howMany) == ERROR)
		{
		MrBayesPrint ("%s   Problem adding taxset %s to list\n", spacer, tempSetName);
		return (ERROR);
		}
	if (howMany != numTaxaSets + 1)
		{
		MrBayesPrint ("%s   Problem adding taxset %s to list\n", spacer, tempSetName, numTaxaSets);
		if (RemoveLastFromString (taxaSetNames) == ERROR)
			return (ERROR);
		return (ERROR);
		}

	/* add set to tempSet */
	if (fromI >= 0 && toJ < 0)
		{
		if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
			{
			if (RemoveLastFromString (taxaSetNames) == ERROR)
				return (ERROR);
			return (ERROR);
			}
		}
	else if (fromI >= 0 && toJ >= 0 && everyK < 0)
		{
		if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
			{
			if (RemoveLastFromString (taxaSetNames) == ERROR)
				return (ERROR);
			return (ERROR);
			}
		}
	else if (fromI >= 0 && toJ >= 0 && everyK >= 0)
		{
		if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
			{
			if (RemoveLastFromString (taxaSetNames) == ERROR)
				return (ERROR);
			return (ERROR);
			}
		}
		
	/* merge tempSet with taxaSet */
	for (i=0; i<numTaxa; i++)
		{
		if (tempSet[i] == 1)
			taxaInfo[i].taxaSet[numTaxaSets] = 1;
		}
	
	/* increment number of char sets */
	numTaxaSets++;
	
	/* show taxset (for debugging) */
#	if 0
	for (i=0; i<numTaxa; i++)
		MrBayesPrint ("%4d  %4d\n", i+1, taxaSet[i]);
#	endif

	return (NO_ERROR);
	
}





int DoTaxasetParm (char *parmName, char *tkn)

{

	int		i, howMany, tempInt;
	
	if (defMatrix == NO)
		{
		MrBayesPrint ("%s   A matrix must be specified before taxsets can be defined\n", spacer);
		return (ERROR);
		}

	if (expecting == Expecting(PARAMETER))
		{
		if (!strcmp(parmName, "Xxxxxxxxxx"))
			{
			/* check size of taxset name */
			if (strlen(tkn) > 99)
				{
				MrBayesPrint ("%s   Taxset name is too long\n", spacer);
				return (ERROR);
				}
				
			/* check to see if the name has already been used as a taxset */
			if (numTaxaSets > 0)
				{
				if (CheckString (tkn, taxaSetNames, &howMany) == ERROR)
					{
					/* if the taxset name has not been used, then we should have an ERROR returned */
					/* we _want_ to be here */

					}
				else
					{
					MrBayesPrint ("%s   Taxset name has been used previously\n", spacer);
					return (ERROR);
					}
				}
			else if (numTaxaSets > 30)
				{
				MrBayesPrint ("%s   You cannot define more than 30 taxsets\n", spacer);
				return (ERROR);
				}
				
			/* add the name to the taxa set */
			strcpy (tempSetName, tkn);
			
			/* clear tempSet */
			for (i=0; i<numTaxa; i++)
				tempSet[i] = 0;
			
			fromI = toJ = everyK = -1;
			foundDash = foundSlash = NO;
			MrBayesPrint ("%s   Defining taxset called %s\n", spacer, tkn);
			expecting = Expecting(EQUALSIGN);
			}
		else
			return (ERROR);
		}
	else if (expecting == Expecting(EQUALSIGN))
		{
		expecting  = Expecting(ALPHA);
		expecting |= Expecting(NUMBER);
		}
	else if (expecting == Expecting(ALPHA))
		{
		/* We are defining a taxon set in terms of another (called tkn, here) or we are referring to
		   the taxon name. We should be able to find tkn in the list of character set names or in the list
		   of taxon names. If we cannot, then we have a problem and return an error. */
		if (CheckString (tkn, taxaNames, &howMany) == ERROR)
			{
			if (numTaxaSets < 1)
				{
				MrBayesPrint ("%s   Could not find a taxset called %s\n", spacer, tkn);
				return (ERROR);
				}
			if (CheckString (tkn, taxaSetNames, &howMany) == ERROR)
				{
				MrBayesPrint ("%s   Could not find a taxset called %s\n", spacer, tkn);
				return (ERROR);
				}
			/* add taxa from taxset tkn to new tempSet */
			for (i=0; i<numTaxa; i++)
				{
				if (taxaInfo[i].taxaSet[howMany-1] == 1)
					tempSet[i] = 1;
				}
			}
		else
			{
			tempSet[howMany-1] = 1;
			}
		fromI = toJ = everyK = -1;

		expecting  = Expecting(ALPHA);
		expecting |= Expecting(NUMBER);
		expecting |= Expecting(SEMICOLON);
		}
	else if (expecting == Expecting(NUMBER))
		{
		if (strlen(tkn) == 1 && !strcmp(tkn, "."))
			{
			tempInt = numTaxa;
			}
		else
			{
			sscanf (tkn, "%d", &tempInt);
			if (tempInt <= 0 || tempInt > numTaxa)
				{
				MrBayesPrint ("%s   Taxon number %d is out of range (should be between %d and %d)\n", spacer, tempInt, 1, numTaxa);
				return (ERROR);
				}
			}
		tempInt--;
		if (foundDash == YES)
			{
			if (fromI >= 0)
				toJ = tempInt;
			else
				{
				MrBayesPrint ("%s   Improperly formatted taxset\n", spacer);
				return (ERROR);
				}
			foundDash = NO;
			}
		else if (foundSlash == YES)
			{
			tempInt++;
			if (tempInt <= 1)
				{
				MrBayesPrint ("%s   Improperly formatted taxset\n", spacer);
				return (ERROR);
				}
			if (fromI >= 0 && toJ >= 0 && fromI < toJ)
				everyK = tempInt;
			else
				{
				MrBayesPrint ("%s   Improperly formatted taxset\n", spacer);
				return (ERROR);
				}
			foundSlash = NO;
			}
		else
			{
			if (fromI >= 0 && toJ < 0)
				{
				if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
					return (ERROR);
				fromI = tempInt;
				}
			else if (fromI < 0 && toJ < 0)
				{
				fromI = tempInt;
				}
			else if (fromI >= 0 && toJ >= 0 && everyK < 0)
				{
				if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
					return (ERROR);
				fromI = tempInt;
				toJ = everyK = -1;
				}
			else if (fromI >= 0 && toJ >= 0 && everyK >= 0)
				{
				if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
					return (ERROR);
				fromI = tempInt;
				toJ = everyK = -1;
				}
			else
				{
				MrBayesPrint ("%s   Improperly formatted taxset\n", spacer);
					{
					return (ERROR);
					}
				}
			}
		
		expecting  = Expecting(ALPHA);
		expecting |= Expecting(NUMBER);
		expecting |= Expecting(SEMICOLON);
		expecting |= Expecting(DASH);
		expecting |= Expecting(BACKSLASH);
		}
	else if (expecting == Expecting(DASH))
		{
		foundDash = YES;
		expecting = Expecting(NUMBER);
		}
	else if (expecting == Expecting(BACKSLASH))
		{
		foundSlash = YES;
		expecting = Expecting(NUMBER);
		}
	else
		return (ERROR);

	return (NO_ERROR);

}





int DoTaxaStat (void)

{

	int			i, j, maxLen, nameLen, nIncludedTaxa;
	char		tempName[100];
	
	if (defMatrix == NO)
		{
		MrBayesPrint ("%s   A character matrix must be defined first\n", spacer);
		return (ERROR);
		}
		
	/* find maximum length of taxon name */
	maxLen = nIncludedTaxa = 0;
	for (i=0; i<numTaxa; i++)
		{
		if (GetNameFromString (taxaNames, tempName, i+1) == ERROR)
			{
			MrBayesPrint ("%s   Could not find taxon %d\n", spacer, i+1);
			return (ERROR);
			}
		if ((int)strlen(tempName) > maxLen)
			maxLen = (int) strlen(tempName);
		if (taxaInfo[i].isDeleted == NO)
			nIncludedTaxa++;
		}
			
	MrBayesPrint ("%s   Showing taxon status:\n\n", spacer);
	if (nIncludedTaxa == numTaxa)
		MrBayesPrint ("%s     Number of taxa        = %d (all of which are included)\n", spacer, numTaxa);
	else
		MrBayesPrint ("%s     Number of taxa        = %d (of which %d are included)\n", spacer, numTaxa, nIncludedTaxa);
	MrBayesPrint ("%s     Number of constraints = %d\n\n", spacer, numDefinedConstraints);
	
	if (numDefinedConstraints > 0)
		{
		for (j=0; j<numDefinedConstraints; j++)
			{
			if (GetNameFromString (constraintNames, tempName, j+1) == ERROR)
				{
				MrBayesPrint ("%s   Error getting constraint names \n", spacer);
				return (ERROR);
				}

			/* for now, ignore the probability */
			MrBayesPrint ("%s     %2d -- Trees with constraint \"%s\" are infinitely\n", spacer, j+1, tempName);
			MrBayesPrint ("%s           more probable than those without \n", spacer);

			/* 
			if (relConstraintProbs[j] > 0.0)
				{
				if (relConstraintProbs[j] > 900000.0)
					{
					MrBayesPrint ("%s     %2d -- Trees with constraint \"%s\" are infinitely\n", spacer, j+1, tempName);
					MrBayesPrint ("%s           more probable than those without \n", spacer);
					}
				else
					{
					MrBayesPrint ("%s     %2d -- Trees with constraint \"%s\" are exp(%1.2lf)\n", spacer, j+1, tempName, relConstraintProbs[j]);
					MrBayesPrint ("%s           times more probable than those without \n", spacer);
					}
				}
			else
				{
				if (relConstraintProbs[j] < -900000.0)
					{
					MrBayesPrint ("%s     %2d -- Trees with constraint \"%s\" are infinitely\n", spacer, j+1, tempName);
					MrBayesPrint ("%s           less probable than those without \n", spacer);
					}
				else
					{
					MrBayesPrint ("%s     %2d -- Trees with constraint \"%s\" are exp(%1.2lf)\n", spacer, j+1, tempName, -relConstraintProbs[j]);
					MrBayesPrint ("%s           times less probable than those without \n", spacer);
					}
				}
			*/
			}
		MrBayesPrint ("\n");
		for (j=0; j<maxLen; j++)
			MrBayesPrint (" ");
		MrBayesPrint ("                             Constraints\n");
		}
	MrBayesPrint ("%s     Taxon  ", spacer);
	for (j=0; j<maxLen; j++)
		MrBayesPrint (" ");
	MrBayesPrint ("   Inclusion");
	MrBayesPrint ("   ");
	for (j=0; j<numDefinedConstraints; j++)
		MrBayesPrint (" %2d", j+1);
	MrBayesPrint ("\n");
	MrBayesPrint ("%s   -------", spacer);
	for (j=0; j<maxLen; j++)
		MrBayesPrint ("-");
	MrBayesPrint ("--------------");
	
	if (numDefinedConstraints > 0)
		{
		MrBayesPrint ("----");
		for (j=0; j<numDefinedConstraints; j++)
			MrBayesPrint ("---");
		}
	MrBayesPrint ("\n");
	for (i=0; i<numTaxa; i++)
		{
		if (GetNameFromString (taxaNames, tempName, i+1) == ERROR)
			{
			MrBayesPrint ("%s   Could not find taxon %d\n", spacer, i+1);
			return (ERROR);
			}
		nameLen = (int) strlen(tempName);
		
		if (i == outGroupNum)
			MrBayesPrint ("%s ->%4d (%s) ", spacer, i+1, tempName);
		else
			MrBayesPrint ("%s   %4d (%s) ", spacer, i+1, tempName);
		for (j=0; j<(maxLen-nameLen); j++)
			MrBayesPrint (" ");
		MrBayesPrint (" -- ");
		
		if (taxaInfo[i].isDeleted == YES)
			MrBayesPrint ("Deleted ");
		else
			MrBayesPrint ("Included");
			
		MrBayesPrint ("    ");
			
		for (j=0; j<numDefinedConstraints; j++)
			{
			if (taxaInfo[i].constraints[j] == 0)
				MrBayesPrint ("  .");
			else
				MrBayesPrint ("  *");
			}
		MrBayesPrint ("\n");
		}
		
	MrBayesPrint ("\n");
	MrBayesPrint ("%s   Arrow indicates current outgroup\n", spacer);

	return (NO_ERROR);

}





int DoUserTree (void)

{

	int			i;
	
	if (pPtr != &userTree->nodes[0])
		{
		MrBayesPrint ("\n   ERROR: Something wrong with user tree\n");
		return (ERROR);
		}
	if (pPtr->anc == NULL)
		userTree->isRooted = YES;
	else
		{
		userTree->isRooted = NO;
		/* adjust number of nodes, since we originally thought the tree was rooted */
		userTree->nIntNodes -= 1;
		userTree->nNodes -= 2;
		}
	if (userTree->isRooted == YES)
		{
		MrBayesPrint ("%s   Defining a rooted user tree\n", spacer);
		pPtr = &userTree->nodes[nextAvailableNode];
		nextAvailableNode++;
		pPtr->left = qPtr;
		qPtr->anc = pPtr;
		userTree->root = pPtr;
		}
	else
		{
		MrBayesPrint ("%s   Defining an unrooted user tree\n", spacer);
		}
	i = numTaxa;
	if (userTree->isRooted == YES)
		FinishTree (userTree->root, &i, YES);
	else
		FinishTree (userTree->root, &i, NO);

	/* we are not going to worry here about the rooting of the user tree
	   but we may want to move the calculation root if the tree is unrooted */

	GetDownPass (userTree);
	isUserTreeDefined = YES;
	if (userTree->isRooted == NO)
		MoveCalculationRoot (userTree, outGroupNum);

	/* for debugging */
	/*
	if (userTree->isRooted == YES)
		MrBayesPrint ("\n   Tree is rooted:\n\n");
	else
		MrBayesPrint ("\n   Tree is unrooted:\n\n");
	ShowNodes (userTree->root, 0, userTree->isRooted);
	getchar();
	*/

	return (NO_ERROR);
	
}





int DoUserTreeParm (char *parmName, char *tkn)

{

	int			i, tempInt, howMany;
	MrBFlt		tempD;
	char		tempName[100];
	
	if (defMatrix == NO)
		{
		MrBayesPrint ("%s   A matrix must be specified before a user tree can be defined\n", spacer);
		return (ERROR);
		}

	if (expecting == Expecting(EQUALSIGN))
		{
		/* delete old user tree if one is present */
		if (isUserTreeDefined == YES)
			{
			FreeTree (userTree);
			}
		isUserTreeDefined = NO;
		memAllocs[ALLOC_USERTREE] = NO;
		if ((userTree = AllocateTree (numTaxa, YES)) == NULL)
			{
			MrBayesPrint ("%s   Could not allocate space for user tree\n", spacer);
			return (ERROR);
			}
		else
			memAllocs[ALLOC_USERTREE] = YES;

		/* start recording user tree */
		isFirstNode = YES;
		pPtr = qPtr = &userTree->nodes[0];
		nextAvailableNode = 0;
		userBrlensDef = NO;
		foundColon = NO;
		for (i=0; i<numTaxa; i++)
			tempSet[i] = NO;
		expecting  = Expecting(LEFTPAR);
		}
	else if (expecting == Expecting(LEFTPAR))
		{
		if (isFirstNode == YES)
			{
			pPtr = &userTree->nodes[nextAvailableNode];
			nextAvailableNode++;
			isFirstNode = NO;
			userTree->root = pPtr;
			}
		else
			{
			if (nextAvailableNode+1 >= 2*numTaxa)
				{
				MrBayesPrint ("%s   Too many nodes on user tree\n", spacer);
				return (ERROR);
				}
			if (pPtr->left == NULL)
				{
				pPtr = &userTree->nodes[nextAvailableNode];
				nextAvailableNode++;
				qPtr->left = pPtr;
				pPtr->anc = qPtr;
				qPtr = pPtr;
				}
			else if (pPtr->right == NULL)
				{
				pPtr = &userTree->nodes[nextAvailableNode];
				nextAvailableNode++;
				qPtr->right = pPtr;
				pPtr->anc = qPtr;
				qPtr = pPtr;
				}
			else if (pPtr->anc == NULL)
				{
				pPtr = &userTree->nodes[nextAvailableNode];
				nextAvailableNode++;
				qPtr->anc = pPtr;
				pPtr->left = qPtr;
				qPtr = pPtr;
				userTree->root = pPtr;
				pPtr->marked = YES;
				}
			else
				{
				MrBayesPrint ("\n   ERROR: Tree is not bifurcating\n");
				return (ERROR);
				}
			}
		expecting  = Expecting(ALPHA);
		expecting |= Expecting(NUMBER);
		expecting |= Expecting(LEFTPAR);
		}
	else if (expecting == Expecting(ALPHA))
		{
		if (nextAvailableNode+1 >= 2*numTaxa)
			{
			MrBayesPrint ("%s   Too many nodes on user tree\n", spacer);
			return (ERROR);
			}
		/* Check to see if the name is in the list of taxon names. */
		if (CheckString (tkn, taxaNames, &howMany) == ERROR)
			{
			MrBayesPrint ("%s   Could not find taxon %s in list of taxa\n", spacer, tkn);
			return (ERROR);
			}
		if (tempSet[howMany-1] == YES)
			{
			MrBayesPrint ("%s   Taxon name %s already used in tree\n", spacer, tkn);
			return (ERROR);
			}
		else
			tempSet[howMany-1] = YES;
		if (pPtr->left == NULL)
			{
			pPtr = &userTree->nodes[nextAvailableNode];
			strcpy (pPtr->label, tkn);
			pPtr->index = howMany - 1;
			nextAvailableNode++;
			qPtr->left = pPtr;
			pPtr->anc = qPtr;
			qPtr = pPtr;
			}
		else if (pPtr->right == NULL)
			{
			pPtr = &userTree->nodes[nextAvailableNode];
			strcpy (pPtr->label, tkn);
			pPtr->index = howMany - 1;
			nextAvailableNode++;
			qPtr->right = pPtr;
			pPtr->anc = qPtr;
			qPtr = pPtr;
			}
		else if (pPtr->anc == NULL)
			{
			pPtr = &userTree->nodes[nextAvailableNode];
			strcpy (pPtr->label, tkn);
			pPtr->index = howMany - 1;
			nextAvailableNode++;
			qPtr->anc = pPtr;
			pPtr->left = qPtr;
			qPtr = pPtr;
			userTree->root = pPtr;
			pPtr->marked = YES;
			}
		else
			{
			MrBayesPrint ("%s   Tree is not bifurcating\n", spacer);
			return (ERROR);
			}
		
		expecting  = Expecting(COMMA);
		expecting |= Expecting(COLON);
		expecting |= Expecting(RIGHTPAR);
		}
	else if (expecting == Expecting(RIGHTPAR))
		{
		if (pPtr->marked == NO)
			{
			if (pPtr->anc != NULL)
				{
				pPtr = pPtr->anc;
				qPtr = pPtr;
				}
			else
				{
				MrBayesPrint ("%s   Cannot go down\n", spacer);
				return (ERROR);
				}
			}
		else
			{
			if (pPtr->left != NULL)
				{
				pPtr = pPtr->left;
				qPtr = pPtr;
				}
			else
				{
				MrBayesPrint ("%s   Cannot go down\n", spacer);
				return (ERROR);
				}
			}
		expecting  = Expecting(COMMA);
		expecting |= Expecting(COLON);
		expecting |= Expecting(RIGHTPAR);
		expecting |= Expecting(SEMICOLON);
		}
	else if (expecting == Expecting(COLON))
		{
		foundColon = YES;
		expecting  = Expecting(NUMBER);
		}
	else if (expecting == Expecting(COMMA))
		{
		if (pPtr->marked == NO)
			{
			if (pPtr->anc != NULL)
				{
				pPtr = pPtr->anc;
				qPtr = pPtr;
				}
			else
				{
				MrBayesPrint ("%s   Cannot go down\n", spacer);
				return (ERROR);
				}
			}
		else
			{
			if (pPtr->left != NULL)
				{
				pPtr = pPtr->left;
				qPtr = pPtr;
				}
			else
				{
				MrBayesPrint ("%s   Cannot go down\n", spacer);
				return (ERROR);
				}
			}
		expecting  = Expecting(ALPHA);
		expecting |= Expecting(NUMBER);
		expecting |= Expecting(LEFTPAR);
		}
	else if (expecting == Expecting(NUMBER))
		{
		if (foundColon == YES)
			{
			/* branch length */
			sscanf (token, "%lf", &tempD);
			if (pPtr->marked == NO)
				pPtr->length = tempD;
			else
				{
				if (pPtr->left != NULL)
					pPtr->left->length = tempD;
				else
					{
					MrBayesPrint ("%s   Cannot assign branch length to left node\n", spacer);
					return (ERROR);
					}
				}
			userBrlensDef = YES;
			foundColon = NO;
			expecting  = Expecting(COMMA);
			expecting |= Expecting(RIGHTPAR);
			}
		else
			{
			/* taxon number */
			sscanf (token, "%d", &tempInt);

			if (nextAvailableNode+1 >= 2 * numTaxa)
				{
				MrBayesPrint ("%s   Too many nodes on user tree\n", spacer);
				return (ERROR);
				}
			/* Check to see if the name is in the list of taxon names. */
			if (CheckString (tkn, taxaNames, &howMany) == ERROR)
				{
				/* The number could not be found as a taxon name in the list of taxon names. We will
				   assume that the user has then input taxa as numbers and not the names. */
				if (tempSet[tempInt-1] == YES)
					{
					MrBayesPrint ("%s   Taxon name %d has already been used in tree\n", spacer, tempInt);
					return (ERROR);
					}
				else
					tempSet[tempInt-1] = YES;
				tempInt--;
				}
			else
				{
				howMany--;
				/* The taxon name is in the list of taxon names */
				if (howMany < 0 || howMany >= numTaxa)
					{
					MrBayesPrint ("%s   Taxon number is out of range\n", spacer);
					return (ERROR);
					}
				if (tempSet[howMany] == YES)
					{
					MrBayesPrint ("%s   Taxon %d has already been used in tree\n", spacer, howMany+1);
					return (ERROR);
					}
				else
					tempSet[howMany] = YES;
				tempInt = howMany;
				}
						
			if (GetNameFromString (taxaNames, tempName, tempInt+1) == ERROR)
				{
				MrBayesPrint ("%s   Error getting taxon name\n", spacer);
				return (ERROR);
				}
						
			if (pPtr->left == NULL)
				{
				pPtr = &userTree->nodes[nextAvailableNode];
				strcpy (pPtr->label, tempName);
				pPtr->index = tempInt;
				nextAvailableNode++;
				qPtr->left = pPtr;
				pPtr->anc = qPtr;
				qPtr = pPtr;
				}
			else if (pPtr->right == NULL)
				{
				pPtr = &userTree->nodes[nextAvailableNode];
				strcpy (pPtr->label, tempName);
				pPtr->index = tempInt;
				nextAvailableNode++;
				qPtr->right = pPtr;
				pPtr->anc = qPtr;
				qPtr = pPtr;
				}
			else if (pPtr->anc == NULL)
				{
				pPtr = &userTree->nodes[nextAvailableNode];
				strcpy (pPtr->label, tempName);
				pPtr->index = tempInt;
				nextAvailableNode++;
				qPtr->anc = pPtr;
				pPtr->left = qPtr;
				qPtr = pPtr;
				userTree->root = pPtr;
				pPtr->marked = YES;
				}
			else
				{
				MrBayesPrint ("%s   Tree is not bifurcating\n", spacer);
				return (ERROR);
				}

			expecting  = Expecting(COMMA);
			expecting |= Expecting(COLON);
			expecting |= Expecting(RIGHTPAR);
			}
		}

	return (NO_ERROR);
	MrBayesPrint ("%s", parmName); /* keep compiler happy */

}





int DoVersion (void)

{

	MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
	MrBayesPrint ("   Version                                                                       \n");
    MrBayesPrint ("                                                                                 \n");
	MrBayesPrint ("   MrBayes v%s                                                                   \n", VERSION_NUMBER);
	MrBayesPrint ("   ---------------------------------------------------------------------------   \n");

	return (NO_ERROR);
	
}





unsigned long int Expecting (int y)

{

	unsigned long int x;
	
	x = (unsigned long int)pow(2.0, (MrBFlt)y);
	
	return (x);

}

#ifdef USE_READLINE
/* This function is for commandline substitution: first word is always a command */
char *command_generator(const char *text, int state) {
	static int list_index, len;
	char *command;

	if(state==0) 
		{
		list_index=0;
		len=strlen(text);
		}
	while ((command=commands[list_index].string)!=NULL) 
		{
		list_index++;
		if (strncasecmp(command,text,len)==0) 
			/* memory is freed by the readline library so we need a strdup here */ 
			return strdup(command);
		}
	return (char *)NULL;
}
#endif



int FindValidCommand (char *tk, int *numMatches)

{

	int				i, j, tkLen, targetLen, numDiff;
	CmdType			*p;

	p = commands + 0;
	tkLen = (int) strlen(tk);

	(*numMatches) = 0;
	for (i=0; i<NUMCOMMANDS; i++)
		{
		targetLen = (int) strlen(p->string);
		if (tkLen <= targetLen)
			{
			for (j=0, numDiff=0; j<tkLen; j++)
				{
				if (ChangeCase(tk[j]) != ChangeCase(p->string[j]))
					numDiff++;
				}
			if (numDiff == 0)
				{
				(*numMatches)++;
				commandPtr = p;
				if (tkLen == targetLen)
					break;
				}
			}
		p++;
		}

	inValidCommand = NO;
	if (*numMatches == 1)
		{
		inValidCommand = YES;
		return (NO_ERROR);
		}
	else
		return (ERROR);
	
}





int FindValidParam (char *tk, int *numMatches)

{

	int			i, j, tkLen, targetLen, numDiff;
	CmdType		*p;
	ParmInfoPtr	q;

	if (commandPtr)
		p = commandPtr;
	else
		{
		MrBayesPrint ("%s   Command pointer is NULL\n", spacer);
		return (ERROR);
		}
	tkLen = (int) strlen(tk);

	*numMatches = 0;
	for (i=0; i<p->numParms; i++)
		{
		q = paramTable + (p->parmList[i]);
		targetLen = (int) strlen(q->string);
		/*printf ("%s %d (%s %d)\n", q->string, targetLen, tk, p->numParms);*/
		if (!strcmp(q->string, "Xxxxxxxxxx"))
			{
			(*numMatches)++;
			paramPtr = q;
			}
		if (tkLen <= targetLen)
			{
			for (j=0, numDiff=0; j<tkLen; j++)
				{
				if (ChangeCase(tk[j]) != ChangeCase(q->string[j]))
					numDiff++;
				}
			if (numDiff == 0)
				{
				(*numMatches)++;
				paramPtr = q;
				if (tkLen == targetLen)
					break;
				}
			}	
		}
	
	if (*numMatches == 1)
		return (NO_ERROR);
	else
		return (ERROR);
	
}





void FinishTree (TreeNode *p, int *i, int isThisTreeRooted)

{

	/* We only reindex the internal nodes of the tree. We
	   assume that the tip nodes have already been indexed
	   0, 1, 2, ..., numTaxa-1. */
	   
	if (p != NULL)
		{
		FinishTree (p->left,  i, isThisTreeRooted);
		FinishTree (p->right, i, isThisTreeRooted);
		p->marked = NO;
		if (p->left == NULL && p->right == NULL && p->anc != NULL)
			{
			}
		else if (p->left != NULL && p->right == NULL && p->anc == NULL)
			{
			if (isThisTreeRooted == YES)
				p->index = (*i)++;
			}
		else
			{
			p->index = (*i)++;
			}
		}
		
}





int FreeMatrix (void)

{

	int		memoryLetFree;
	
	memoryLetFree = NO;

	if (memAllocs[ALLOC_MATRIX] == YES)
		{
		free (matrix);
		memAllocs[ALLOC_MATRIX] = NO;
		memoryLetFree = YES;
		}
	if (memAllocs[ALLOC_CHARINFO] == YES)
		{
		free (charInfo);
		memAllocs[ALLOC_CHARINFO] = NO;
		memoryLetFree = YES;
		}
	if (memAllocs[ALLOC_TAXAINFO] == YES)
		{
		free (taxaInfo);
		memAllocs[ALLOC_TAXAINFO] = NO;
		memoryLetFree = YES;
		}
	if (memAllocs[ALLOC_PARTITIONNAMES] == YES)
		{
		free (partitionNames);
		memAllocs[ALLOC_PARTITIONNAMES] = NO;
		memoryLetFree = YES;
		}
	if (memAllocs[ALLOC_CHARSETNAMES] == YES)
		{
		free (charSetNames);
		memAllocs[ALLOC_CHARSETNAMES] = NO;
		memoryLetFree = YES;
		}
	if (memAllocs[ALLOC_TAXANAMES] == YES)
		{
		free (taxaNames);
		memAllocs[ALLOC_TAXANAMES] = NO;
		memoryLetFree = YES;
		}
	if (memAllocs[ALLOC_TMPSET] == YES)
		{
		free (tempSet);
		memAllocs[ALLOC_TMPSET] = NO;
		memoryLetFree = YES;
		}
	if (memAllocs[ALLOC_TAXASETNAMES] == YES)
		{
		free (taxaSetNames);
		memAllocs[ALLOC_TAXASETNAMES] = NO;
		memoryLetFree = YES;
		}
	if (memAllocs[ALLOC_CONSTRAINTNAMES] == YES)
		{
		free (constraintNames);
		memAllocs[ALLOC_CONSTRAINTNAMES] = NO;
		memoryLetFree = YES;
		}
	if (memAllocs[ALLOC_USERTREE] == YES)
		{
		FreeTree (userTree);
		memAllocs[ALLOC_USERTREE] = NO;
		memoryLetFree = YES;
		}
	if (memAllocs[ALLOC_TRANSFROM] == YES)
		{
		free (transFrom);
		memAllocs[ALLOC_TRANSFROM] = NO;
		memoryLetFree = YES;
		}
	if (memAllocs[ALLOC_TRANSTO] == YES)
		{
		free (transTo);
		memAllocs[ALLOC_TRANSTO] = NO;
		memoryLetFree = YES;
		}
	
	if (memoryLetFree == YES)
		MrBayesPrint ("%s   Deleting matrix\n", spacer);

	/* reinitialize everything dependent on a matrix */
	if (ReinitializeMrBayes() == ERROR)
		return (ERROR);

	return (NO_ERROR);

}





int GetNameFromString (char *s, char *tkn, int n)

{

	int		i, j, startI=0, numPrev;
	   	
	i = numPrev = 0;
	while (s[i] != '\0')
		{
		if (s[i] == '|')
			numPrev++;
		i++;
		}
	if (numPrev < n)
		{
		return (ERROR);
		}
		
	if (n == 1)
		startI = 0;
	else
		{
		i = j = 0;
		while (s[i] != '\0')
			{
			if (s[i] == '|')
				j++;
			i++;
			if (j == n - 1)
				{
				startI = i;
				break;
				}
			}
		}
		
	if (s[startI] == '\0')
		{
		MrBayesPrint ("%s   String is too full\n", spacer);
		return (ERROR);
		}
	
	i = startI;
	j = 0;
	while(s[i] != '\0' && s[i] != '|')
		{
		tkn[j++] = s[i++];
		if (s[i] == '\0')
			{
			MrBayesPrint ("%s   String is too full\n", spacer);
			return (ERROR);
			}
		}
	tkn[j] = '\0';

	return (NO_ERROR);
	
}





int GetNumPartDivisions (int n)

{

	int			i, divFound[MAX_NUM_DIVS], numDivs;
	
	for (i=0; i<MAX_NUM_DIVS; i++)
		divFound[i] = NO;
	
	for (i=0; i<numChar; i++)
		{
		divFound[charInfo[i].partitionId[n-1] - 1] = YES;
		}
		
	numDivs = 0;
	for (i=0; i<MAX_NUM_DIVS; i++)
		if (divFound[i] == YES)
			numDivs++;
		
	return (numDivs);
	
}





void GetToken (int *tokenType)

{
		
	int				allNumbers;
	register char	*temp;
	
	(*tokenType) = 0;
	temp = token;
	
	while (IsWhite(*tokenP) == 1 || IsWhite(*tokenP) == 2)
		{
		if (IsWhite(*tokenP) == 2)
			{
			*tokenType = RETURNSYMBOL;
			/*foundNewLine = YES;*/
			/* MrBayesPrint ("RETURN\n"); */
			}
		++tokenP;
		}
	
	if (readWord == YES && *tokenP != '"')
		{
		while (isgraph(*tokenP) && *tokenP!=';')
			*temp++ = *tokenP++;
		*temp = '\0';
		*tokenType = ALPHA;
		readWord = NO;
		return;
		}
	*tokenType = UNKNOWN_TOKEN_TYPE;
	if (IsIn(*tokenP,"="))
		{
		*temp++ = *tokenP++;
		*tokenType = EQUALSIGN;
		}
	else if (IsIn(*tokenP,";"))
		{
		*temp++ = *tokenP++;
		*tokenType = SEMICOLON;
		}
	else if (IsIn(*tokenP,":"))
		{
		*temp++ = *tokenP++;
		*tokenType = COLON;
		}
	else if (IsIn(*tokenP,","))
		{
		*temp++ = *tokenP++;
		*tokenType = COMMA;
		}
	else if (IsIn(*tokenP,"#"))
		{
		*temp++ = *tokenP++;
		*tokenType = POUNDSIGN;
		}
	else if (IsIn(*tokenP,"("))
		{
		*temp++ = *tokenP++;
		*tokenType = LEFTPAR;
		}
	else if (IsIn(*tokenP,")"))
		{
		*temp++ = *tokenP++;
		*tokenType = RIGHTPAR;
		}
	else if (IsIn(*tokenP,"{"))
		{
		*temp++ = *tokenP++;
		*tokenType = LEFTCURL;
		}
	else if (IsIn(*tokenP,"}"))
		{
		*temp++ = *tokenP++;
		*tokenType = RIGHTCURL;
		}
	else if (IsIn(*tokenP,"["))
		{
		*temp++ = *tokenP++;
		*tokenType = LEFTCOMMENT;
		}
	else if (IsIn(*tokenP,"]"))
		{
		*temp++ = *tokenP++;
		*tokenType = RIGHTCOMMENT;
		}
	else if (IsIn(*tokenP,"?"))
		{
		*temp++ = *tokenP++;
		*tokenType = QUESTIONMARK;
		}
	else if (IsIn(*tokenP,"-"))
		{
		*temp++ = *tokenP++;
		*tokenType = DASH;
		}
	else if (IsIn(*tokenP,"\"") && readWord == YES)
		{
		tokenP++;
		while(*tokenP != '"' && *tokenP != '\0')
			{
			*temp++ = *tokenP++;
			}
		*tokenType = ALPHA;
		*tokenP++;
		}
	else if (IsIn(*tokenP,"abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ_0123456789."))
		{
		if (IsIn(*tokenP,"0123456789."))
			allNumbers = TRUE;
		else
			allNumbers = FALSE;
		*temp++ = *tokenP++;
		while(IsIn(*tokenP,"abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ_0123456789.-"))
			{
			if(allNumbers == TRUE && *tokenP=='-') break;
			if (!IsIn(*tokenP,"0123456789."))
				allNumbers = FALSE;
			*temp++ = *tokenP++;
			}
		if (allNumbers == TRUE)
			*tokenType = NUMBER;
		else
			*tokenType = ALPHA;
		}
	else if (IsIn(*tokenP,"*"))
		{
		*temp++ = *tokenP++;
		*tokenType = ASTERISK;
		}
	else if (IsIn(*tokenP,"/"))
		{
		*temp++ = *tokenP++;
		*tokenType = FORWARDSLASH;
		}
	else if (IsIn(*tokenP,"'\\'"))
		{
		*temp++ = *tokenP++;
		*tokenType = BACKSLASH;
		}
	else if (IsIn(*tokenP,"!"))
		{
		*temp++ = *tokenP++;
		*tokenType = EXCLAMATIONMARK;
		}
	else if (IsIn(*tokenP,"%"))
		{
		*temp++ = *tokenP++;
		*tokenType = PERCENT;
		}
	else if (IsIn(*tokenP,"\""))
		{
		*temp++ = *tokenP++;
		*tokenType = QUOTATIONMARK;
		}
	else if (IsIn(*tokenP,"&~+^$@|{}`><"))
		{
		*temp++ = *tokenP++;
		*tokenType = WEIRD;
		}

	*temp = '\0';
	
}





void GetUserDownPass (TreeNode *p, TreeNode **x, int *y)

{

	if (p != NULL)
		{
		GetUserDownPass (p->left, x, y);
		GetUserDownPass (p->right, x, y);
		x[(*y)] = p;
		(*y)++;
		}
		
}





int GetUserHelp (char *helpTkn)

{

	int			i, j, tempInt;
	char		yesNoStr[4];
	
	if (!strcmp(helpTkn, "Begin"))
		{
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		MrBayesPrint ("   Begin                                                                         \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   This command is used to format data or commands in the program. The correct   \n");
	    MrBayesPrint ("   usage is                                                                      \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("      begin <data or mrbayes>;                                                   \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   The two valid uses of the \"begin\" command, then, are                        \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("      begin data;                                                                \n");
	    MrBayesPrint ("      begin mrbayes;                                                             \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   The \"data\" specifier is used to specify the beginning of a data block; your \n");
	    MrBayesPrint ("   character data should follow. For example, the following is an example of     \n");
	    MrBayesPrint ("   a data block for four taxa and ten DNA sites:                                 \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("      begin data;                                                                \n");
	    MrBayesPrint ("         dimensions ntax=4 nchar=10;                                             \n");
	    MrBayesPrint ("         format datatype=dna;                                                    \n");
	    MrBayesPrint ("         matrix                                                                  \n");
	    MrBayesPrint ("         taxon_1  AACGATTCGT                                                     \n");
	    MrBayesPrint ("         taxon_2  AAGGATTCCA                                                     \n");
	    MrBayesPrint ("         taxon_3  AACGACTCCT                                                     \n");
	    MrBayesPrint ("         taxon_4  AAGGATTCCT                                                     \n");
	    MrBayesPrint ("         ;                                                                       \n");
	    MrBayesPrint ("      end;                                                                       \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   The other commands -- dimensions, format, and matrix -- are discussed         \n");
	    MrBayesPrint ("   in the appropriate help menu. The only thing to note here is that the         \n");
	    MrBayesPrint ("   block begins with a \"begin data\" command. The \"mrbayes\" command is        \n");
	    MrBayesPrint ("   used to enter commands specific to the MrBayes program into the file.         \n");
	    MrBayesPrint ("   This allows you to automatically process commands on execution of the         \n");
	    MrBayesPrint ("   program. The following is a simple mrbayes block:                             \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("      begin mrbayes;                                                             \n");
	    MrBayesPrint ("         charset first  = 1-10\\3;                                               \n");
	    MrBayesPrint ("         charset second = 2-10\\3;                                               \n");
	    MrBayesPrint ("         charset third  = 3-10\\3;                                               \n");
	    MrBayesPrint ("      end;                                                                       \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   This mrbayes block sets off the three \"charset\" commands, used to           \n");
	    MrBayesPrint ("   predefine some blocks of characters. The mrbayes block can be very useful.    \n");
	    MrBayesPrint ("   For example, in this case, it would save you the time of typing the char-     \n");
	    MrBayesPrint ("   acter sets each time you executed the file. Also, note that every             \n");
	    MrBayesPrint ("   \"begin <data or mrbayes>\" command ends with an \"end\". Finally, you can    \n");
	    MrBayesPrint ("   have so-called foreign blocks in the file. An example of a foreign block      \n");
	    MrBayesPrint ("   would be \"begin paup\". The program will simply skip this block. This is     \n");
	    MrBayesPrint ("   useful because it means that you can use the same file for MrBayes, PAUP*     \n");
	    MrBayesPrint ("   or MacClade (although it isn't clear why you would want to use those other    \n");
	    MrBayesPrint ("   programs).                                                                    \n");
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		}
	else if (!strcmp(helpTkn, "End"))
		{
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		MrBayesPrint ("   End                                                                           \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   This command is used to terminate a data or mrbayes block. The correct        \n");
	    MrBayesPrint ("   usage is                                                                      \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("      end;                                                                       \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   For more information on this, check the help for the \"begin\" command.       \n");
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		}
	else if (!strcmp(helpTkn, "Endblock"))
		{
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		MrBayesPrint ("   Endblock                                                                      \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   This is an older, deprecated version of \"End\", see that command.            \n");
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		}
	else if (!strcmp(helpTkn, "Plot"))
		{
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		MrBayesPrint ("   Plot                                                                          \n");
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("   This command plots specified parameters in the .p file created by the         \n");
		MrBayesPrint ("   program. The program prints two files during a MCMC analysis: a tree file     \n");
		MrBayesPrint ("   and a parameter file. The parameter file has the extension \".p\".            \n");
		MrBayesPrint ("   This command, plot, makes an x-y graph of the parameter over the course       \n");
		MrBayesPrint ("   of the chain. The command can be useful for visually diagnosing convergence   \n");
		MrBayesPrint ("   for many of the parameters of the phylogenetic model. The parameter to be     \n");
		MrBayesPrint ("   plotted is specified by the \"parameter\" option. Several parameters can be   \n");
		MrBayesPrint ("   plotted at once by using the \"match\" option, which has a default value of   \n");
		MrBayesPrint ("   \"perfect\". For example, if you were to set \"parameter = pi\" and           \n");
		MrBayesPrint ("   \"match = consistentwith\", then all of the state frequency parameters would  \n");
		MrBayesPrint ("   be plotted. You can also set \"match=all\", in which case all of the          \n");
		MrBayesPrint ("   parameters are plotted.                                                       \n");
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("   Parameter       Options                      Current Setting                  \n");
		MrBayesPrint ("   ------------------------------------------------------------                  \n");
		MrBayesPrint ("   Filename        <name>                       %s                               \n", plotParams.plotFileName);
		MrBayesPrint ("   Burnin          <number>                     %d                               \n", plotParams.plotBurnIn);
		MrBayesPrint ("   Parameter       <name>                       %s                               \n", plotParams.parameter);
		MrBayesPrint ("   Match           Perfect/Consistentwith/All   %s                               \n", plotParams.match);
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		}
	else if (!strcmp(helpTkn, "Dimensions"))
		{
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		MrBayesPrint ("   Dimensions                                                                    \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   This command is used in a data block to define the number of taxa and         \n");
	    MrBayesPrint ("   characters. The correct usage is                                              \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("      dimensions ntax=<number> nchar=<number>                                    \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   The dimensions must be the first command in a data block. The following       \n");
	    MrBayesPrint ("   provides an example of the proper use of this command:                        \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("      begin data;                                                                \n");
	    MrBayesPrint ("         dimensions ntax=4 nchar=10;                                             \n");
	    MrBayesPrint ("         format datatype=dna;                                                    \n");
	    MrBayesPrint ("         matrix                                                                  \n");
	    MrBayesPrint ("         taxon_1  AACGATTCGT                                                     \n");
	    MrBayesPrint ("         taxon_2  AAGGATTCCA                                                     \n");
	    MrBayesPrint ("         taxon_3  AACGACTCCT                                                     \n");
	    MrBayesPrint ("         taxon_4  AAGGATTCCT                                                     \n");
	    MrBayesPrint ("         ;                                                                       \n");
	    MrBayesPrint ("      end;                                                                       \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   Here, the dimensions command tells MrBayes to expect a matrix with four       \n");
	    MrBayesPrint ("   taxa and 10 characters.                                                       \n");
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		}
	else if (!strcmp(helpTkn, "Format"))
		{
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		MrBayesPrint ("   Format                                                                        \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   This command is used in a data block to define the format of the char-        \n");
	    MrBayesPrint ("   acter matrix. The correct usage is                                            \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("      format datatype=<name> ... <parameter>=<option>                            \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   The format command must be the second command in a data block. The following  \n");
	    MrBayesPrint ("   provides an example of the proper use of this command:                        \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("      begin data;                                                                \n");
	    MrBayesPrint ("         dimensions ntax=4 nchar=10;                                             \n");
	    MrBayesPrint ("         format datatype=dna gap=-;                                              \n");
	    MrBayesPrint ("         matrix                                                                  \n");
	    MrBayesPrint ("         taxon_1  AACGATTCGT                                                     \n");
	    MrBayesPrint ("         taxon_2  AAGGAT--CA                                                     \n");
	    MrBayesPrint ("         taxon_3  AACGACTCCT                                                     \n");
	    MrBayesPrint ("         taxon_4  AAGGATTCCT                                                     \n");
	    MrBayesPrint ("         ;                                                                       \n");
	    MrBayesPrint ("      end;                                                                       \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   Here, the format command tells MrBayes to expect a matrix with DNA char-      \n");
	    MrBayesPrint ("   acters and with gaps coded as \"-\".                                          \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   The following are valid options for format:                                   \n");
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("   Datatype   -- This parameter MUST BE INCLUDED in the format command. More-    \n");
		MrBayesPrint ("                 over, it must be the first parameter in the line. The           \n");
		MrBayesPrint ("                 datatype command specifies what type of characters are          \n");
		MrBayesPrint ("                 in the matrix. The following are valid options:                 \n");
		MrBayesPrint ("                    Datatype = Dna: DNA states (A,C,G,T,R,Y,M,K,S,W,H,B,         \n");
		MrBayesPrint ("                               V,D,N)                                            \n");
		MrBayesPrint ("                    Datatype = Rna: DNA states (A,C,G,U,R,Y,M,K,S,W,H,B,         \n");
		MrBayesPrint ("                               V,D,N)                                            \n");
		MrBayesPrint ("                    Datatype = Protein: Amino acid states (A,R,N,D,C,Q,E,        \n");
		MrBayesPrint ("                               G,H,I,L,K,M,F,P,S,T,W,Y,V)                        \n");
		MrBayesPrint ("                    Datatype = Restriction: Restriction site (0,1) states        \n");
		MrBayesPrint ("                    Datatype = Standard: Morphological (0,1) states              \n");
		MrBayesPrint ("                    Datatype = Continuous: Real number valued states             \n");
		MrBayesPrint ("                    Datatype = Mixed(<type>:<range>,...,<type>:<range>): A       \n");
		MrBayesPrint ("                               mixture of the above datatypes. For example,      \n");
		MrBayesPrint ("                               \"datatype=mixed(dna:1-100,protein:101-200)\"     \n");
		MrBayesPrint ("                               would specify a mixture of DNA and amino acid     \n");
		MrBayesPrint ("                               characters with the DNA characters occupying      \n");
		MrBayesPrint ("                               the first 100 sites and the amino acid char-      \n");
		MrBayesPrint ("                               acters occupying the last 100 sites.              \n");
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("   Interleave -- This parameter specifies whether the data matrix is in          \n");
		MrBayesPrint ("                 interleave format. The valid options are \"Yes\" or \"No\",     \n");
		MrBayesPrint ("                 with \"No\" as the default. An interleaved matrix looks like    \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("                    format datatype=dna gap=- interleave=yes;                    \n");
	    MrBayesPrint ("                    matrix                                                       \n");
	    MrBayesPrint ("                    taxon_1  AACGATTCGT                                          \n");
	    MrBayesPrint ("                    taxon_2  AAGGAT--CA                                          \n");
	    MrBayesPrint ("                    taxon_3  AACGACTCCT                                          \n");
	    MrBayesPrint ("                    taxon_4  AAGGATTCCT                                          \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("                    taxon_1  CCTGGTAC                                            \n");
	    MrBayesPrint ("                    taxon_2  CCTGGTAC                                            \n");
	    MrBayesPrint ("                    taxon_3  ---GGTAG                                            \n");
	    MrBayesPrint ("                    taxon_4  ---GGTAG                                            \n");
	    MrBayesPrint ("                    ;                                                            \n");
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("   Gap        -- This parameter specifies the format for gaps. Note that         \n");
		MrBayesPrint ("                 gap character can only be a single character and that it        \n");
		MrBayesPrint ("                 cannot correspond to a standard state (e.g., A,C,G,T,R,Y,       \n");
		MrBayesPrint ("                 M,K,S,W,H,B,V,D,N for nucleotide data).                         \n");
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("   Missing    -- This parameter specifies the format for missing data. Note      \n");
		MrBayesPrint ("                 that the missing character can only be a single character and   \n");
		MrBayesPrint ("                 cannot correspond to a standard state (e.g., A,C,G,T,R,Y,       \n");
		MrBayesPrint ("                 M,K,S,W,H,B,V,D,N for nucleotide data). This is often an        \n");
		MrBayesPrint ("                 unnecessary parameter to set because many data types, such      \n");
		MrBayesPrint ("                 as nucleotide or amino acid, already have a missing char-       \n");
		MrBayesPrint ("                 acter specified. However, for morphological or restriction      \n");
		MrBayesPrint ("                 site data, \"missing=?\" is often used to specify ambiguity     \n");
		MrBayesPrint ("                 or unobserved data.                                             \n");
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("   Matchchar  -- This parameter specifies the matching character for the         \n");
		MrBayesPrint ("                 matrix. For example,                                            \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("                    format datatype=dna gap=- matchchar=.;                       \n");
	    MrBayesPrint ("                    matrix                                                       \n");
	    MrBayesPrint ("                    taxon_1  AACGATTCGT                                          \n");
	    MrBayesPrint ("                    taxon_2  ..G...--CA                                          \n");
	    MrBayesPrint ("                    taxon_3  .....C..C.                                          \n");
	    MrBayesPrint ("                    taxon_4  ..G.....C.                                          \n");
	    MrBayesPrint ("                    ;                                                            \n");
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("                 is equivalent to                                                \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("                    format datatype=dna gap=-;                                   \n");
	    MrBayesPrint ("                    matrix                                                       \n");
	    MrBayesPrint ("                    taxon_1  AACGATTCGT                                          \n");
	    MrBayesPrint ("                    taxon_2  AAGGAT--CA                                          \n");
	    MrBayesPrint ("                    taxon_3  AACGACTCCT                                          \n");
	    MrBayesPrint ("                    taxon_4  AAGGATTCCT                                          \n");
	    MrBayesPrint ("                    ;                                                            \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   The only non-standard NEXUS format option is the use of the \"mixed\",        \n");
	    MrBayesPrint ("   \"restriction\", \"standard\" and \"continuous\" datatypes. Hence, if         \n");
	    MrBayesPrint ("   you use any of these datatype specifiers, a program like PAUP* or             \n");
	    MrBayesPrint ("   MacClade will report an error (as they should because MrBayes is not          \n");
	    MrBayesPrint ("   strictly NEXUS compliant).                                                    \n");
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		}
	else if (!strcmp(helpTkn, "Matrix"))
		{
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		MrBayesPrint ("   Matrix                                                                        \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   This command specifies the actual data for the phylogenetic analysis.         \n");
	    MrBayesPrint ("   The character matrix should follow the dimensions and format commands         \n");
	    MrBayesPrint ("   in a data block. The matrix can have all of the characters for a taxon        \n");
	    MrBayesPrint ("   on a single line:                                                             \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("      begin data;                                                                \n");
	    MrBayesPrint ("         dimensions ntax=4 nchar=10;                                             \n");
	    MrBayesPrint ("         format datatype=dna gap=-;                                              \n");
	    MrBayesPrint ("         matrix                                                                  \n");
	    MrBayesPrint ("         taxon_1  AACGATTCGT                                                     \n");
	    MrBayesPrint ("         taxon_2  AAGGAT--CA                                                     \n");
	    MrBayesPrint ("         taxon_3  AACGACTCCT                                                     \n");
	    MrBayesPrint ("         taxon_4  AAGGATTCCT                                                     \n");
	    MrBayesPrint ("         ;                                                                       \n");
	    MrBayesPrint ("      end;                                                                       \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   or be in \"interleaved\" format:                                              \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("      begin data;                                                                \n");
	    MrBayesPrint ("         dimensions ntax=4 nchar=20;                                             \n");
	    MrBayesPrint ("         format datatype=dna gap=- interleave=yes;                               \n");
	    MrBayesPrint ("         matrix                                                                  \n");
	    MrBayesPrint ("         taxon_1  AACGATTCGT                                                     \n");
	    MrBayesPrint ("         taxon_2  AAGGAT--CA                                                     \n");
	    MrBayesPrint ("         taxon_3  AACGACTCCT                                                     \n");
	    MrBayesPrint ("         taxon_4  AAGGATTCCT                                                     \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("         taxon_1  TTTTCGAAGC                                                     \n");
	    MrBayesPrint ("         taxon_2  TTTTCGGAGC                                                     \n");
	    MrBayesPrint ("         taxon_3  TTTTTGATGC                                                     \n");
	    MrBayesPrint ("         taxon_4  TTTTCGGAGC                                                     \n");
	    MrBayesPrint ("         ;                                                                       \n");
	    MrBayesPrint ("      end;                                                                       \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   Note that the taxon names must not have spaces. If you really want to         \n");
	    MrBayesPrint ("   indicate a space in a taxon name (perhaps between a genus and species         \n");
	    MrBayesPrint ("   name), then you might use an underline (\"_\"). There should be at            \n");
	    MrBayesPrint ("   least a single space after the taxon name, separating the name from           \n");
	    MrBayesPrint ("   the actual data on that line. There can be spaces between the char-           \n");
	    MrBayesPrint ("   acters.                                                                       \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   If you have mixed data, then you specify all of the data in the same          \n");
	    MrBayesPrint ("   matrix. Here is an example that includes two different data types:            \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("      begin data;                                                                \n");
	    MrBayesPrint ("         dimensions ntax=4 nchar=20;                                             \n");
	    MrBayesPrint ("         format datatype=mixed(dna:1-10,standard:21-30) interleave=yes;          \n");
	    MrBayesPrint ("         matrix                                                                  \n");
	    MrBayesPrint ("         taxon_1  AACGATTCGT                                                     \n");
	    MrBayesPrint ("         taxon_2  AAGGAT--CA                                                     \n");
	    MrBayesPrint ("         taxon_3  AACGACTCCT                                                     \n");
	    MrBayesPrint ("         taxon_4  AAGGATTCCT                                                     \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("         taxon_1  0001111111                                                     \n");
	    MrBayesPrint ("         taxon_2  0111110000                                                     \n");
	    MrBayesPrint ("         taxon_3  1110000000                                                     \n");
	    MrBayesPrint ("         taxon_4  1000001111                                                     \n");
	    MrBayesPrint ("         ;                                                                       \n");
	    MrBayesPrint ("      end;                                                                       \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   The matrix command is terminated by a semicolon.                              \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   Finally, just a note on data presentation. It is much easier for others       \n");
	    MrBayesPrint ("   to (1) understand your data and (2) repeat your analyses if you make          \n");
	    MrBayesPrint ("   your data clean, comment it liberally (using the square brackets), and        \n");
	    MrBayesPrint ("   embed the commands you used in a publication in the mrbayes block.            \n");
	    MrBayesPrint ("   Remember that the data took a long time for you to collect. You might         \n");
	    MrBayesPrint ("   as well spend a little time making the data file look nice and clear to       \n");
	    MrBayesPrint ("   any that may later request the data for further analysis.                     \n");
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		}
	else if (!strcmp(helpTkn, "Pairs"))
		{
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		MrBayesPrint ("   Pairs                                                                         \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   This command is used to specify pairs of nucleotides. For example, your       \n");
	    MrBayesPrint ("   data may be RNA sequences with a known secondary structure of stems and       \n");
	    MrBayesPrint ("   loops. Substitutions in nucleotides involved in a Watson-Crick pairing        \n");
	    MrBayesPrint ("   in stems are not strictly independent; a change in one changes the prob-      \n");
	    MrBayesPrint ("   ability of a change in the partner. A solution to this problem is to          \n");
	    MrBayesPrint ("   expand the model around the pair of nucleotides in the stem. This             \n");
	    MrBayesPrint ("   command allows you to do this. The correct usage is:                          \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("      pairs <NUC1>:<NUC2>, <NUC1>:<NUC2>,..., <NUC1>:<NUC2>;                     \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   For example,                                                                  \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("      pairs 30:56, 31:55, 32:54, 33:53, 34:52, 35:51, 36:50;                     \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   specifies pairings between nucleotides 30 and 56, 31 and 55, etc. Only        \n");
	    MrBayesPrint ("   nucleotide data (DNA or RNA) may be paired using this command. Note that      \n");
	    MrBayesPrint ("   in order for the program to actually implement a \"doublet\" model            \n");
	    MrBayesPrint ("   involving a 16 X 16 rate matrix, you must specify that the structure of       \n");
	    MrBayesPrint ("   the model is 16 X 16 using \"lset nucmodel=doublet\".                         \n");
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		}
	else if (!strcmp(helpTkn, "Databreaks"))
		{
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		MrBayesPrint ("   Databreaks                                                                    \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   This command is used to specify breaks in your input data matrix. Your        \n");
	    MrBayesPrint ("   data may be a mixture of genes or a mixture of different types of data.       \n");
	    MrBayesPrint ("   Some of the models implemented by MrBayes account for nonindependence at      \n");
	    MrBayesPrint ("   adjacent characters. The autocorrelated gamma model, for example, allows      \n");
	    MrBayesPrint ("   rates at adjacent sites to be correlated. However, there is no way for        \n");
	    MrBayesPrint ("   such a model to tell whether two sites, adjacent in the matrix, are           \n");
	    MrBayesPrint ("   actually separated by many kilobases or megabases in the genome. The          \n");
	    MrBayesPrint ("   databreaks command allows you to specify such breaks. The correct             \n");
	    MrBayesPrint ("   usage is:                                                                     \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("      databreaks <break 1> <break 2> <break 3> ...                               \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   For example, say you have a data matrix of 3204 characters that include       \n");
	    MrBayesPrint ("   nucleotide data from three genes. The first gene covers characters 1 to       \n");
	    MrBayesPrint ("   970, the second gene covers characters 971 to 2567, and the third gene        \n");
	    MrBayesPrint ("   covers characters 2568 to 3204. Also, let's assume that the genes are         \n");
	    MrBayesPrint ("   not directly adjacent to one another in the genome, as might be likely        \n");
	    MrBayesPrint ("   if you have mitochondrial sequences. In this case, you can specify            \n");
	    MrBayesPrint ("   breaks between the genes using:                                               \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("      databreaks 970 2567;                                                       \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   The first break, between genes one and two, is after character 970 and        \n");
	    MrBayesPrint ("   the second break, between genes two and three, is after character 2567.       \n");
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		}
	else if (!strcmp(helpTkn, "Acknowledgments"))
		{
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		MrBayesPrint ("   Acknowledgments                                                               \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   This command shows the authors' acknowledgments.                              \n");
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		}
	else if (!strcmp(helpTkn, "About"))
		{
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		MrBayesPrint ("   About                                                                         \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   This command provides some general information about the program.             \n");
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		}
	else if (!strcmp(helpTkn, "Version"))
		{
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		MrBayesPrint ("   Version                                                                       \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   This command shows the release version of the program.                        \n");
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		}
	else if (!strcmp(helpTkn, "Citations"))
		{
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		MrBayesPrint ("   Citations                                                                     \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   This command shows a thorough list of citations you may consider using        \n");
	    MrBayesPrint ("   when publishing the results of a MrBayes analysis.                            \n");
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		}
	else if (!strcmp(helpTkn, "Showmatrix"))
		{
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		MrBayesPrint ("   Showmatrix                                                                    \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   This command shows the character matrix currently in memory.                  \n");
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		}
	else if (!strcmp(helpTkn, "Constraint"))
		{
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		MrBayesPrint ("   Constraint                                                                    \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   This command defines a tree constraint. The format for the constraint         \n");
	    MrBayesPrint ("   command is                                                                    \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("      constraint <constraint name> <probability> = <list of taxa>                \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   A list of taxa can be specified using a taxset, taxon names, or taxon         \n");
	    MrBayesPrint ("   numbers. A probability must also be specified. For now, MrBayes ignores       \n");
	    MrBayesPrint ("   this probability value and treats the constraint as an absolute requirement   \n");
	    MrBayesPrint ("   of trees. That is, trees that are not compatible with the constraint have     \n");
	    MrBayesPrint ("   zero prior (and hence zero posterior) probability.                            \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   Future releases of MrBayes will use the probability value to determine how    \n");
	    MrBayesPrint ("   much more probable a tree is that contains the constraint than a tree without \n");
	    MrBayesPrint ("   the constraint. For example, the following command                            \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("      constraint example 100 = taxon_2 taxon_3                                   \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   defines a constraint called \"example\" that includes two taxa. In future     \n");
	    MrBayesPrint ("   releases of MrBayes, trees that contain a clade with those two taxa together  \n");
	    MrBayesPrint ("   will have a prior probability that is 100 times that of trees without the     \n");
	    MrBayesPrint ("   constraint. In the current version, the probability value will be ignored     \n");
	    MrBayesPrint ("   and trees without the two taxa together will not be sampled.                  \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   If you are interested in inferring ancestral states for a particular node,    \n");
	    MrBayesPrint ("   you need to constrain that node first using the 'constraint' command. For     \n");
	    MrBayesPrint ("   more information on how to infer ancestral states, see the help for the       \n");
	    MrBayesPrint ("   'report' command.                                                             \n");

		/*
		MrBayesPrint ("   An alternative way to specify the probability is to use \"exp(<number>)\"        \n");
	    MrBayesPrint ("   instead. For example, the following command                                   \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("      constraint example exp(10) = taxon_2 taxon_3                               \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   defines a constraint called \"example\" that includes two taxa. Trees         \n");
	    MrBayesPrint ("   that contain the clade are exp(10), or 22026.46, times more probable than     \n");
	    MrBayesPrint ("   trees not containing the clade.                                               \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   In principal, the use of constraints allows you to specify very compl-        \n");
	    MrBayesPrint ("   icated prior probabilities on trees. It is hoped that this prior more         \n");
	    MrBayesPrint ("   faithfully represents your prior beliefs in trees. Using several con-         \n");
	    MrBayesPrint ("   straints simultaneously can create a very complicated prior on trees. For     \n");
	    MrBayesPrint ("   example, consider the following two constraints:                              \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("      constraint A 100 =  3 4 5                                                  \n");
	    MrBayesPrint ("      constraint B 50 = 4 5                                                      \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   which can be summarized in a table (if there are 10 taxa total)               \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("       1  2  3  4  5  6  7  8  9 10  Times more probable                         \n");
	    MrBayesPrint ("       ----------------------------  -------------------                         \n");
	    MrBayesPrint ("    A  .  .  +  +  +  .  .  .  .  .         100                                  \n");
	    MrBayesPrint ("    B  .  .  .  +  +  .  .  .  .  .          50                                  \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   Trees that contain taxa 3, 4, and 5 are 100 times more probable than          \n");
	    MrBayesPrint ("   trees not containing the group. Similarly, trees that contain taxa 4 and      \n");
	    MrBayesPrint ("   5 are 50 times more probable than those without. Calculating the rela-        \n");
	    MrBayesPrint ("   tive probabilties of any two trees is relatively easy. For example, say       \n");
	    MrBayesPrint ("   we have two trees, T1 and T2, that have the following properties:             \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("      Tree  Constraint A    Constraint B                                         \n");
	    MrBayesPrint ("      ----------------------------------                                         \n");
	    MrBayesPrint ("      T1    Consistent      Consistent                                           \n");
	    MrBayesPrint ("      T2    Inconsistent    Consistent                                           \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   Both trees are consistent with constraint B but T1 is consistent with         \n");
	    MrBayesPrint ("   constraint A whereas T2 is not. Therefore, T1 is 100 times more probable      \n");
	    MrBayesPrint ("   than T2. We only need to know the relative probabilies of the two trees       \n");
	    MrBayesPrint ("   (i.e., the probabilities up to a constant) because the MCMC will allow        \n");
	    MrBayesPrint ("   us to normalize the probabilities. Life is somewhat more complicated          \n");
	    MrBayesPrint ("   when the trees differ at two or more constraints. For example, say the        \n");
	    MrBayesPrint ("   two trees have the following properties:                                      \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("      Tree  Constraint A    Constraint B                                         \n");
	    MrBayesPrint ("      ----------------------------------                                         \n");
	    MrBayesPrint ("      T1    Consistent      Consistent                                           \n");
	    MrBayesPrint ("      T2    Inconsistent    Inconsistent                                         \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   Now T1 is 50 X 100 = 5000 times more probable than T2. You can see that       \n");
	    MrBayesPrint ("   by simply defining two constraints, we have defined a fairly complicated      \n");
	    MrBayesPrint ("   prior probability distribution on trees:                                      \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("       |                                                                         \n");
	    MrBayesPrint ("       |                             11              Trees with A and B          \n");
	    MrBayesPrint ("       |                                                                         \n");
	    MrBayesPrint ("    R  |                                                                         \n");
	    MrBayesPrint ("    e  |                                                                         \n");
	    MrBayesPrint ("    l  |                         2222  2222          Trees with A and not B      \n");
	    MrBayesPrint ("       |                                                                         \n");
	    MrBayesPrint ("    P  |                 33333333          33333333  Trees with B but not A      \n");
	    MrBayesPrint ("    r  |                                                                         \n");
	    MrBayesPrint ("    o  |                                                                         \n");
	    MrBayesPrint ("    b  |                                                                         \n");
	    MrBayesPrint ("       | 4444444444444444                          4444444444444444              \n");
	    MrBayesPrint ("       |                                                                         \n");
	    MrBayesPrint ("       |                                                                         \n");
	    MrBayesPrint ("       |_____________________________________________________________            \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("                                    Trees                                        \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   There are four levels on this landscape of trees. The peak (1) contains       \n");
	    MrBayesPrint ("   trees that are consistent with both constraints. The valley (4) has trees     \n");
	    MrBayesPrint ("   that contain neither constraint. There are two plateaus of trees that         \n");
	    MrBayesPrint ("   contain one but not both constraints. Note that the relative probab-          \n");
	    MrBayesPrint ("   ilities between levels are as follows:                                        \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("      T1 is 5000 times more probable than T4                                     \n");
	    MrBayesPrint ("      T1 is  100 times more probable than T3                                     \n");
	    MrBayesPrint ("      T1 is   50 times more probable than T2                                     \n");
	    MrBayesPrint ("      T2 is  100 times more probable than T4                                     \n");
	    MrBayesPrint ("      T2 is    2 times more probable than T3                                     \n");
	    MrBayesPrint ("      T3 is   50 times more probable than T4                                     \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   Conversely,                                                                   \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("      T4 is 1/5000th as probable as T1                                           \n");
	    MrBayesPrint ("      T3 is  1/100th as probable as T1                                           \n");
	    MrBayesPrint ("      T2 is   1/50th as probable as T1                                           \n");
	    MrBayesPrint ("      T4 is  1/100th as probable as T2                                           \n");
	    MrBayesPrint ("      T3 is      1/2 as probable as T2                                           \n");
	    MrBayesPrint ("      T4 is   1/50th as probable as T3                                           \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   To get the relative probabilities between levels 1 and 4, you can             \n");
	    MrBayesPrint ("   simply multiply the probabilities to successive plateaus:                     \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("       P(1,4) = P(1,2) X P(2,3) X P(3,4) = 50 X 2 X 50 = 5000                    \n");
	    */
		MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   It is important to note that simply defining a constraint using this          \n");
	    MrBayesPrint ("   command is not sufficient for the program to actually implement the           \n");
	    MrBayesPrint ("   constraint in an analysis. You must also specify the constraints using        \n");
	    MrBayesPrint ("   'prset topologypr = constraints (<name of constraint>)'. For more infor-      \n");
	    MrBayesPrint ("   mation, see the help on the 'prset' command.                                  \n");
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		}
	else if (!strcmp(helpTkn, "Showmodel"))
		{
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		MrBayesPrint ("   Showmodel                                                                     \n");
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("   This command shows the current model settings. The correct usage is           \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("      showmodel                                                                  \n");
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("   After typing \"showmodel\", the modelling assumptions are shown on a          \n");
		MrBayesPrint ("   partition-by-partition basis.                                                 \n");
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		}
	else if (!strcmp(helpTkn, "Execute"))
		{
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		MrBayesPrint ("   Execute                                                                       \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   This command executes a file called <file name>. The correct usage is:        \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("      execute <file name>                                                        \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   For example,                                                                  \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("      execute replicase.nex                                                      \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   would execute the file named \"replicase.nex\". This file must be in the      \n");
	    MrBayesPrint ("   same directory as the executable.                                             \n");
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		}
	else if (!strcmp(helpTkn, "Lset"))
		{
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		MrBayesPrint ("   Lset                                                                          \n");
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("   This command sets the parameters of the likelihood model. The likelihood      \n");
		MrBayesPrint ("   function is the probability of observing the data conditional on the phylo-   \n");
		MrBayesPrint ("   genetic model. In order to calculate the likelihood, you must assume a        \n");
		MrBayesPrint ("   model of character change. This command lets you tailor the biological        \n");
		MrBayesPrint ("   assumptions made in the phylogenetic model. The correct usage is              \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("      lset <parameter>=<option> ... <parameter>=<option>                         \n");
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("   For example, \"lset nst=6 rates=gamma\" would set the model to a general      \n");
		MrBayesPrint ("   model of DNA substition (the GTR) with gamma-distributed rate variation       \n");
		MrBayesPrint ("   across sites.                                                                 \n");
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("   Options:                                                                      \n");
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("   Applyto   -- This option allows you to apply the lset commands to specific    \n");
		MrBayesPrint ("                partitions. This command should be the first in the list of      \n");
		MrBayesPrint ("                commands specified in lset. Moreover, it only makes sense to     \n");
		MrBayesPrint ("                be using this command if the data have been partitioned. A       \n");
		MrBayesPrint ("                default partition is set on execution of a matrix. If the data   \n");
		MrBayesPrint ("                are homogeneous (i.e., all of the same data type), then this     \n");
		MrBayesPrint ("                partition will not subdivide the characters. Up to 30 other      \n");
		MrBayesPrint ("                partitions can be defined, and you can switch among them using   \n");
		MrBayesPrint ("                \"set partition=<partition name>\". Now, you may want to         \n");
		MrBayesPrint ("                specify different models to different partitions of the data.    \n");
		MrBayesPrint ("                Applyto allows you to do this. For example, say you have         \n");
		MrBayesPrint ("                partitioned the data by codon position, and you want to apply    \n");
		MrBayesPrint ("                a nst=2 model to the first two partitions and nst=6 to the       \n");
		MrBayesPrint ("                last. This could be implemented in two uses of lset:             \n");
		MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("                   lset applyto=(1,2) nst=2                                      \n");
		MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("                   lset applyto=(3) nst=6                                        \n");
		MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("                The first applies the parameters after \"applyto\" to the        \n");
		MrBayesPrint ("                first and second partitions. The second lset applies nst=6       \n");
		MrBayesPrint ("                to the third partition. You can also use applyto=(all), which    \n");
		MrBayesPrint ("                attempts to apply the parameter settings to all of the data      \n");
		MrBayesPrint ("                partitions. Importantly, if the option is not consistent with    \n");
		MrBayesPrint ("                the data in the partition, the program will not apply the        \n");
		MrBayesPrint ("                lset option to that partition.                                   \n");
		MrBayesPrint ("   Nucmodel  -- This specifies the general form of the nucleotide substitution   \n");
		MrBayesPrint ("                model. The options are \"4by4\" [the standard model of DNA       \n");
		MrBayesPrint ("                substitution in which there are only four states (A,C,G,T/U)],   \n");
		MrBayesPrint ("                \"doublet\" (a model appropriate for modelling the stem regions  \n");
		MrBayesPrint ("                of ribosomal genes where the state space is the 16 doublets of   \n");
		MrBayesPrint ("                nucleotides), and \"codon\" (the substitution model is expanded  \n");
		MrBayesPrint ("                around triplets of nucleotides--a codon).                        \n");
		MrBayesPrint ("   Nst       -- Sets the number of substitution types: \"1\" constrains all of   \n");
		MrBayesPrint ("                the rates to be the same (e.g., a JC69 or F81 model); \"2\" all- \n");
		MrBayesPrint ("                ows transitions and transversions to have potentially different  \n");
		MrBayesPrint ("                rates (e.g., a K80 or HKY85 model); \"6\" allows all rates to    \n");
		MrBayesPrint ("                be different, subject to the constraint of time-reversibility    \n");
		MrBayesPrint ("                (e.g., a GTR model).                                             \n");
		MrBayesPrint ("   Code      -- Enforces the use of a particular genetic code. The default       \n");
		MrBayesPrint ("                is the universal code. Other options include \"vertmt\" for      \n");
		MrBayesPrint ("                vertebrate mitocondrial DNA, \"mycoplasma\", \"yeast\",          \n");
		MrBayesPrint ("                \"ciliates\", and \"metmt\" (for metazoan mitochondrial DNA      \n");
		MrBayesPrint ("                except vertebrates).                                             \n");
		MrBayesPrint ("   Ploidy    -- Specifies the ploidy of the organism. Options are \"Haploid\"    \n");
		MrBayesPrint ("                or \"Diploid\". This option is used when a coalescence prior     \n");
		MrBayesPrint ("                is used on trees.                                                \n");
		MrBayesPrint ("   Rates     -- Sets the model for among-site rate variation. In general, the    \n");
		MrBayesPrint ("                rate at a site is considered to be an unknown random variable.   \n");
		MrBayesPrint ("                The valid options are:                                           \n");
		MrBayesPrint ("                * equal    -- No rate variation across sites.                    \n");
		MrBayesPrint ("                * gamma    -- Gamma-distributed rates across sites. The rate     \n");
		MrBayesPrint ("                              at a site is drawn from a gamma distribution.      \n");
		MrBayesPrint ("                              The gamma distribution has a single parameter      \n");
		MrBayesPrint ("                              that describes how much rates vary.                \n");
		MrBayesPrint ("                * adgamma  -- Autocorrelated rates across sites. The marg-       \n");
		MrBayesPrint ("                              inal rate distribution is gamma, but adjacent      \n");
		MrBayesPrint ("                              sites have correlated rates.                       \n");
		MrBayesPrint ("                * propinv  -- A proportion of the sites are invariable.          \n");
		MrBayesPrint ("                * invgamma -- A proportion of the sites are invariable while     \n");
		MrBayesPrint ("                              the rate for the remaining sites are drawn from    \n");
		MrBayesPrint ("                              a gamma distribution.                              \n");
		MrBayesPrint ("                Note that MrBayes versions 2.0 and earlier supported options     \n");
		MrBayesPrint ("                that allowed site specific rates (e.g., ssgamma). In versions    \n");
		MrBayesPrint ("                3.0 and later, site specific rates are allowed, but set using    \n");
		MrBayesPrint ("                the 'prset ratepr' command for each partition.                   \n");
		MrBayesPrint ("   Ngammacat -- Sets the number of rate categories for the gamma distribution.   \n");
		MrBayesPrint ("                The gamma distribution is continuous. However, it is virtually   \n");
		MrBayesPrint ("                impossible to calculate likelihoods under the continuous gamma   \n");
		MrBayesPrint ("                distribution. Hence, an approximation to the continuous gamma    \n");
		MrBayesPrint ("                is used; the gamma distribution is broken into ncat categories   \n");
		MrBayesPrint ("                of equal weight (1/ncat). The mean rate for each category rep-   \n");
		MrBayesPrint ("                resents the rate for the entire cateogry. This option allows     \n");
		MrBayesPrint ("                you to specify how many rate categories to use when approx-      \n");
		MrBayesPrint ("                imating the gamma. The approximation is better as ncat is inc-   \n");
		MrBayesPrint ("                reased. In practice, \"ncat=4\" does a reasonable job of         \n");
		MrBayesPrint ("                approximating the continuous gamma.                              \n");
		MrBayesPrint ("   Nbetacat  -- Sets the number of rate categories for the beta distribution.    \n");
		MrBayesPrint ("                A symmetric beta distribution is used to model the station-      \n");
		MrBayesPrint ("                ary frequencies when morphological data are used. This option    \n");
		MrBayesPrint ("                specifies how well the beta distribution will be approx-         \n");
		MrBayesPrint ("                imated.                                                          \n");
		MrBayesPrint ("   Omegavar  -- Allows the nonsynonymous/synonymous rate ratio (omega) to vary   \n");
		MrBayesPrint ("                across codons. Ny98 assumes that there are three classes, with   \n");
		MrBayesPrint ("                potentially different omega values (omega1, omega2, omega3):     \n");
		MrBayesPrint ("                omega2 = 1; 0 < omega1 <1; and omega3 > 1. Like the Ny98 model,  \n");
		MrBayesPrint ("                the M3 model has three omega classes. However, their values are  \n");
		MrBayesPrint ("                less constrained, with omega1 < omega2 < omega3. The default     \n");
		MrBayesPrint ("                (omegavar = equal) has no variation on omega across sites.       \n");
		MrBayesPrint ("   Covarion  -- This forces the use of a covarion-like model of substitution     \n");
		MrBayesPrint ("                for nucleotide or amino acid data. The valid options are \"yes\" \n");
		MrBayesPrint ("                and \"no\". The covarion model allows the rate at a site to      \n");
		MrBayesPrint ("                change over its evolutionary history. Specifically, the site     \n");
		MrBayesPrint ("                is either on or off. When it is off, no substitutions are poss-  \n");
		MrBayesPrint ("                ible. When the process is on, substitutions occur according to   \n");
		MrBayesPrint ("                a specified substitution model (specified using the other        \n");
		MrBayesPrint ("                lset options).                                                   \n");
		MrBayesPrint ("   Coding    -- This specifies how characters were sampled. If all site pat-     \n");
		MrBayesPrint ("                terns had the possibility of being sampled, then \"all\" should  \n");
		MrBayesPrint ("                be specified (the default). Otherwise \"variable\" (only var-    \n");
		MrBayesPrint ("                iable characters had the possibility of being sampled),          \n");
		MrBayesPrint ("                \"noabsence\" (characters for which all taxa were coded as       \n");
		MrBayesPrint ("                absent were not sampled), and \"nopresence\" (characters for     \n");
		MrBayesPrint ("                which all taxa were coded as present were not sampled. \"All\"   \n");
		MrBayesPrint ("                works for all data types. However, the others only work for      \n");
		MrBayesPrint ("                morphological (all/variable) or restriction site (all/variable/  \n");
		MrBayesPrint ("                noabsence/nopresence) data.                                      \n");
		MrBayesPrint ("   Parsmodel -- This forces calculation under the so-called parsimony model      \n");
		MrBayesPrint ("                described by Tuffley and Steel (1998). The options are \"yes\"   \n");
		MrBayesPrint ("                or \"no\". Note that the biological assumptions of this model    \n");
		MrBayesPrint ("                are anything but parsimonious. In fact, this model assumes many  \n");
		MrBayesPrint ("                more parameters than the next most complicated model imple-      \n");
		MrBayesPrint ("                mented in this program. If you really believe that the pars-     \n");
		MrBayesPrint ("                imony model makes the biological assumptions described by        \n");
		MrBayesPrint ("                Tuffley and Steel, then the parsimony method is miss-named.      \n");
		/*MrBayesPrint ("   Augment   -- This allows the chain to consider the missing entries of         \n");
		MrBayesPrint ("                the data matrix as random variables. A Gibbs sampler is          \n");
		MrBayesPrint ("                used to sample states.                                           \n");*/
	    MrBayesPrint ("                                                                                 \n");
	    if (numCurrentDivisions == 0)
	    	tempInt = 1;
	    else
	    	tempInt = numCurrentDivisions;
	    for (i=0; i<tempInt; i++)
	    	{
		    if (numCurrentDivisions == 0)
				MrBayesPrint ("   Default model settings:                                                       \n");
			else
				MrBayesPrint ("   Model settings for partition %d:                                              \n", i+1);
	    	MrBayesPrint ("                                                                                 \n");
			MrBayesPrint ("   Parameter    Options                               Current Setting            \n");
			MrBayesPrint ("   ------------------------------------------------------------------            \n");		
			MrBayesPrint ("   Nucmodel     4by4/Doublet/Codon                    %s                         \n", modelParams[i].nucModel);
			MrBayesPrint ("   Nst          1/2/6                                 %s                         \n", modelParams[i].nst);
			MrBayesPrint ("   Code         Universal/Vertmt/Mycoplasma/                                     \n");
			MrBayesPrint ("                Yeast/Ciliates/Metmt                  %s                         \n", modelParams[i].geneticCode);
			MrBayesPrint ("   Ploidy       Haploid/Diploid                       %s                         \n", modelParams[i].ploidy);
			MrBayesPrint ("   Rates        Equal/Gamma/Propinv/Invgamma/Adgamma  %s                         \n", modelParams[i].ratesModel);
			MrBayesPrint ("   Ngammacat    <number>                              %d                         \n", modelParams[i].numGammaCats);
			MrBayesPrint ("   Nbetacat     <number>                              %d                         \n", modelParams[i].numBetaCats);
			MrBayesPrint ("   Omegavar     Equal/Ny98/M3                         %s                         \n", modelParams[i].omegaVar);
			MrBayesPrint ("   Covarion     No/Yes                                %s                         \n", modelParams[i].covarionModel);
			MrBayesPrint ("   Coding       All/Variable/Noabsencesites/                                     \n");
			MrBayesPrint ("                Nopresencesites                       %s                         \n", modelParams[i].coding);
			MrBayesPrint ("   Parsmodel    No/Yes                                %s                         \n", modelParams[i].parsModel);
			/*MrBayesPrint ("   Augment      No/Yes                                %s                         \n", modelParams[i].augmentData);*/
			MrBayesPrint ("   ------------------------------------------------------------------            \n");		
	    	MrBayesPrint ("                                                                                 \n");
			}
		}
	else if (!strcmp(helpTkn, "Prset"))
		{
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		MrBayesPrint ("   Prset                                                                         \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   This command sets the priors for the phylogenetic model. Remember that        \n");
	    MrBayesPrint ("   in a Bayesian analysis, you must specify a prior probability distribution     \n");
	    MrBayesPrint ("   for the parameters of the likelihood model. The prior distribution rep-       \n");
	    MrBayesPrint ("   resents your prior beliefs about the parameter before observation of the      \n");
	    MrBayesPrint ("   data. This command allows you to tailor your prior assumptions to a large     \n");
	    MrBayesPrint ("   extent.                                                                       \n");
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("   Options:                                                                      \n");
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("   Applyto       -- This option allows you to apply the prset commands to        \n");
		MrBayesPrint ("                    specific partitions. This command should be the first        \n");
		MrBayesPrint ("                    in the list of commands specified in prset. Moreover, it     \n");
		MrBayesPrint ("                    only makes sense to be using this command if the data        \n");
		MrBayesPrint ("                    have been partitioned. A default partition is set on         \n");
		MrBayesPrint ("                    execution of a matrix. If the data are homogeneous           \n");
		MrBayesPrint ("                    (i.e., all of the same data type), then this partition       \n");
		MrBayesPrint ("                    will not subdivide the characters. Up to 30 other part-      \n");
		MrBayesPrint ("                    itions can be defined, and you can switch among them using   \n");
		MrBayesPrint ("                    \"set partition=<partition name>\". Now, you may want to     \n");
		MrBayesPrint ("                    specify different priors to different partitions of the      \n");
		MrBayesPrint ("                    data. Applyto allows you to do this. For example, say        \n");
		MrBayesPrint ("                    you have partitioned the data by codon position, and         \n");
		MrBayesPrint ("                    you want to fix the statefreqs to equal for the first two    \n");
		MrBayesPrint ("                    partitions but apply a flat Dirichlet prior to the state-    \n");
		MrBayesPrint ("                    freqs of the last. This could be implemented in two uses of  \n");
		MrBayesPrint ("                    prset:                                                       \n");
		MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("                       prset applyto=(1,2) statefreqs=fixed(equal)               \n");
		MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("                       prset applyto=(3) statefreqs=dirichlet(1,1,1,1)           \n");
		MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("                    The first applies the parameters after \"applyto\"           \n");
		MrBayesPrint ("                    to the first and second partitions. The second prset         \n");
		MrBayesPrint ("                    applies a flat Dirichlet to the third partition. You can     \n");
		MrBayesPrint ("                    also use applyto=(all), which attempts to apply the para-    \n");
		MrBayesPrint ("                    meter settings to all of the data partitions. Importantly,   \n");
		MrBayesPrint ("                    if the option is not consistent with the data in the part-   \n");
		MrBayesPrint ("                    ition, the program will not apply the prset option to        \n");
		MrBayesPrint ("                    that partition.                                              \n");
		MrBayesPrint ("   Tratiopr      -- This parameter sets the prior for the transition/trans-      \n");
		MrBayesPrint ("                    version rate ratio (tratio). The options are:                \n");
		MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("                       prset tratiopr = beta(<number>, <number>)                 \n");
		MrBayesPrint ("                       prset tratiopr = fixed(<number>)                          \n");
		MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("                    The program assumes that the transition and transversion     \n");
		MrBayesPrint ("                    rates are independent gamma-distributed random variables     \n");
		MrBayesPrint ("                    with the same scale parameter when beta is selected. If you  \n");
		MrBayesPrint ("                    want a diffuse prior that puts equal emphasis on transition/ \n");
		MrBayesPrint ("                    transversion rate ratios above 1.0 and below 1.0, then use a \n");
		MrBayesPrint ("                    flat Beta, beta(1,1), which is the default. If you wish to   \n");
		MrBayesPrint ("                    concentrate this distribution more in the equal-rates region,\n");
		MrBayesPrint ("                    then use a prior of the type beta(x,x), where the magnitude  \n");
		MrBayesPrint ("                    of x determines how much the prior is concentrated in the    \n");
		MrBayesPrint ("                    equal rates region. For instance, a beta(20,20) puts more    \n");
		MrBayesPrint ("                    probability on rate ratios close to 1.0 than a beta(1,1). If \n");
		MrBayesPrint ("                    you think it is likely that the transition/transversion rate \n");
		MrBayesPrint ("                    ratio is 2.0, you can use a prior of the type beta(2x,x),    \n");
		MrBayesPrint ("                    where x determines how strongly the prior is concentrated on \n");
		MrBayesPrint ("                    tratio values near 2.0. For instance, a beta(2,1) is much    \n");
		MrBayesPrint ("                    more diffuse than a beta(80,40) but both have the expected   \n");
		MrBayesPrint ("                    tratio 2.0 in the absence of data. The parameters of the     \n");
		MrBayesPrint ("                    Beta can be interpreted as counts: if you have observed x    \n");
		MrBayesPrint ("                    transitions and y transversions, then a beta(x+1,y+1) is a   \n");
		MrBayesPrint ("                    good representation of this information. The fixed option    \n");
		MrBayesPrint ("                    allows you to fix the tratio to a particular value.          \n");
		MrBayesPrint ("   Revmatpr      -- This parameter sets the prior for the substitution rates     \n");
		MrBayesPrint ("                    of the GTR model for nucleotide data. The options are:       \n");
		MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("                       prset revmatpr = dirichlet(<number>,<number>,...,<number>)\n");
		MrBayesPrint ("                       prset revmatpr = fixed(<number>,<number>,...,<number>)    \n");
		MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("                    The program assumes that the six substitution rates          \n");
		MrBayesPrint ("                    are independent gamma-distributed random variables with the  \n");
		MrBayesPrint ("                    same scale parameter when dirichlet is selected. The six     \n");
		MrBayesPrint ("                    numbers in brackets each corresponds to a particular substi- \n");
		MrBayesPrint ("                    tution type. Together, they determine the shape of the prior.\n");
		MrBayesPrint ("                    The six rates are in the order A<->C, A<->G, A<->T, C<->G,   \n");
		MrBayesPrint ("                    C<->T, and G<->T. If you want an uninformative prior you can \n");
		MrBayesPrint ("                    use dirichlet(1,1,1,1,1,1), also referred to as a 'flat'     \n");
		MrBayesPrint ("                    Dirichlet. This is the default setting. If you wish a prior  \n");
		MrBayesPrint ("                    where the C<->T rate is 5 times and the A<->G rate 2 times   \n");
		MrBayesPrint ("                    higher, on average, than the transversion rates, which are   \n");
		MrBayesPrint ("                    all the same, then you should use a prior of the form        \n");
		MrBayesPrint ("                    dirichlet(x,2x,x,x,5x,x), where x determines how much the    \n");
		MrBayesPrint ("                    prior is focused on these particular rates. For more info,   \n");
		MrBayesPrint ("                    see tratiopr. The fixed option allows you to fix the substi- \n");
		MrBayesPrint ("                    tution rates to particular values.                           \n");
		MrBayesPrint ("   Aamodelpr     -- This parameter sets the rate matrix for amino acid data.     \n");
		MrBayesPrint ("                    You can either fix the model by specifying aamodelpr=        \n");
		MrBayesPrint ("                    fixed(<model name>), where <model name> is 'poisson' (a      \n");
		MrBayesPrint ("                    glorified Jukes-Cantor model), 'jones', 'dayhoff', 'mtrev',  \n");
		MrBayesPrint ("                    'mtmam', 'wag', 'rtrev', 'cprev', 'vt', 'blosum', 'equalin'  \n");
		MrBayesPrint ("                    (a glorified Felsenstein 1981 model), or 'gtr'. You can also \n");
		MrBayesPrint ("                    average over the first ten models by specifying aamodelpr=   \n");
		MrBayesPrint ("                    mixed. If you do so, the Markov chain will sample each model \n");
		MrBayesPrint ("                    according to its probability. The sampled model is reported  \n");
		MrBayesPrint ("                    as an index: poisson(0), jones(1), dayhoff(2), mtrev(3),     \n");
		MrBayesPrint ("                    mtmam(4), wag(5), rtrev(6), cprev(7), vt(8), or blosum(9).   \n");
		MrBayesPrint ("                    The 'Sump' command summarizes the MCMC samples and calculates\n");
		MrBayesPrint ("                    the posterior probability estimate for each of these models. \n");
		MrBayesPrint ("   Aarevmatpr    -- This parameter sets the prior for the substitution rates     \n");
		MrBayesPrint ("                    of the GTR model for amino acid data. The options are:       \n");
		MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("                       prset revmatpr = dirichlet(<number>,<number>,...,<number>)\n");
		MrBayesPrint ("                       prset revmatpr = fixed(<number>,<number>,...,<number>)    \n");
		MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("                    The options are the same as those for 'Revmatpr' except that \n");
		MrBayesPrint ("                    they are defined over the 190 rates of the time-reversible   \n");
		MrBayesPrint ("                    GTR model for amino acids instead of over the 6 rates of the \n");
		MrBayesPrint ("                    GTR model for nucleotides. The rates are in the order A<->R, \n");
		MrBayesPrint ("                    A<->N, etc to Y<->V. In other words, amino acids are listed  \n");
		MrBayesPrint ("                    in alphabetic order based on their full name. The first amino\n");
		MrBayesPrint ("                    acid (Alanine) is then combined in turn with all amino acids \n");
		MrBayesPrint ("                    following it in the list, starting with amino acid 2 (Argi-  \n");
		MrBayesPrint ("                    nine) and finishing with amino acid 20 (Valine). The second  \n");
		MrBayesPrint ("                    amino acid (Arginine) is then combined in turn with all amino\n");
		MrBayesPrint ("                    acids following it, starting with amino acid 3 (Asparagine)  \n");
		MrBayesPrint ("                    and finishing with amino acid 20 (Valine), and so on.        \n");
		MrBayesPrint ("   Omegapr       -- This parameter specifies the prior on the nonsynonymous/     \n");
		MrBayesPrint ("                    synonymous rate ratio. The options are:                      \n");
		MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("                       prset omegapr = uniform(<number>,<number>)                \n");
		MrBayesPrint ("                       prset omegapr = exponential(<number>)                     \n");
		MrBayesPrint ("                       prset omegapr = fixed(<number>)                           \n");
		MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("                    This parameter is only in effect if the nucleotide sub-      \n");
		MrBayesPrint ("                    stitution model is set to codon using the lset command       \n");
		MrBayesPrint ("                    (lset nucmodel=codon). Moreover, it only applies to the      \n");
		MrBayesPrint ("                    case when there is no variation in omega across sites (i.e., \n");
		MrBayesPrint ("                    \"lset omegavar=equal\").                                    \n");
		MrBayesPrint ("   Ny98omega1pr  -- This parameter specifies the prior on the nonsynonymous/     \n");
		MrBayesPrint ("                    synonymous rate ratio for sites under purifying selection.   \n");
		MrBayesPrint ("                    The options are:                                             \n");
		MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("                       prset Ny98omega1pr = beta(<number>,<number>)              \n");
		MrBayesPrint ("                       prset Ny98omega1pr = fixed(<number>)                      \n");
		MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("                    This parameter is only in effect if the nucleotide sub-      \n");
		MrBayesPrint ("                    stitution model is set to codon using the lset command       \n");
		MrBayesPrint ("                    (lset nucmodel=codon). Moreover, it only applies to the      \n");
		MrBayesPrint ("                    case where omega varies across sites using the model of      \n");
		MrBayesPrint ("                    Nielsen and Yang (1998) (i.e., \"lset omegavar=ny98\"). If   \n");
		MrBayesPrint ("                    fixing the parameter, you must specify a number between      \n");
		MrBayesPrint ("                    0 and 1.                                                     \n");
		MrBayesPrint ("   Ny98omega3pr  -- This parameter specifies the prior on the nonsynonymous/     \n");
		MrBayesPrint ("                    synonymous rate ratio for positively selected sites. The     \n");
		MrBayesPrint ("                    options are:                                                 \n");
		MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("                       prset Ny98omega3pr = uniform(<number>,<number>)           \n");
		MrBayesPrint ("                       prset Ny98omega3pr = exponential(<number>)                \n");
		MrBayesPrint ("                       prset Ny98omega3pr = fixed(<number>)                      \n");
		MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("                    This parameter is only in effect if the nucleotide sub-      \n");
		MrBayesPrint ("                    stitution model is set to codon using the lset command       \n");
		MrBayesPrint ("                    (lset nucmodel=codon). Moreover, it only applies to the      \n");
		MrBayesPrint ("                    case where omega varies across sites according to the        \n");
		MrBayesPrint ("                    NY98 model. Note that if the NY98 model is specified         \n");
		MrBayesPrint ("                    that this parameter must be greater than 1, so you should    \n");
		MrBayesPrint ("                    not specify a uniform(0,10) prior, for example.              \n");
		MrBayesPrint ("   M3omegapr     -- This parameter specifies the prior on the nonsynonymous/     \n");
		MrBayesPrint ("                    synonymous rate ratios for all three classes of sites for    \n");
		MrBayesPrint ("                    the M3 model. The options are:                               \n");
		MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("                       prset M3omegapr = exponential                             \n");
		MrBayesPrint ("                       prset M3omegapr = fixed(<number>,<number>,<number>)       \n");
		MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("                    This parameter is only in effect if the nucleotide sub-      \n");
		MrBayesPrint ("                    stitution model is set to codon using the lset command       \n");
		MrBayesPrint ("                    (lset nucmodel=codon). Moreover, it only applies to the      \n");
		MrBayesPrint ("                    case where omega varies across sites using the M3 model of   \n");
		MrBayesPrint ("                    Yang et al. (2000) (i.e., \"lset omegavar=M3\"). Under the   \n");
		MrBayesPrint ("                    exponential prior, the four rates (dN1, dN2, dN3, and dS)    \n");
		MrBayesPrint ("                    are all considered to be independent draws from the same     \n");
		MrBayesPrint ("                    exponential distribution (the parameter of the exponential   \n");
		MrBayesPrint ("                    does not matter, and so you don't need to specify it). The   \n");
		MrBayesPrint ("                    rates dN1, dN2, and dN3 are taken to be the order statistics \n");
		MrBayesPrint ("                    with dN1 < dN2 < dN3. These three rates are all scaled to    \n");
		MrBayesPrint ("                    the same synonymous rate, dS. The other option is to simply  \n");
		MrBayesPrint ("                    fix the three rate ratios to some values.                    \n");
		MrBayesPrint ("   Codoncatfreqs -- This parameter specifies the prior on frequencies of sites   \n");
		MrBayesPrint ("                    under purifying, neutral, and positive selection. The        \n");
		MrBayesPrint ("                    options are:                                                 \n");
		MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("                       prset codoncatfreqs = dirichlet(<num>,<num>,<num>)        \n");
		MrBayesPrint ("                       prset codoncatfreqs = fixed(<number>,<number>,<number>)   \n");
		MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("                    This parameter is only in effect if the nucleotide sub-      \n");
		MrBayesPrint ("                    stitution model is set to codon using the lset command       \n");
		MrBayesPrint ("                    (lset nucmodel=codon). Moreover, it only applies to the      \n");
		MrBayesPrint ("                    case where omega varies across sites using the models of     \n");
		MrBayesPrint ("                    Nielsen and Yang (1998) (i.e., \"lset omegavar=ny98\")       \n");
		MrBayesPrint ("                    or Yang et al. (2000) (i.e., \"lset omegavar=M3\")           \n");
		MrBayesPrint ("                    Note that the sum of the three frequencies must be 1.        \n");
		MrBayesPrint ("   Statefreqpr   -- This parameter specifies the prior on the state freq-        \n");
		MrBayesPrint ("                    uencies. The options are:                                    \n");
		MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("                       prset statefreqpr = dirichlet(<number>)                   \n");
		MrBayesPrint ("                       prset statefreqpr = dirichlet(<number>,...,<number>)      \n");
		MrBayesPrint ("                       prset statefreqpr = fixed(equal)                          \n");
		MrBayesPrint ("                       prset statefreqpr = fixed(empirical)                      \n");
		MrBayesPrint ("                       prset statefreqpr = fixed(<number>,...,<number>)          \n");
		MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("                    For the dirichlet, you can specify either a single number    \n");
		MrBayesPrint ("                    or as many numbers as there are states. If you specify a     \n");
		MrBayesPrint ("                    single number, then the prior has all states equally         \n");
		MrBayesPrint ("                    probable with a variance related to the single parameter     \n");
		MrBayesPrint ("                    passed in.                                                   \n");
		MrBayesPrint ("   Shapepr       -- This parameter specifies the prior for the gamma shape       \n");
		MrBayesPrint ("                    parameter for among-site rate variation. The options are:    \n");
		MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("                       prset shapepr = uniform(<number>,<number>)                \n");
		MrBayesPrint ("                       prset shapepr = exponential(<number>)                     \n");
		MrBayesPrint ("                       prset shapepr = fixed(<number>)                           \n");
		MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("   Pinvarpr      -- This parameter specifies the prior for the proportion of     \n");
		MrBayesPrint ("                    invariable sites. The options are:                           \n");
		MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("                       prset pinvarpr = uniform(<number>,<number>)               \n");
		MrBayesPrint ("                       prset pinvarpr = fixed(<number>)                          \n");
		MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("                    Note that the valid range for the parameter is between 0     \n");
		MrBayesPrint ("                    and 1. Hence, \"prset pinvarpr=uniform(0,0.8)\" is valid     \n");
		MrBayesPrint ("                    while \"prset pinvarpr=uniform(0,10)\" is not. The def-      \n");
		MrBayesPrint ("                    ault setting is \"prset pinvarpr=uniform(0,1)\".             \n");
		MrBayesPrint ("   Ratecorrpr    -- This parameter specifies the prior for the autocorrelation   \n");
		MrBayesPrint ("                    parameter of the autocorrelated gamma distribution for       \n");
		MrBayesPrint ("                    among-site rate variation. The options are:                  \n");
		MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("                       prset ratecorrpr = uniform(<number>,<number>)             \n");
		MrBayesPrint ("                       prset ratecorrpr = fixed(<number>)                        \n");
		MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("                    Note that the valid range for the parameter is between -1    \n");
		MrBayesPrint ("                    and 1. Hence, \"prset ratecorrpr=uniform(-1,1)\" is valid    \n");
		MrBayesPrint ("                    while \"prset ratecorrpr=uniform(-11,10)\" is not. The       \n");
		MrBayesPrint ("                    default setting is \"prset ratecorrpr=uniform(-1,1)\".       \n");
		MrBayesPrint ("   Covswitchpr   -- This option sets the prior for the covarion switching        \n");
		MrBayesPrint ("                    rates. The options are:                                      \n");
		MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("                       prset covswitchpr = uniform(<number>,<number>)            \n");
		MrBayesPrint ("                       prset covswitchpr = exponential(<number>)                 \n");
		MrBayesPrint ("                       prset covswitchpr = fixed(<number>,<number>)              \n");
		MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("                    The covarion model has two rates: a rate from on to off      \n");
		MrBayesPrint ("                    and a rate from off to on. The rates are assumed to have     \n");
		MrBayesPrint ("                    independent priors that individually are either uniformly    \n");
		MrBayesPrint ("                    or exponentially distributed. The other option is to         \n");
		MrBayesPrint ("                    fix the switching rates, in which case you must specify      \n");
		MrBayesPrint ("                    both rates. (The first number is off->on and the second      \n");
		MrBayesPrint ("                    is on->off).                                                 \n");
		MrBayesPrint ("   Symdirihyperpr - This option sets the prior for the stationary frequencies   \n");
		MrBayesPrint ("                    of the states for morphological (standard) data. There can   \n");
		MrBayesPrint ("                    be as many as 10 states for standard data. However, the      \n");
		MrBayesPrint ("                    labelling of the states is somewhat arbitrary. For example,  \n");
		MrBayesPrint ("                    the state \"1\" for different characters does not have the   \n");
		MrBayesPrint ("                    same meaning. This is not true for DNA characters, for ex-   \n");
		MrBayesPrint ("                    ample, where a \"G\" has the same meaning across characters. \n");
		MrBayesPrint ("                    The fact that the labelling of morphological characters is   \n");
		MrBayesPrint ("                    arbitrary makes it difficult to allow unequal character-     \n");
		MrBayesPrint ("                    state frequencies. MrBayes gets around this problem by       \n");
		MrBayesPrint ("                    assuming that the states have a dirichlet prior, with all    \n");
		MrBayesPrint ("                    states having equal frequency. The variation in the diri-    \n");
		MrBayesPrint ("                    chlet can be controlled by this parameter--symdirihyperpr.   \n");
		MrBayesPrint ("                    Symdirihyperpr specifies the distribution on the variance    \n");
		MrBayesPrint ("                    parameter of the dirichlet. The valid options are:           \n");
		MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("                       prset Symdirihyperpr = uniform(<number>,<number>)         \n");
		MrBayesPrint ("                       prset Symdirihyperpr = exponential(<number>)              \n");
		MrBayesPrint ("                       prset Symdirihyperpr = fixed(<number>)                    \n");
		MrBayesPrint ("                       prset Symdirihyperpr = fixed(infinity)                    \n");
		MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("                    If \"fixed(infinity)\" is chosen, the dirichlet prior is     \n");
		MrBayesPrint ("                    fixed such that all character states have equal frequency.   \n");
		MrBayesPrint ("   Topologypr    -- This parameter specifies the prior probabilities of          \n");
		MrBayesPrint ("                    phylogenies. The options are:                                \n");
		MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("                       prset topologypr = uniform                                \n");
		MrBayesPrint ("                       prset topologypr = constraints(<list>)                    \n");
		MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("                    If the prior is selected to be \"uniform\", the default,     \n");
		MrBayesPrint ("                    then all possible trees are considered a priori equally      \n");
		MrBayesPrint ("                    probable. The constraints option allows you to specify       \n");
		MrBayesPrint ("                    complicated prior probabilities on trees (constraints        \n");
		MrBayesPrint ("                    are discussed more fully in \"help constraint\"). Note       \n");
		MrBayesPrint ("                    that you must specify a list of constraints that you         \n");
		MrBayesPrint ("                    wish to be obeyed. The list can be either the constraints'   \n");
		MrBayesPrint ("                    name or number. Also, note that the constraints simply       \n");
		MrBayesPrint ("                    tell you how much more (or less) probable individual         \n");
		MrBayesPrint ("                    trees are that possess the constraint than trees not         \n");
		MrBayesPrint ("                    possessing the constraint.                                   \n");
		MrBayesPrint ("   Brlenspr      -- This parameter specifies the prior probability dist-         \n");
		MrBayesPrint ("                    ribution on branch lengths. The options are:                 \n");
		MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("                       prset brlenspr = unconstrained:uniform(<num>,<num>)       \n");
		MrBayesPrint ("                       prset brlenspr = unconstrained:exponential(<number>)      \n");
		MrBayesPrint ("                       prset brlenspr = clock:uniform                            \n");
		MrBayesPrint ("                       prset brlenspr = clock:birthdeath                         \n");
		MrBayesPrint ("                       prset brlenspr = clock:coalescence                        \n");
		MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("                    Trees with unconstrained branch lengths are unrooted         \n");
		MrBayesPrint ("                    whereas clock-constrained trees are rooted. The option       \n");
		MrBayesPrint ("                    after the colon specifies the details of the probability     \n");
		MrBayesPrint ("                    density of branch lengths. If you choose a birth-death       \n");
		MrBayesPrint ("                    or coalescence prior, you may want to modify the details     \n");
		MrBayesPrint ("                    of the parameters of those processes.                        \n");
		MrBayesPrint ("   Treeheightpr  -- This parameter specifies the prior probability dist-         \n");
		MrBayesPrint ("                    ribution on the tree height, when a clock model is           \n");
		MrBayesPrint ("                    specified. The options are:                                  \n");
		MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("                       prset treeheightpr = Gamma(<num>,<num>)                   \n");
		MrBayesPrint ("                       prset treeheightpr = Exponential(<number>)                \n");
		MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("                    (And, yes, we know the exponential is a special case of      \n");
		MrBayesPrint ("                    the gamma distribution.) The tree height is the expected     \n");
		MrBayesPrint ("                    number of substitutions on a single branch that extends      \n");
		MrBayesPrint ("                    from the root of the tree to the tips. This parameter        \n");
		MrBayesPrint ("                    does not come into play for the coalescence prior. It        \n");
		MrBayesPrint ("                    insures that the prior probability distribution for          \n");
		MrBayesPrint ("                    unconstrained and birth-death models is proper.              \n");
		MrBayesPrint ("   Ratepr        -- This parameter allows you to specify the site specific       \n");
		MrBayesPrint ("                    rates model. First, you must have defined a partition of     \n");
		MrBayesPrint ("                    the characters. For example, you may define a partition      \n");
		MrBayesPrint ("                    that divides the characters by codon position, if you        \n");
		MrBayesPrint ("                    have DNA data. Second, you must make that partition the      \n");
		MrBayesPrint ("                    active one using the set command. For example, if your       \n");
		MrBayesPrint ("                    partition is called \"by_codon\", then you make that the     \n");
		MrBayesPrint ("                    active partition using \"set partition=by_codon\". Now       \n");
		MrBayesPrint ("                    that you have defined and activated a partition, you can     \n");
		MrBayesPrint ("                    specify the rate multipliers for the various partitions.     \n");
		MrBayesPrint ("                    The options are:                                             \n");
		MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("                       prset ratepr = fixed                                      \n");
		MrBayesPrint ("                       prset ratepr = variable                                   \n");
		MrBayesPrint ("                       prset ratepr = dirichlet(<number>,<number>,...,<number>)  \n");
		MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("                    If you specify \"fixed\", then the rate multiplier for       \n");
		MrBayesPrint ("                    that partition is set to 1 (i.e., the rate is fixed to       \n");
		MrBayesPrint ("                    the average rate across partitions). On the other hand,      \n");
		MrBayesPrint ("                    if you specify \"variable\", then the rate is allowed to     \n");
		MrBayesPrint ("                    vary across partitions subject to the constraint that the    \n");
		MrBayesPrint ("                    average rate of substitution across the partitions is 1.     \n");
		MrBayesPrint ("                    You must specify a variable rate prior for at least two      \n");
		MrBayesPrint ("                    partitions, otherwise the option is not activated when       \n");
		MrBayesPrint ("                    calculating likelihoods. The variable option automatically   \n");
		MrBayesPrint ("                    associates the partition rates with a dirichlet(1,...,1)     \n");
		MrBayesPrint ("                    prior. The dirichlet option is an alternative way of setting \n");
		MrBayesPrint ("                    a partition rate to be variable, and also gives accurate     \n");
		MrBayesPrint ("                    control of the shape of the prior. The parameters of the     \n");
		MrBayesPrint ("                    Dirichlet are listed in the order of the partitions that the \n");
		MrBayesPrint ("                    ratepr is applied to. For instance, \"prset applyto=(1,3,4)  \n");
		MrBayesPrint ("                    ratepr = dirichlet(10,40,15)\" would set the Dirichlet para- \n");
		MrBayesPrint ("                    meter 10 to partition 1, 40 to partition 3, and 15 to parti- \n");
		MrBayesPrint ("                    tion 4.                                                      \n");
		MrBayesPrint ("   Speciationpr  -- This parameter sets the prior on the speciation rate. The    \n");
		MrBayesPrint ("                    options are:                                                 \n");
		MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("                       prset speciationpr = uniform(<number>,<number>)           \n");
		MrBayesPrint ("                       prset speciationpr = exponential(<number>)                \n");
		MrBayesPrint ("                       prset speciationpr = fixed(<number>)                      \n");
		MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("                    This parameter is only relevant if the birth-death           \n");
		MrBayesPrint ("                    process is selected as the prior on branch lengths.          \n");
		MrBayesPrint ("   Extinctionpr  -- This parameter sets the prior on the extinction rate. The    \n");
		MrBayesPrint ("                    options are:                                                 \n");
		MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("                       prset extinctionpr = uniform(<number>,<number>)           \n");
		MrBayesPrint ("                       prset extinctionpr = exponential(<number>)                \n");
		MrBayesPrint ("                       prset extinctionpr = fixed(<number>)                      \n");
		MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("                    This parameter is only relevant if the birth-death           \n");
		MrBayesPrint ("                    process is selected as the prior on branch lengths.          \n");
		MrBayesPrint ("   Sampleprob    -- This parameter sets the fraction of species that are         \n");
		MrBayesPrint ("                    sampled in the analysis. This is used with the birth-        \n");
		MrBayesPrint ("                    death prior on trees (see Yang and Rannala, 1997).           \n");
		MrBayesPrint ("   Thetapr       -- This parameter sets the prior on the coalescence para-       \n");
		MrBayesPrint ("                    meter. The options are:                                      \n");
		MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("                       prset thetapr = uniform(<number>,<number>)                \n");
		MrBayesPrint ("                       prset thetapr = exponential(<number>)                     \n");
		MrBayesPrint ("                       prset thetapr = fixed(<number>)                           \n");
		MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("                    This parameter is only relevant if the coalescence           \n");
		MrBayesPrint ("                    process is selected as the prior on branch lengths.          \n");
	/*	MrBayesPrint ("   Growthpr      -- This parameter sets the prior on the exponential growth      \n");
		MrBayesPrint ("                    parameter of the coalescence process. The options are:       \n");
		MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("                       prset growthpr = uniform(<number>,<number>)               \n");
		MrBayesPrint ("                       prset growthpr = exponential(<number>)                    \n");
		MrBayesPrint ("                       prset growthpr = fixed(<number>)                          \n");
		MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("                    This parameter is only relevant if the coalescence           \n");
		MrBayesPrint ("                    process is selected as the prior on branch lengths.          \n");*/
	    MrBayesPrint ("                                                                                 \n");
	    if (numCurrentDivisions == 0)
	    	tempInt = 1;
	    else
	    	tempInt = numCurrentDivisions;
	    for (i=0; i<tempInt; i++)
	    	{
		    if (numCurrentDivisions == 0)
				MrBayesPrint ("   Default model settings:                                                       \n");
			else
				MrBayesPrint ("   Model settings for partition %d:                                              \n", i+1);
	    	MrBayesPrint ("                                                                                 \n");
			MrBayesPrint ("   Parameter        Options                      Current Setting                 \n");
			MrBayesPrint ("   ------------------------------------------------------------------            \n");		

			MrBayesPrint ("   Tratiopr         Beta/Fixed                   %s", modelParams[i].tRatioPr);
			if (!strcmp(modelParams[i].tRatioPr, "Beta"))
				MrBayesPrint ("(%1.1lf,%1.1lf)\n", modelParams[i].tRatioDir[0], modelParams[i].tRatioDir[1]);
			else
				MrBayesPrint ("(%1.1lf)\n", modelParams[i].tRatioFix);

			MrBayesPrint ("   Revmatpr         Dirichlet/Fixed              %s", modelParams[i].revMatPr);
			if (!strcmp(modelParams[i].revMatPr, "Dirichlet"))
				MrBayesPrint ("(%1.1lf,%1.1lf,%1.1lf,%1.1lf,%1.1lf,%1.1lf)\n", modelParams[i].revMatDir[0],
				modelParams[i].revMatDir[1], modelParams[i].revMatDir[2], modelParams[i].revMatDir[3],
				modelParams[i].revMatDir[4], modelParams[i].revMatDir[5]);
			else
				MrBayesPrint ("(%1.1lf,%1.1lf,%1.1lf,%1.1lf,%1.1lf,%1.1lf)\n", modelParams[i].revMatFix[0],
				modelParams[i].revMatFix[1], modelParams[i].revMatFix[2], modelParams[i].revMatFix[3],
				modelParams[i].revMatFix[4], modelParams[i].revMatFix[5]);

			MrBayesPrint ("   Aamodelpr        Fixed/Mixed                  %s", modelParams[i].aaModelPr);
			if (!strcmp(modelParams[i].aaModelPr, "Fixed"))
				MrBayesPrint ("(%s)\n", modelParams[i].aaModel);
			else
				MrBayesPrint ("\n");

			MrBayesPrint ("   Aarevmatpr       Dirichlet/Fixed              %s", modelParams[i].aaRevMatPr);
			if (!strcmp(modelParams[i].revMatPr, "Dirichlet"))
				{
				for (j=1; j<190; j++)
					if (AreDoublesEqual (modelParams[i].aaRevMatDir[0], modelParams[i].aaRevMatDir[j], 0.00001) == NO)
						break;
				if (j==190)
					MrBayesPrint ("(%1.1lf,%1.1lf,...)\n", modelParams[i].revMatDir[0], modelParams[i].revMatDir[0]);
				else
					MrBayesPrint (" (use 'Showmodel' to see values set by user)\n");
				}
			else
				{
				for (j=1; j<190; j++)
					if (AreDoublesEqual (modelParams[i].aaRevMatFix[0], modelParams[i].aaRevMatFix[j], 0.00001) == NO)
						break;
				if (j==190)
					MrBayesPrint ("(%1.1lf,%1.1lf,...)\n", modelParams[i].revMatFix[0], modelParams[i].revMatFix[0]);
				else
					MrBayesPrint (" (use 'Showmodel' to see values set by user)\n");
				}

			MrBayesPrint ("   Omegapr          Dirichlet/Fixed              %s", modelParams[i].omegaPr);
			if (!strcmp(modelParams[i].omegaPr, "Dirichlet"))
				MrBayesPrint ("(%1.1lf,%1.1lf)\n", modelParams[i].omegaDir[0], modelParams[i].omegaDir[1]);
			else
				MrBayesPrint ("(%1.1lf)\n", modelParams[i].omegaFix);

			MrBayesPrint ("   Ny98omega1pr     Beta/Fixed                   %s", modelParams[i].ny98omega1pr);
			if (!strcmp(modelParams[i].ny98omega1pr, "Beta"))
				MrBayesPrint ("(%1.1lf,%1.1lf)\n", modelParams[i].ny98omega1Beta[0], modelParams[i].ny98omega1Beta[1]);
			else if (!strcmp(modelParams[i].ny98omega1pr, "Fixed"))
				MrBayesPrint ("(%1.1lf)\n", modelParams[i].ny98omega1Fixed);
				
			MrBayesPrint ("   Ny98omega3pr     Uniform/Exponential/Fixed    %s", modelParams[i].ny98omega3pr);
			if (!strcmp(modelParams[i].ny98omega3pr, "Uniform"))
				MrBayesPrint ("(%1.1lf,%1.1lf)\n", modelParams[i].ny98omega3Uni[0], modelParams[i].ny98omega3Uni[1]);
			else if (!strcmp(modelParams[i].ny98omega3pr, "Exponential"))
				MrBayesPrint ("(%1.1lf)\n", modelParams[i].ny98omega3Exp);
			else
				MrBayesPrint ("(%1.1lf)\n", modelParams[i].ny98omega3Fixed);

			MrBayesPrint ("   M3omegapr        Exponential/Fixed            %s", modelParams[i].m3omegapr);
			if (!strcmp(modelParams[i].m3omegapr, "Exponential"))
				MrBayesPrint ("\n");
			else if (!strcmp(modelParams[i].m3omegapr, "Fixed"))
				MrBayesPrint ("(%1.1lf,%1.1lf,%1.1lf)\n", modelParams[i].m3omegaFixed[0], modelParams[i].m3omegaFixed[1], modelParams[i].m3omegaFixed[2]);
				
			MrBayesPrint ("   Codoncatfreqs    Dirichlet/Fixed              %s", modelParams[i].codonCatFreqPr);
			if (!strcmp(modelParams[i].codonCatFreqPr, "Dirichlet"))
				MrBayesPrint ("(%1.1lf,%1.1lf,%1.1lf)\n", modelParams[i].codonCatDir[0], modelParams[i].codonCatDir[1], modelParams[i].codonCatDir[2]);
			else
				MrBayesPrint ("(%1.1lf,%1.1lf,%1.1lf)\n", modelParams[i].codonCatFreqFix[0], modelParams[i].codonCatFreqFix[1], modelParams[i].codonCatFreqFix[2]);

			MrBayesPrint ("   Statefreqpr      Dirichlet/Fixed              %s", modelParams[i].stateFreqPr);
			if (!strcmp(modelParams[i].stateFreqPr, "Dirichlet"))
				{
				if (modelParams[i].dataType == DNA || modelParams[i].dataType == RNA)
					{
					if (!strcmp(modelParams[i].nucModel, "4by4"))
						MrBayesPrint ("(%1.1lf,%1.1lf,%1.1lf,%1.1lf)\n", modelParams[i].stateFreqsDir[0], modelParams[i].stateFreqsDir[1],
							modelParams[i].stateFreqsDir[2], modelParams[i].stateFreqsDir[3]);
					else
						MrBayesPrint ("\n");
					}
				else if (modelParams[i].dataType == RESTRICTION)
					{
					MrBayesPrint ("(%1.1lf,%1.1lf)\n", modelParams[i].stateFreqsDir[0], modelParams[i].stateFreqsDir[1]);
					}
				else
					MrBayesPrint ("\n");
				}
			else if (!strcmp(modelParams[i].stateFreqPr, "Fixed"))
				{
				if (modelParams[i].dataType == DNA || modelParams[i].dataType == RNA)
					{
					if (!strcmp(modelParams[i].nucModel, "4by4"))
						MrBayesPrint ("(%1.1lf,%1.1lf,%1.1lf,%1.1lf)\n", modelParams[i].stateFreqsFix[0], modelParams[i].stateFreqsFix[1],
							modelParams[i].stateFreqsFix[2], modelParams[i].stateFreqsFix[3]);
					else
						MrBayesPrint ("\n");
					}
				else if (modelParams[i].dataType == RESTRICTION)
					{
					MrBayesPrint ("(%1.1lf,%1.1lf)\n", modelParams[i].stateFreqsFix[0], modelParams[i].stateFreqsFix[1]);
					}
				else
					MrBayesPrint ("\n");
				}


			MrBayesPrint ("   Treeheightpr     Exponential/Gamma            ");
			if (!strcmp(modelParams[i].treeHeightPr, "Exponential"))
				MrBayesPrint ("Exponential(%1.1lf)\n", modelParams[i].treeHeightExp);
			else
				MrBayesPrint ("Gamma(%1.1lf,%1.1lf)\n", modelParams[i].treeHeightGamma[0],modelParams[i].treeHeightGamma[1]);

			MrBayesPrint ("   Ratepr           Fixed/Variable=Dirichlet     %s", modelParams[i].ratePr);
			if (!strcmp(modelParams[i].ratePr, "Dirichlet"))
				MrBayesPrint ("(...,%1.1lf,...)\n", modelParams[i].ratePrDir);
			else
				MrBayesPrint ("\n");

			MrBayesPrint ("   Shapepr          Uniform/Exponential/Fixed    %s", modelParams[i].shapePr);
			if (!strcmp(modelParams[i].shapePr, "Uniform"))
				MrBayesPrint ("(%1.1lf,%1.1lf)\n", modelParams[i].shapeUni[0], modelParams[i].shapeUni[1]);
			else if (!strcmp(modelParams[i].shapePr, "Exponential"))
				MrBayesPrint ("(%1.1lf)\n", modelParams[i].shapeExp);
			else
				MrBayesPrint ("(%1.1lf)\n", modelParams[i].shapeFix);

			MrBayesPrint ("   Ratecorrpr       Uniform/Fixed                %s", modelParams[i].adGammaCorPr);
			if (!strcmp(modelParams[i].adGammaCorPr, "Uniform"))
				MrBayesPrint ("(%1.1lf,%1.1lf)\n", modelParams[i].corrUni[0], modelParams[i].corrUni[1]);
			else
				MrBayesPrint ("(%1.1lf)\n", modelParams[i].corrFix);

			MrBayesPrint ("   Pinvarpr         Uniform/Fixed                %s", modelParams[i].pInvarPr);
			if (!strcmp(modelParams[i].pInvarPr, "Uniform"))
				MrBayesPrint ("(%1.1lf,%1.1lf)\n", modelParams[i].pInvarUni[0], modelParams[i].pInvarUni[1]);
			else
				MrBayesPrint ("(%1.1lf)\n", modelParams[i].pInvarFix);

			MrBayesPrint ("   Covswitchpr      Uniform/Exponential/Fixed    %s", modelParams[i].covSwitchPr);
			if (!strcmp(modelParams[i].covSwitchPr, "Uniform"))
				MrBayesPrint ("(%1.1lf,%1.1lf)\n", modelParams[i].covswitchUni[0], modelParams[i].covswitchUni[1]);
			else if (!strcmp(modelParams[i].covSwitchPr, "Exponential"))
				MrBayesPrint ("(%1.1lf)\n", modelParams[i].covswitchExp);
			else
				MrBayesPrint ("(%1.1lf,%1.1lf)\n", modelParams[i].covswitchFix[0], modelParams[i].covswitchFix[1]);

			MrBayesPrint ("   Symdirihyperpr   Uniform/Exponential/Fixed    %s", modelParams[i].symPiPr);
			if (!strcmp(modelParams[i].symPiPr, "Uniform"))
				MrBayesPrint ("(%1.1lf,%1.1lf)\n", modelParams[i].symBetaUni[0], modelParams[i].symBetaUni[1]);
			else if (!strcmp(modelParams[i].covSwitchPr, "Exponential"))
				MrBayesPrint ("(%1.1lf)\n", modelParams[i].symBetaExp);
			else
				{
				if (modelParams[i].symBetaFix < 0)
					MrBayesPrint ("(Infinity)\n");
				else
					MrBayesPrint ("(%1.1lf)\n", modelParams[i].symBetaFix);
				}
			
			MrBayesPrint ("   Topologypr       Uniform/Constraints          %s", modelParams[i].topologyPr);
			if (!strcmp(modelParams[i].topologyPr, "Constraints"))
				{
				MrBayesPrint ("(");
				for (j=0; j<modelParams[i].numActiveConstraints; j++)
					{
					if (modelParams[i].activeConstraints[j] == YES)
						{
						if (j+1 == modelParams[i].numActiveConstraints)
							MrBayesPrint ("%d)\n", j+1);
						else
							MrBayesPrint ("%d,", j+1);
						}
					}
				}
			else
				MrBayesPrint ("\n");
				
			MrBayesPrint ("   Brlenspr         Unconstrained/Clock          %s:", modelParams[i].brlensPr);
			if (!strcmp(modelParams[i].brlensPr, "Unconstrained"))
				{
				if (!strcmp(modelParams[i].unconstrainedPr, "Uniform"))
					MrBayesPrint ("Uni(%1.1lf,%1.1lf)\n", modelParams[i].brlensUni[0], modelParams[i].brlensUni[1]);
				else
					MrBayesPrint ("Exp(%1.1lf)\n", modelParams[i].brlensExp);
				}
			else
				{
				MrBayesPrint ("%s\n", modelParams[i].clockPr);
				}
			
			MrBayesPrint ("   Speciationpr     Uniform/Exponential/Fixed    %s", modelParams[i].speciationPr);
			if (!strcmp(modelParams[i].speciationPr, "Uniform"))
				MrBayesPrint ("(%1.1lf,%1.1lf)\n", modelParams[i].speciationUni[0], modelParams[i].speciationUni[1]);
			else if (!strcmp(modelParams[i].speciationPr, "Exponential"))
				MrBayesPrint ("(%1.1lf)\n", modelParams[i].speciationExp);
			else
				MrBayesPrint ("(%1.1lf)\n", modelParams[i].speciationFix);
			
			MrBayesPrint ("   Extinctionpr     Uniform/Exponential/Fixed    %s", modelParams[i].extinctionPr);
			if (!strcmp(modelParams[i].extinctionPr, "Uniform"))
				MrBayesPrint ("(%1.1lf,%1.1lf)\n", modelParams[i].extinctionUni[0], modelParams[i].extinctionUni[1]);
			else if (!strcmp(modelParams[i].extinctionPr, "Exponential"))
				MrBayesPrint ("(%1.1lf)\n", modelParams[i].extinctionExp);
			else
				MrBayesPrint ("(%1.1lf)\n", modelParams[i].extinctionFix);
			MrBayesPrint ("   Sampleprob       <number>                     %1.2lf\n", modelParams[i].sampleProb);
			
			MrBayesPrint ("   Thetapr          Uniform/Exponential/Fixed    %s", modelParams[i].thetaPr);
			if (!strcmp(modelParams[i].thetaPr, "Uniform"))
				MrBayesPrint ("(%1.1lf,%1.1lf)\n", modelParams[i].thetaUni[0], modelParams[i].thetaUni[1]);
			else if (!strcmp(modelParams[i].thetaPr, "Exponential"))
				MrBayesPrint ("(%1.1lf)\n", modelParams[i].thetaExp);
			else
				MrBayesPrint ("(%1.1lf)\n", modelParams[i].thetaFix);

			/*
			MrBayesPrint ("   Growthpr         Uniform/Exponential/         \n");
			MrBayesPrint ("                    Fixed/Normal                 %s", modelParams[i].growthPr);
			if (!strcmp(modelParams[i].growthPr, "Uniform"))
				MrBayesPrint ("(%1.1lf,%1.1lf)\n", modelParams[i].growthUni[0], modelParams[i].growthUni[1]);
			else if (!strcmp(modelParams[i].growthPr, "Exponential"))
				MrBayesPrint ("(%1.1lf)\n", modelParams[i].growthExp);
			else if (!strcmp(modelParams[i].growthPr, "Normal"))
				MrBayesPrint ("(%1.1lf,%1.1lf)\n", modelParams[i].growthNorm[0], modelParams[i].growthNorm[1]);
			else
				MrBayesPrint ("(%1.1lf)\n", modelParams[i].growthFix); 
			*/

			MrBayesPrint ("   ------------------------------------------------------------------            \n");		
	    	MrBayesPrint ("                                                                                 \n");
			}
		}
	else if (!strcmp(helpTkn, "Ctype"))
		{
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		MrBayesPrint ("   Ctype                                                                         \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   This command sets the character ordering for standard-type data. The          \n");
	    MrBayesPrint ("   correct usage is:                                                             \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("      ctype <ordering>:<characters>                                              \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   The available options for the <ordering> specifier are:                       \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("     unordered    -- Movement directly from one state to another is              \n");
	    MrBayesPrint ("                     allowed in an instant of time.                              \n");
	    MrBayesPrint ("     ordered      -- Movement is only allowed between adjacent characters.       \n");
	    MrBayesPrint ("                     For example, perhaps only between 0 <-> 1 and 1 <-> 2       \n");
	    MrBayesPrint ("                     for a three state character ordered as 0 - 1 - 2.           \n");
	    MrBayesPrint ("     irreversible -- Rates of change for losses are 0.                           \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   The characters to which the ordering is applied is specified in manner        \n");
	    MrBayesPrint ("   that is identical to commands such as \"include\" or \"exclude\". For         \n");
	    MrBayesPrint ("   example,                                                                      \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("      ctype ordered: 10 23 45                                                    \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   defines charactes 10, 23, and 45 to be of type ordered. Similarly,            \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("      ctype irreversible: 54 - 67  71-92                                         \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   defines characters 54 to 67 and characters 71 to 92 to be of type             \n");
	    MrBayesPrint ("   irreversible. You can use the \".\" to denote the last character, and         \n");
	    MrBayesPrint ("   \"all\" to denote all of the characters. Finally, you can use the             \n");
	    MrBayesPrint ("   specifier \"\\\" to apply the ordering to every n-th character or             \n");
	    MrBayesPrint ("   you can use predefined charsets to specify the character.                     \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   Only one ordering can be used on any specific application of ctype.           \n");
	    MrBayesPrint ("   If you want to apply different orderings to different characters, then        \n");
	    MrBayesPrint ("   you need to use ctype multiple times. For example,                            \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("      ctype ordered: 1-50                                                        \n");
	    MrBayesPrint ("      ctype irreversible: 51-100                                                 \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   sets characters 1 to 50 to be ordered and characters 51 to 100 to be          \n");
	    MrBayesPrint ("   irreversible.                                                                 \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   The ctype command is only sensible with morphological (here called            \n");
	    MrBayesPrint ("   \"standard\") characters. The program ignores attempts to apply char-         \n");
	    MrBayesPrint ("   acter orderings to other types of characters, such as DNA characters.         \n");

		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		}
	else if (!strcmp(helpTkn, "Props"))
		{
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		MrBayesPrint ("   Props                                                                         \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   This command allows the user to change the details of the MCMC mechanism      \n");
	    MrBayesPrint ("   that updates the state of the chain. The useage is:                           \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("      props                                                                      \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   On typing \"props\", you will get a list of parameters to change. The         \n");
	    MrBayesPrint ("   program works as follows: On typing \"mcmc\", MrBayes figures out which       \n");
	    MrBayesPrint ("   model parameters need to be updated. For example, if you include a trans-     \n");
	    MrBayesPrint ("   ition/transversion rate parameter, then the program needs to update this      \n");
	    MrBayesPrint ("   parameter along with others, such as the tree and branch lengths. Once        \n");
	    MrBayesPrint ("   MrBayes figures out which moves are needed, it figures out the probability    \n");
	    MrBayesPrint ("   of making each move on every cycle of the chain. MrBayes updates param-       \n");
	    MrBayesPrint ("   eters in blocks; it decides which parameter to update, changes the param-     \n");
	    MrBayesPrint ("   eter, and then accepts or rejects the move according to the Metropolis-       \n");
	    MrBayesPrint ("   Hastings equation. The probability of making a move is calculated as the      \n");
	    MrBayesPrint ("   proposal rate for the move divided by the sum of the proposal rates for       \n");
	    MrBayesPrint ("   all of the other parameters that need to be updated. This command also        \n");
	    MrBayesPrint ("   allows you to change the details of each proposal mechanism. Many of the      \n");
	    MrBayesPrint ("   moves change parameters using sliding windows centered on the current         \n");
	    MrBayesPrint ("   value of the parameter. If you increase or decrease the window size, you      \n");
	    MrBayesPrint ("   will respectively decrease or increase the acceptance rate of the move.       \n");
	    MrBayesPrint ("   Some of the other moves update using a dirichlet or beta distribution,        \n");
	    MrBayesPrint ("   centered on the current values. You can change the variance parameter of      \n");
	    MrBayesPrint ("   the dirichlet or beta distribution. Finally, a few of the topology moves      \n");
	    MrBayesPrint ("   have a tuning parameter which influences the degree to which branch           \n");
	    MrBayesPrint ("   lengths are modified. If you increase this tuning parameter, you will         \n");
	    MrBayesPrint ("   make more radical changes to the branch lengths.                              \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   One word of warning: You should be extremely careful when modifying any       \n");
	    MrBayesPrint ("   of the chain parameters using \"props\". It is quite possible to completely   \n");
	    MrBayesPrint ("   wreck any hope of achieving convergence by inappropriately setting the        \n");
	    MrBayesPrint ("   chain parameters. Please exercise this command with caution.                  \n");
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		}
	else if (!strcmp(helpTkn, "Log"))
		{
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		MrBayesPrint ("   Log                                                                           \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   This command allows output to the screen to also be output to a file.         \n");
	    MrBayesPrint ("   The useage is:                                                                \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("      log start/stop filename=<name> append/replace                              \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   The options are:                                                              \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   Start/Stop     -- Starts or stops logging of output to file.                  \n");
	    MrBayesPrint ("   Append/Replace -- Either append to or replace existing file.                  \n");
	    MrBayesPrint ("   Filename       -- Name of log file (currently, the name of the log            \n");
	    MrBayesPrint ("                     file is \"%s\").\n", logFileName);
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		}
	else if (!strcmp(helpTkn, "Translate"))
		{
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		MrBayesPrint ("   Translate                                                                     \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   This command is used by MrBayes to specify the mapping between taxon names    \n");
	    MrBayesPrint ("   and taxon numbers in a Nexus tree file. For instance,						    \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("      translate                                                                  \n");
	    MrBayesPrint ("         1 Homo,                                                                 \n");
	    MrBayesPrint ("         2 Pan,                                                                  \n");
	    MrBayesPrint ("         3 Gorilla,                                                              \n");
	    MrBayesPrint ("         4 Hylobates;                                                            \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   establishes that the taxon labeled 1 in the trees that follow is Homo, the    \n");
	    MrBayesPrint ("   taxon labeled 2 is Pan, etc.                                                  \n");
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		}
	else if (!strcmp(helpTkn, "Usertree"))
		{
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		MrBayesPrint ("   Usertree                                                                      \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   This command allows you to specify a user tree. The user tree can then be     \n");
	    MrBayesPrint ("   used as a starting tree for a MCMC analysis. The format for the command is    \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("      usertree = <tree in Newick format>                                         \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   For example,                                                                  \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("      usertree = (A,B,(C,D))                                                     \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   specifies an unrooted tree of four species. Note that the program re-         \n");
	    MrBayesPrint ("   quires that trees are binary (i.e., strictly bifurcating). Hence, there       \n");
	    MrBayesPrint ("   can be only one three-way split, as shown in the example. If the tree         \n");
	    MrBayesPrint ("   is not binary, the program will return an error.                              \n");
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		}
	else if (!strcmp(helpTkn, "Mcmc"))
		{
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		MrBayesPrint ("   Mcmc                                                                          \n");
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("   This command starts the Markov chain Monte Carlo (MCMC) analysis. The         \n");
		MrBayesPrint ("   posterior probability of phylogenetic trees (and other parameters of the      \n");
		MrBayesPrint ("   substitution model) cannot be determined analytically. Instead, MCMC is       \n");
		MrBayesPrint ("   used to approximate the posterior probabilities of trees by drawing           \n");
		MrBayesPrint ("   (dependent) samples from the posterior distribution. This program can         \n");
		MrBayesPrint ("   implement a variant of MCMC called \"Metropolis-coupled Markov chain Monte    \n");
		MrBayesPrint ("   Carlo\", or MCMCMC for short. Basically, \"Nchains\" are run, with            \n");
		MrBayesPrint ("   Nchains - 1 of them heated. The chains are labelled 1, 2, ..., Nchains.       \n");
		MrBayesPrint ("   The heat that is applied to the i-th chain is B = 1 / (1 + temp X i). B       \n");
		MrBayesPrint ("   is the power to which the posterior probability is raised. When B = 0, all    \n");
		MrBayesPrint ("   trees have equal probability and the chain freely visits trees. B = 1 is      \n");
		MrBayesPrint ("   the \"cold\" chain (or the distribution of interest). MCMCMC can mix          \n");
		MrBayesPrint ("   better than ordinary MCMC; after all of the chains have gone through          \n");
		MrBayesPrint ("   one cycle, two chains are chosen at random and an attempt is made to          \n");
		MrBayesPrint ("   swap the states (with the probability of a swap being determined by the       \n");
		MrBayesPrint ("   Metropolis et al. equation). This allows the chain to potentially jump        \n");
		MrBayesPrint ("   a valley in a single bound. The correct usage is                              \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("      mcmc <parameter> = <value> ... <parameter> = <value>                       \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   For example,                                                                  \n");
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("      mcmc ngen=100000 nchains=4 temp=0.5                                        \n");
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("   performs a MCMCMC analysis with four chains with the temperature set to       \n");
		MrBayesPrint ("   0.5. The chains would be run for 100,000 cycles.                              \n");
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("   Options:                                                                      \n");
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("   Seed         -- Sets the seed number for the random number generator. The     \n");
		MrBayesPrint ("                   random number seed is initialized haphazardly at the beg-     \n");
		MrBayesPrint ("                   inning of each MrBayes session. This option allows you to     \n");
		MrBayesPrint ("                   set the seed to some specific value, thereby allowing you     \n");
		MrBayesPrint ("                   to exactly repeat an analysis. If the analysis uses swapping  \n");
		MrBayesPrint ("                   between cold and heated chains, you must also set the swap    \n");
		MrBayesPrint ("                   seed (see below) to exactly repeat the analysis.              \n");
		MrBayesPrint ("   Swapseed     -- Sets the seed used for generating the swapping sequence       \n");
		MrBayesPrint ("                   when Metropolis-coupled heated chains are used. This seed     \n");
		MrBayesPrint ("                   is initialized haphazardly at the beginning of each MrBayes   \n");
		MrBayesPrint ("                   session. This option allows you to set the seed to some       \n");
		MrBayesPrint ("                   specific value, thereby allowing you to exactly repeat a      \n");
		MrBayesPrint ("                   swap sequence. See also the 'Seed' option.                    \n");
		MrBayesPrint ("   Ngen         -- This option sets the number of cycles for the MCMC alg-       \n");
		MrBayesPrint ("                   orithm. This should be a big number as you want the chain     \n");
		MrBayesPrint ("                   to first reach stationarity, and then remain there for        \n");
		MrBayesPrint ("                   enough time to take lots of samples.                          \n");
		MrBayesPrint ("   Nruns        -- How many independent analyses are started simultaneously.     \n");
		MrBayesPrint ("   Nchains      -- How many chains are run for each analysis for the MCMCMC      \n");
		MrBayesPrint ("                   variant. The default is 4: 1 cold chain and 3 heated chains.  \n");
		MrBayesPrint ("                   If Nchains is set to 1, MrBayes will use regular MCMC sam-    \n");
		MrBayesPrint ("                   pling, without heating.                                       \n");
		MrBayesPrint ("   Temp         -- The temperature parameter for heating the chains. The higher  \n");
		MrBayesPrint ("                   the temperature, the more likely the heated chains are to     \n");
		MrBayesPrint ("                   move between isolated peaks in the posterior distribution.    \n");
		MrBayesPrint ("                   However, excessive heating may lead to very low acceptance    \n");
		MrBayesPrint ("                   rates for swaps between different chains. Before changing the \n");
		MrBayesPrint ("                   default setting, however, note that the acceptance rates of   \n");
		MrBayesPrint ("                   swaps tend to fluctuate during the burn-in phase of the run.  \n");
		MrBayesPrint ("   Reweight     -- Here, you specify three numbers, that respectively represent  \n");
		MrBayesPrint ("                   the percentage of characters to decrease in weight, the       \n");
		MrBayesPrint ("                   percentage of characters to increase in weight, and the       \n");
		MrBayesPrint ("                   increment. An increase/decrease in weight is acheived by      \n");
		MrBayesPrint ("                   replicating/removing a character in the matrix. This is       \n");
		MrBayesPrint ("                   only done to non-cold chains. The format for this parameter   \n");
		MrBayesPrint ("                   is \"reweight=(<number>,<number>)\" or \"reweight=(<number>,\"\n");
		MrBayesPrint ("                   <number>,<number>)\".                                         \n");
		MrBayesPrint ("   Swapfreq     -- This specifies how often swaps of states between chains are   \n");
		MrBayesPrint ("                   attempted. You must be running at least two chains for this   \n");
		MrBayesPrint ("                   option to be relevant. The default is Swapfreq=1, resulting   \n");
		MrBayesPrint ("                   in Nswaps (see below) swaps being tried each generation of    \n");
		MrBayesPrint ("                   the run. If Swapfreq is set to 10, then Nswaps swaps will be  \n");
		MrBayesPrint ("                   tried every tenth generation of the run.                      \n");
		MrBayesPrint ("   Nswaps       -- The number of swaps tried for each swapping generation of the \n");
		MrBayesPrint ("                   chain (see also Swapfreq).                                    \n");
		MrBayesPrint ("   Samplefreq   -- This specifies how often the Markov chain is sampled. You     \n");
		MrBayesPrint ("                   can sample the chain every cycle, but this results in very    \n");
		MrBayesPrint ("                   large output files. Thinning the chain is a way of making     \n");
		MrBayesPrint ("                   these files smaller and making the samples more independent.  \n");
		MrBayesPrint ("   Printfreq    -- This specifies how often information about the chain is       \n");
		MrBayesPrint ("                   printed to the screen.                                        \n");
		MrBayesPrint ("   Printall     -- If set to NO, only cold chains in a MCMC analysis are printed \n");
		MrBayesPrint ("                   to screen. If set to YES, both cold and heated chains will be \n");
		MrBayesPrint ("                   output. This setting only affects the printing to screen, it  \n");
		MrBayesPrint ("                   does not change the way values are written to file.           \n");
		MrBayesPrint ("   Printmax     -- The maximum number of chains to print to screen.              \n");
		MrBayesPrint ("   Mcmcdiagn    -- Determines whether acceptance ratios of moves and swaps will  \n");
		MrBayesPrint ("                   be printed to file. The file will be named similarly to the   \n");
		MrBayesPrint ("                   '.p' and '.t' files, but will have the ending '.mcmc'. If     \n");
		MrBayesPrint ("                   more than one independent analysis is run simultaneously (see \n");
		MrBayesPrint ("                   Nruns below), convergence diagnostics for tree topology will  \n");
		MrBayesPrint ("                   also be printed to this file. The convergence diagnostic used \n");
		MrBayesPrint ("                   is the average standard deviation in partition frequency      \n");
		MrBayesPrint ("                   values across independent analyses. The Burnin setting (see   \n");
		MrBayesPrint ("                   below) determines how many samples will be discarded as burnin\n");
		MrBayesPrint ("                   before calculating the partition frequencies. The Minpartfreq \n");
		MrBayesPrint ("                   setting (see below) determines the minimum partition frequency\n");
		MrBayesPrint ("                   required for a partition to be included in the calculation. As\n");
		MrBayesPrint ("                   the independent analyses approach stationarity (converge), the\n");
		MrBayesPrint ("                   value of the diagnostic is expected to approach zero.         \n");
		MrBayesPrint ("   Diagnfreq    -- The number of generations between the calculation of MCMC     \n");
		MrBayesPrint ("                   diagnostics (see Mcmcdiagn above).                            \n");
		MrBayesPrint ("   Minpartfreq  -- The minimum frequency required for a partition to be included \n");
		MrBayesPrint ("                   in the calculation of the topology convergence diagnostic. The\n");
		MrBayesPrint ("                   partition is included if the minimum frequency is reached in  \n");
		MrBayesPrint ("                   at least one of the independent tree samples that are com-    \n");
		MrBayesPrint ("                   pared.                                                        \n");
		MrBayesPrint ("   Allchains    -- If this option is set to YES, acceptance ratios for moves are \n");
		MrBayesPrint ("                   recorded for all chains, cold or heated. By default, only the \n");
		MrBayesPrint ("                   acceptance ratios for the cold chain are recorded.            \n");
		MrBayesPrint ("   Allcomps     -- If this option is set to YES, topological convergence diag-   \n");
		MrBayesPrint ("                   nostics are calculated over all pairwise comparisons of runs. \n");
		MrBayesPrint ("                   If it is set to NO, only the overall value is reported.       \n");
		MrBayesPrint ("   Relburnin    -- If this option is set to YES, then a proportion of the sampled\n");
		MrBayesPrint ("                   values will be discarded as burnin when calculating the con-  \n");
		MrBayesPrint ("                   vergence diagnostic. The proportion to be discarded is set    \n");
		MrBayesPrint ("                   with Burninfrac (see below). By default, the Relburnin option \n");
		MrBayesPrint ("                   is set to NO, resulting in a specific number of samples being \n");
		MrBayesPrint ("                   discarded instead. This number is set by Burnin (see below).  \n");
		MrBayesPrint ("   Burnin       -- Determines the number of samples (not generations) that will  \n");
		MrBayesPrint ("                   be discarded when convergence diagnostics are calculated.     \n");
		MrBayesPrint ("                   The value of this option is only relevant when Relburnin is   \n");
		MrBayesPrint ("                   set to NO.                                                    \n");
		MrBayesPrint ("   BurninFrac   -- Determines the fraction of samples that will be discarded     \n");
		MrBayesPrint ("                   when convergence diagnostics are calculated. The value of     \n");
		MrBayesPrint ("                   this option is only relevant when Relburnin is set to YES.    \n");
		MrBayesPrint ("                   Example: A value for this option of 0.25 means that 25 % of   \n");
		MrBayesPrint ("                   the samples will be discarded.                                \n");
		MrBayesPrint ("   Stoprule     -- If this option is set to NO, then the chain is run the number \n");
		MrBayesPrint ("                   of generations determined by Ngen. If it is set to YES, and   \n");
		MrBayesPrint ("                   topological convergence diagnostics are calculated (Mcmcdiagn \n");
		MrBayesPrint ("                   is set to YES), then the chain will be stopped before the pre-\n");
		MrBayesPrint ("                   determined number of generations if the convergence diagnostic\n");
		MrBayesPrint ("                   falls below the stop value.                                   \n");
		MrBayesPrint ("   Stopval      -- The critical value for the topological convergence diagnostic.\n");
		MrBayesPrint ("                   Only used when Stoprule and Mcmcdiagn are set to yes, and     \n");
		MrBayesPrint ("                   more than one analysis is run simultaneously (Nruns > 1).     \n");
		MrBayesPrint ("   Filename     -- The name of the files that will be generated. Two files       \n");
		MrBayesPrint ("                   are generated: \"<Filename>.t\" and \"<Filename>.p\".         \n");
		MrBayesPrint ("                   The .t file contains the trees whereas the .p file con-       \n");
		MrBayesPrint ("                   tains the sampled values of the parameters.                   \n");
		MrBayesPrint ("   Startingtree -- The starting tree for the chain can either be randomly        \n");
		MrBayesPrint ("                   selected or user-defined. It might be a good idea to          \n");
		MrBayesPrint ("                   start from randomly chosen trees; convergence seems           \n");
		MrBayesPrint ("                   likely if independently run chains, each of which             \n");
		MrBayesPrint ("                   started from different random trees, converge to the same     \n");
		MrBayesPrint ("                   answer.                                                       \n");
		MrBayesPrint ("   Nperts       -- This is the number of random perturbations to apply to the    \n");
		MrBayesPrint ("                   user starting tree. This allows you to have something         \n");
		MrBayesPrint ("                   between completely random and user-defined trees start        \n");
		MrBayesPrint ("                   the chain.                                                    \n");
		MrBayesPrint ("   Savebrlens   -- This specifies whether branch length information is           \n");
		MrBayesPrint ("                   saved on the trees.                                           \n");
		/*
		MrBayesPrint ("   Data         -- When Data is set to NO, the chain is run without data. This   \n");
		MrBayesPrint ("                   should be used only for debugging Mcmc proposals or for       \n");
		MrBayesPrint ("                   examining induced priors.                                     \n");
		*/
		MrBayesPrint ("   Ordertaxa    -- Determines whether taxa should be ordered before trees are    \n");
		MrBayesPrint ("                   printed to file. If set to 'Yes', terminals in the sampled    \n");
		MrBayesPrint ("                   trees will be reordered to match the order of the taxa in the \n");
		MrBayesPrint ("                   data matrix as closely as possible. By default, trees will be \n");
		MrBayesPrint ("                   printed without reordering of taxa.                           \n");
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("   Parameter       Options               Current Setting                         \n");
		MrBayesPrint ("   -----------------------------------------------------                         \n");
		MrBayesPrint ("   Seed            <number>              %ld                                     \n", chainParams.chainSeed);
		MrBayesPrint ("   Swapseed        <number>              %ld                                     \n", swapSeed);
		MrBayesPrint ("   Ngen            <number>              %d                                      \n", chainParams.numGen);
		MrBayesPrint ("   Nruns           <number>              %d                                      \n", chainParams.numRuns);
		MrBayesPrint ("   Nchains         <number>              %d                                      \n", chainParams.numChains);
		MrBayesPrint ("   Temp            <number>              %lf                                     \n", chainParams.chainTemp);
		MrBayesPrint ("   Reweight        <number>,<number>     %1.2lf v %1.2lf ^                       \n", chainParams.weightScheme[0], chainParams.weightScheme[1]);
		MrBayesPrint ("   Swapfreq        <number>              %d                                      \n", chainParams.swapFreq);
		MrBayesPrint ("   Nswaps          <number>              %d                                      \n", chainParams.numSwaps);
		MrBayesPrint ("   Samplefreq      <number>              %d                                      \n", chainParams.sampleFreq);
		MrBayesPrint ("   Printfreq       <number>              %d                                      \n", chainParams.printFreq);
		PrintYesNo (chainParams.printAll, yesNoStr);
		MrBayesPrint ("   Printall        Yes/No                %s                                      \n", yesNoStr);
		MrBayesPrint ("   Printmax        <number>              %d                                      \n", chainParams.printMax);
		PrintYesNo (chainParams.mcmcDiagn, yesNoStr);
		MrBayesPrint ("   Mcmcdiagn       Yes/No                %s                                     \n", yesNoStr);
		MrBayesPrint ("   Diagnfreq       <number>              %d                                      \n", chainParams.diagnFreq);
		MrBayesPrint ("   Minpartfreq     <number>              %1.2lf                                  \n", chainParams.minPartFreq);
		PrintYesNo (chainParams.allChains, yesNoStr);
		MrBayesPrint ("   Allchains       Yes/No                %s                                     \n", yesNoStr);
		PrintYesNo (chainParams.allComps, yesNoStr);
		MrBayesPrint ("   Allcomps        Yes/No                %s                                     \n", yesNoStr);
		PrintYesNo (chainParams.relativeBurnin, yesNoStr);
		MrBayesPrint ("   Relburnin       Yes/No                %s                                     \n", yesNoStr);
		MrBayesPrint ("   Burnin          <number>              %d                                      \n", chainParams.chainBurnIn);
		MrBayesPrint ("   Burninfrac      <number>              %1.2lf                                  \n", chainParams.burninFraction);
		PrintYesNo (chainParams.stopRule, yesNoStr);
		MrBayesPrint ("   Stoprule        Yes/No                %s                                     \n", yesNoStr);
		MrBayesPrint ("   Stopval         <number>              %1.2lf                                  \n", chainParams.stopVal);
		MrBayesPrint ("   Filename        <name>                %s.<p/t>\n", chainParams.chainFileName);
		MrBayesPrint ("   Startingtree    Random/User           %s                                      \n", chainParams.chainStartTree);
		MrBayesPrint ("   Nperts          <number>              %d                                      \n", chainParams.numStartPerts);
		PrintYesNo (chainParams.saveBrlens, yesNoStr);
		MrBayesPrint ("   Savebrlens      Yes/No                %s                                     \n", yesNoStr);
		PrintYesNo (chainParams.runWithData, yesNoStr);
		/*
		MrBayesPrint ("   Data            Yes/No                %s                                     \n", yesNoStr);
		*/
		MrBayesPrint ("   Ordertaxa       Yes/No                %s                                     \n", chainParams.orderTaxa == YES? "Yes" : "No");
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		}
	else if (!strcmp(helpTkn, "Mcmcp"))
		{
		PrintYesNo (chainParams.saveBrlens, yesNoStr);
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		MrBayesPrint ("   Mcmcp                                                                         \n");
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("   This command sets the parameters of the Markov chain Monte Carlo (MCMC)       \n");
		MrBayesPrint ("   analysis without actually starting the chain. This command is identical       \n");
		MrBayesPrint ("   in all respects to Mcmc, except that the analysis will not start after        \n");
		MrBayesPrint ("   this command is issued. For more details on the options, check the help       \n");
		MrBayesPrint ("   menu for Mcmc.\n");
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("   Parameter       Options               Current Setting                         \n");
		MrBayesPrint ("   -----------------------------------------------------                         \n");
		MrBayesPrint ("   Seed            <number>              %ld                                     \n", chainParams.chainSeed);
		MrBayesPrint ("   Swapseed        <number>              %ld                                     \n", swapSeed);
		MrBayesPrint ("   Ngen            <number>              %d                                      \n", chainParams.numGen);
		MrBayesPrint ("   Nruns           <number>              %d                                      \n", chainParams.numRuns);
		MrBayesPrint ("   Nchains         <number>              %d                                      \n", chainParams.numChains);
		MrBayesPrint ("   Temp            <number>              %lf                                     \n", chainParams.chainTemp);
		MrBayesPrint ("   Reweight        <number>,<number>     %1.2lf v %1.2lf ^                       \n", chainParams.weightScheme[0], chainParams.weightScheme[1]);
		MrBayesPrint ("   Swapfreq        <number>              %d                                      \n", chainParams.swapFreq);
		MrBayesPrint ("   Nswaps          <number>              %d                                      \n", chainParams.numSwaps);
		MrBayesPrint ("   Samplefreq      <number>              %d                                      \n", chainParams.sampleFreq);
		MrBayesPrint ("   Printfreq       <number>              %d                                      \n", chainParams.printFreq);
		PrintYesNo (chainParams.printAll, yesNoStr);
		MrBayesPrint ("   Printall        Yes/No                %s                                      \n", yesNoStr);
		MrBayesPrint ("   Printmax        <number>              %d                                      \n", chainParams.printMax);
		PrintYesNo (chainParams.mcmcDiagn, yesNoStr);
		MrBayesPrint ("   Mcmcdiagn       Yes/No                %s                                      \n", yesNoStr);
		MrBayesPrint ("   Diagnfreq       <number>              %d                                      \n", chainParams.diagnFreq);
		MrBayesPrint ("   Minpartfreq     <number>              %1.2lf                                  \n", chainParams.minPartFreq);
		PrintYesNo (chainParams.allChains, yesNoStr);
		MrBayesPrint ("   Allchains       Yes/No                %s                                      \n", yesNoStr);
		PrintYesNo (chainParams.allComps, yesNoStr);
		MrBayesPrint ("   Allcomps        Yes/No                %s                                      \n", yesNoStr);
		PrintYesNo (chainParams.relativeBurnin, yesNoStr);
		MrBayesPrint ("   Relburnin       Yes/No                %s                                      \n", yesNoStr);
		MrBayesPrint ("   Burnin          <number>              %d                                      \n", chainParams.chainBurnIn);
		MrBayesPrint ("   Burninfrac      <number>              %1.2lf                                  \n", chainParams.burninFraction);
		PrintYesNo (chainParams.stopRule, yesNoStr);
		MrBayesPrint ("   Stoprule        Yes/No                %s                                      \n", yesNoStr);
		MrBayesPrint ("   Stopval         <number>              %1.2lf                                  \n", chainParams.stopVal);
		MrBayesPrint ("   Filename        <name>                %s.<p/t>\n", chainParams.chainFileName);
		MrBayesPrint ("   Startingtree    Random/User           %s                                      \n", chainParams.chainStartTree);
		MrBayesPrint ("   Nperts          <number>              %d                                      \n", chainParams.numStartPerts);
		PrintYesNo (chainParams.saveBrlens, yesNoStr);
		MrBayesPrint ("   Savebrlens      Yes/No                %s                                      \n", yesNoStr);
		PrintYesNo (chainParams.runWithData, yesNoStr);
		/*
		MrBayesPrint ("   Data            Yes/No                %s                                      \n", yesNoStr);
		*/
		MrBayesPrint ("   Ordertaxa       Yes/No                %s                                     \n", chainParams.orderTaxa == YES? "Yes" : "No");
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		}
	else if (!strcmp(helpTkn, "Set"))
		{
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		MrBayesPrint ("   Set                                                                           \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   This command is used to set some general features of the model or program     \n");
	    MrBayesPrint ("   behavior. The correct usage is                                                \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("      set <parameter>=<value> ... <parameter>=<value>                            \n");
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("   Only five parameters can be changed using \"set\". First, you can set         \n");
		MrBayesPrint ("   the autoclose feature:                                                        \n");
		MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("      set autoclose=<yes/no>                                                     \n");
		MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("   If autoclose is set to yes, then the program will not prompt you during       \n");
		MrBayesPrint ("   the course of executing a file. Second, you can set the partition that        \n");
		MrBayesPrint ("   is in effect:                                                                 \n");
		MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("      set partition=<partition id>                                               \n");
		MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("   A valid partition ID is either a number or a partition name. This command     \n");
		MrBayesPrint ("   enforces use of a specific partitioning of the data. When the program         \n");
		MrBayesPrint ("   executes, a default partition (that may not divide the data at all) is        \n");
		MrBayesPrint ("   created called \"Default\". You can always go back to the original or         \n");
		MrBayesPrint ("   default partition by typing                                                   \n");
		MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("      set partition=default                                                      \n");
		MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("   or                                                                            \n");
		MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("      set partition=1                                                            \n");
		MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("   Third, you can set the nowarnings feature:                                    \n");
		MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("      set nowarnings=<yes/no>                                                    \n");
		MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("   If nowarnings is set to yes, then the program will not prompt you when        \n");
		MrBayesPrint ("   overwriting or appending an output file that is already present. If           \n");
		MrBayesPrint ("   nowarnings=no (the default setting) then the program prompts the user         \n");
		MrBayesPrint ("   before overwriting output files.                                              \n");
		MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("   Fourth, you can set the quitonerror feature:                                  \n");
		MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("      set quitonerror=<yes/no>                                                   \n");
		MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("   If quitonerror is set to yes, then the program will quit when an error is     \n");
		MrBayesPrint ("   encountered, after printing an error message. If quitonerror=no (the default  \n");
		MrBayesPrint ("   setting) then the program will wait for additional commands from the command  \n");
		MrBayesPrint ("   line after printing the error message.                                        \n");
		MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("   Fifth, you can set the autooverwrite feature:                                 \n");
		MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("      set autooverwrite=<yes/no>                                                 \n");
		MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("   When nowarnings is set to yes, by default MrBayes will overwrite output       \n");
		MrBayesPrint ("   files. If autooverwrite=no, output will be appendend if the output file       \n");
		MrBayesPrint ("   already exists. The default is autooverwrite=yes.                             \n"); 
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		}
	else if (!strcmp(helpTkn, "Charset"))
		{
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		MrBayesPrint ("   Charset                                                                       \n");
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("   This command defines a character set. The format for the charset command      \n"); 
		MrBayesPrint ("   is                                                                            \n"); 
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("      charset <name> = <character numbers>                                       \n"); 
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("   For example, \"charset first_pos = 1-720\\3\" defines a character set         \n");
		MrBayesPrint ("   called \"first_pos\" that includes every third site from 1 to 720.            \n");
		MrBayesPrint ("   The character set name cannot have any spaces in it. The slash (\\)           \n");
		MrBayesPrint ("   is a nifty way of telling the program to assign every third (or               \n");
		MrBayesPrint ("   second, or fifth, or whatever) character to the character set.                \n");
		MrBayesPrint ("   This option is best used not from the command line, but rather as a           \n");
		MrBayesPrint ("   line in the mrbayes block of a file. Note that you can use \".\" to           \n");
		MrBayesPrint ("   stand in for the last character (e.g., charset 1-.\\3).                       \n");
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		}
	else if (!strcmp(helpTkn, "Outgroup"))
		{
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		MrBayesPrint ("   Outgroup                                                                      \n");
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("   This command assigns a taxon to the outgroup. The correct usage is:           \n"); 
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("      outgroup <number>/<taxon name>                                             \n"); 
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("   For example, \"outgroup 3\" assigns the third taxon in the matrix to be       \n");
		MrBayesPrint ("   the outgroup. Similarly, \"outgroup Homo_sapiens\" assings the taxon          \n");
		MrBayesPrint ("   \"Homo_sapiens\" to be the outgroup (assuming that there is a taxon named     \n");
		MrBayesPrint ("   \"Homo_sapiens\" in the matrix). Only a single taxon can be assigned to       \n");
		MrBayesPrint ("   be the outgroup.                                                              \n");
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		}
	else if (!strcmp(helpTkn, "Showtree"))
		{
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		MrBayesPrint ("   Showtree                                                                      \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   This command shows the current user tree. The correct usage                   \n");
	    MrBayesPrint ("   is \"showtree\".                                                              \n");
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		}
	else if (!strcmp(helpTkn, "Deroot"))
		{
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		MrBayesPrint ("   Deroot                                                                        \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   This command deroots the user tree. If the tree is already unrooted, a        \n");
	    MrBayesPrint ("   warning is issued. The correct usage is \"deroot\".                           \n");
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		}
	else if (!strcmp(helpTkn, "Root"))
		{
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		MrBayesPrint ("   Root                                                                          \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   This command roots the tree. If the tree is already rooted, a warning         \n");
	    MrBayesPrint ("   is issued. The tree is rooted at the midpoint between the outgroup species    \n");
		MrBayesPrint ("   and the ingroup species. The correct usage is \"root\".                       \n");
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		}
	else if (!strcmp(helpTkn, "Taxset"))
		{
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		MrBayesPrint ("   Taxset                                                                        \n");
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("   This command defines a taxon set. The format for the taxset command           \n"); 
		MrBayesPrint ("   is                                                                            \n"); 
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("      taxset <name> = <taxon names or numbers>                                   \n"); 
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("   For example, \"taxset apes = Homo Pan Gorilla Orang gibbon\" defines a        \n");
		MrBayesPrint ("   taxon set called \"apes\" that includes five taxa (namely, apes).             \n");
		MrBayesPrint ("   You can assign up to 30 taxon sets. This option is best used                  \n");
		MrBayesPrint ("   not from the command line but rather as a line in the mrbayes block           \n");
		MrBayesPrint ("   of a file.                                                                    \n");
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		}
	else if (!strcmp(helpTkn, "Charstat"))
		{
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		MrBayesPrint ("   Charstat                                                                      \n");
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("   This command shows the status of all the characters. The correct usage        \n");
		MrBayesPrint ("   is                                                                            \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("      charstat                                                                   \n");
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("   After typing \"charstat\", the character number, whether it is excluded       \n");
		MrBayesPrint ("   or included, and the partition identity are shown. The output is paused       \n");
		MrBayesPrint ("   every 100 characters. This pause can be turned off by setting autoclose       \n");
		MrBayesPrint ("   to \"yes\" (set autoclose=yes).                                               \n");
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		}
	else if (!strcmp(helpTkn, "Taxastat"))
		{
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		MrBayesPrint ("   Taxastat                                                                      \n");
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("   This command shows the status of all the taxa. The correct usage is           \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("      taxastat                                                                   \n");
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("   After typing \"taxastat\", the taxon number, name, and whether it is          \n");
		MrBayesPrint ("   excluded or included are shown.                                               \n");
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		}
	else if (!strcmp(helpTkn, "Partition"))
		{
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		MrBayesPrint ("   Partition                                                                     \n");
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("   This command allows you to specify a character partition. The format for      \n"); 
		MrBayesPrint ("   this command is                                                               \n"); 
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("      partition <name> = <num parts>:<chars in first>, ...,<chars in last>       \n"); 
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("   For example, \"partition by_codon = 3:1st_pos,2nd_pos,3rd_pos\" specifies     \n"); 
		MrBayesPrint ("   a partition called \"by_codon\" which consists of three parts (first,         \n"); 
		MrBayesPrint ("   second, and third codon positions). Here, we are assuming that the sites      \n"); 
		MrBayesPrint ("   in each partition were defined using the charset command. You can specify     \n"); 
		MrBayesPrint ("   a partition without using charset as follows:                                 \n"); 
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("      partition by_codon = 3:1 4 6 9 12,2 5 7 10 13,3 6 8 11 14                  \n"); 
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("   However, we recommend that you use the charsets to define a set of char-      \n"); 
		MrBayesPrint ("   acters and then use these predefined sets when defining the partition.        \n"); 
		MrBayesPrint ("   Also, it makes more sense to define a partition as a line in the mrbayes      \n"); 
		MrBayesPrint ("   block than to issue the command from the command line (then again, you        \n"); 
		MrBayesPrint ("   may be a masochist, and want to do extra work).                               \n"); 
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		}
	else if (!strcmp(helpTkn, "Exclude"))
		{
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		MrBayesPrint ("   Exclude                                                                       \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   This command excludes characters from the analysis. The correct usage is      \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("      exclude <number> <number> <number>                                         \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   or                                                                            \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("      exclude <number> - <number>                                                \n");
	    MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   or                                                                            \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("      exclude <charset>                                                          \n");
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("   or some combination thereof. Moreover, you can use the specifier \"\\\" to    \n");
		MrBayesPrint ("   exclude every nth character. For example, the following                       \n");
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("      exclude 1-100\\3                                                           \n");
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("   would exclude every third character. As a specific example,                   \n");
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("      exclude 2 3 10-14 22                                                       \n");
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("   excludes sites 2, 3, 10, 11, 12, 13, 14, and 22 from the analysis. Also,      \n");
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("      exclude all                                                                \n");
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("   excludes all of the characters from the analysis. Excluding all characters    \n");
		MrBayesPrint ("   does not leave you much information for inferring phylogeny.                  \n");
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		}
	else if (!strcmp(helpTkn, "Include"))
		{
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		MrBayesPrint ("   Include                                                                       \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   This command includes characters that were previously excluded from the       \n");
	    MrBayesPrint ("   analysis. The correct usage is                                                \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("      include <number> <number> <number>                                         \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   or                                                                            \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("      include <number> - <number>                                                \n");
	    MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   or                                                                            \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("      include <charset>                                                          \n");
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("   or some combination thereof. Moreover, you can use the specifier \"\\\" to    \n");
		MrBayesPrint ("   include every nth character. For example, the following                       \n");
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("      include 1-100\\3                                                           \n");
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("   would include every third character. As a specific example,                   \n");
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("      include 2 3 10-14 22                                                       \n");
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("   includes sites 2, 3, 10, 11, 12, 13, 14, and 22 from the analysis. Also,      \n");
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("      include all                                                                \n");
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("   includes all of the characters in the analysis. Including all of the          \n");
		MrBayesPrint ("   characters (even if many of them are bad) is a very total-evidence-like       \n");
		MrBayesPrint ("   thing to do. Doing this will make a certain group of people very happy.       \n");
		MrBayesPrint ("   On the other hand, simply using this program would make those same people     \n");
		MrBayesPrint ("   unhappy.                                                                      \n");
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		}
	else if (!strcmp(helpTkn, "Delete"))
		{
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		MrBayesPrint ("   Delete                                                                        \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   This command deletes taxa from the analysis. The correct usage is:            \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("      delete <name and/or number and/or taxset> ...                              \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   A list of the taxon names or taxon numbers (labelled 1 to ntax in the order   \n");
	    MrBayesPrint ("   in the matrix) or taxset(s) can be used.  For example, the following:         \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("      delete 1 2 Homo_sapiens                                                    \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   deletes taxa 1, 2, and the taxon labelled Homo_sapiens from the analysis.     \n");
	    MrBayesPrint ("   You can also use \"all\" to delete all of the taxa. For example,              \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("      delete all                                                                 \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   deletes all of the taxa from the analysis. Of course, a phylogenetic anal-    \n");
	    MrBayesPrint ("   ysis that does not include any taxa is fairly uninteresting.                  \n");
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		}
	else if (!strcmp(helpTkn, "Restore"))
		{
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		MrBayesPrint ("   Restore                                                                       \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   This command restores taxa to the analysis. The correct usage is:             \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("      restore <name and/or number and/or taxset> ...                             \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   A list of the taxon names or taxon numbers (labelled 1 to ntax in the order   \n");
	    MrBayesPrint ("   in the matrix) or taxset(s) can be used.  For example, the following:         \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("      restore 1 2 Homo_sapiens                                                   \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   restores taxa 1, 2, and the taxon labelled Homo_sapiens to the analysis.      \n");
	    MrBayesPrint ("   You can also use \"all\" to restore all of the taxa. For example,             \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("      restore all                                                                \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   restores all of the taxa to the analysis.                                     \n");
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		}
	else if (!strcmp(helpTkn, "Quit"))
		{
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		MrBayesPrint ("   Quit                                                                          \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   This command quits the program. The correct usage is:                         \n");
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("      quit                                                                       \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   It is a very easy command to use properly.                                    \n");
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		}
	else if (!strcmp(helpTkn, "Disclaimer"))
		{
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		MrBayesPrint ("   Disclaimer                                                                    \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   This command shows the disclaimer for the program. In short, the disclaimer   \n");
	    MrBayesPrint ("   states that the authors (John Huelsenbeck and Fredrik Ronquist) are not       \n");
	    MrBayesPrint ("   responsible for any silly things you may do to your computer or any           \n");
	    MrBayesPrint ("   unforseen but possibly nasty things the computer program may inadvertently    \n");
	    MrBayesPrint ("   do to you.                                                                    \n");
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		}
	else if (!strcmp(helpTkn, "Unlink"))
		{
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		MrBayesPrint ("   Unlink                                                                         \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   This command unlinks model parameters across partitions of the data. The      \n");
	    MrBayesPrint ("   correct usage is:                                                             \n");
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("      unlink <parameter name> = (<all> or <partition list>)                      \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   A little background is necessary to understand this command. Upon exe-        \n");
	    MrBayesPrint ("   cution of a file, a default partition is set up. This partition refer-        \n");
	    MrBayesPrint ("   enced either by its name (\"default\") or number (0). If your data are        \n");
	    MrBayesPrint ("   all of one type, then this default partition does not actually divide up      \n");
	    MrBayesPrint ("   your characters. However, if your datatype is mixed, then the default         \n");
	    MrBayesPrint ("   partition contains as many divisions as there are datatypes in your           \n");
	    MrBayesPrint ("   character matrix. Of course, you can also define other partitions, and        \n");
	    MrBayesPrint ("   switch among them using the set command (\"set partition=<name/number>\").    \n");
	    MrBayesPrint ("   Importantly, you can also assign model parameters to individual part-         \n");
	    MrBayesPrint ("   itions or to groups of them using the \"applyto\" option in lset and          \n");
	    MrBayesPrint ("   prset. When the program attempts to perform an analysis, the model is         \n");
	    MrBayesPrint ("   set for individual partitions. If the same parameter applies to differ-       \n");
	    MrBayesPrint ("   partitions and if that parameter has the same prior, then the program         \n");
	    MrBayesPrint ("   will link the parameters: that is, it will use a single value for the         \n");
	    MrBayesPrint ("   parameter. The program's default, then, is to strive for parsimony.           \n");
	    MrBayesPrint ("   However, there are lots of cases where you may want unlink a parameter        \n");
	    MrBayesPrint ("   across partitions. For example, you may want a different transition/          \n");
	    MrBayesPrint ("   transversion rate ratio to apply to different partitions. This command        \n");
	    MrBayesPrint ("   allows you to unlink the parameters, or to make them different across         \n");
	    MrBayesPrint ("   partitions. The converse of this command is \"link\", which links to-         \n");
	    MrBayesPrint ("   gether parameters that were previously told to be different. The list         \n");
	    MrBayesPrint ("   of parameters that can be unlinked includes:                                  \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("      Tratio          -- Transition/transversion rate ratio                      \n");
	    MrBayesPrint ("      Revmat          -- Substitution rates of GTR model                         \n");
	    MrBayesPrint ("      Omega           -- Nonsynonymous/synonymous rate ratio                     \n");
	    MrBayesPrint ("      Statefreq       -- Character state frequencies                             \n");
	    MrBayesPrint ("      Shape           -- Gamma shape parameter                                   \n");
	    MrBayesPrint ("      Pinvar          -- Proportion of invariable sites                          \n");
	    MrBayesPrint ("      Correlation     -- Correlation parameter of autodiscrete gamma             \n");
	    MrBayesPrint ("      Switchrates     -- Switching rates for covarion model                      \n");
	    MrBayesPrint ("      Brlens          -- Branch lengths of tree                                  \n");
	    MrBayesPrint ("      Topology        -- Topology of tree                                        \n");
	    MrBayesPrint ("      Speciationrates -- Speciation rates for birth-death process                \n");
	    MrBayesPrint ("      Extinctionrates -- Extinction rates for birth-death process                \n");
	    MrBayesPrint ("      Theta           -- Parameter for coalescence process                       \n");
		MrBayesPrint ("      Growthrate      -- Growth rate of coalescence process                      \n"); 
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   For example,                                                                  \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("      unlink shape=(all)                                                         \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   unlinks the gamma shape parameter across all partitions of the data.          \n");
	    MrBayesPrint ("   You can use \"showmodel\" to see the current linking status of the            \n");
	    MrBayesPrint ("   characters.                                                                   \n");
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		}
	else if (!strcmp(helpTkn, "Link"))
		{
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		MrBayesPrint ("   Link                                                                          \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   This command links model parameters across partitions of the data. The        \n");
	    MrBayesPrint ("   correct usage is:                                                             \n");
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("      link <parameter name> = (<all> or <partition list>)                        \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   The list of parameters that can be linked includes:                           \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("      Tratio          -- Transition/transversion rate ratio                      \n");
	    MrBayesPrint ("      Revmat          -- Substitution rates of GTR model                         \n");
	    MrBayesPrint ("      Omega           -- Nonsynonymous/synonymous rate ratio                     \n");
	    MrBayesPrint ("      Statefreq       -- Character state frequencies                             \n");
	    MrBayesPrint ("      Shape           -- Gamma shape parameter                                   \n");
	    MrBayesPrint ("      Pinvar          -- Proportion of invariable sites                          \n");
	    MrBayesPrint ("      Correlation     -- Correlation parameter of autodiscrete gamma             \n");
	    MrBayesPrint ("      Switchrates     -- Switching rates for covarion model                      \n");
	    MrBayesPrint ("      Brlens          -- Branch lengths of tree                                  \n");
	    MrBayesPrint ("      Topology        -- Topology of tree                                        \n");
	    MrBayesPrint ("      Speciationrates -- Speciation rates for birth-death process                \n");
	    MrBayesPrint ("      Extinctionrates -- Extinction rates for birth-death process                \n");
	    MrBayesPrint ("      Theta           -- Parameter for coalescence process                       \n");
	    MrBayesPrint ("      Growthrate      -- Growth rate of coalescence process                      \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   For example,                                                                  \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("      link shape=(all)                                                           \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   links the gamma shape parameter across all partitions of the data.            \n");
	    MrBayesPrint ("   You can use \"showmodel\" to see the current linking status of the            \n");
	    MrBayesPrint ("   characters. For more information on this command, see the help menu           \n");
	    MrBayesPrint ("   for link's converse, unlink (\"help unlink\");                                \n");
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		}
	else if (!strcmp(helpTkn, "Help"))
		{
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		MrBayesPrint ("   Help                                                                          \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   This command provides useful information on the use of this program. The      \n");
	    MrBayesPrint ("   correct usage is                                                              \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("      help                                                                       \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   which gives a list of all available commands with a brief description of      \n");
	    MrBayesPrint ("   each or                                                                       \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("      help <command>                                                             \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   which gives detailed information on the use of <command>.                     \n");
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		}
	else if (!strcmp(helpTkn, "Sump"))
		{
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		MrBayesPrint ("   Sump                                                                          \n");
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("   During a MCMC analysis, MrBayes prints the sampled parameter values to a tab- \n");
		MrBayesPrint ("   delimited text file. This file has the extension \".p\". The command 'Sump'   \n");
		MrBayesPrint ("   summarizes the information in the parameter file. By default, the name of the \n");
		MrBayesPrint ("   parameter file is assumed to be the name of the last matrix-containing nexus  \n");
		MrBayesPrint ("   file, but with a '.p' extension. You can set 'Sump' to summarize the infor-   \n");
		MrBayesPrint ("   mation in any other parameter file by setting the 'filename' option to the    \n");
		MrBayesPrint ("   appropriate file name. The 'Sump' command does not require a matrix to be     \n");
		MrBayesPrint ("   read in first. When you invoke the 'Sump' command, three items are output:    \n");
		MrBayesPrint ("   (1) a generation plot of the likelihood values; (2) estimates of the mar-     \n");
		MrBayesPrint ("   ginal likelihood of the model; and (3) a table with the mean, variance, and   \n");
		MrBayesPrint ("   95 percent credible interval for the sampled parameters. Each of these items  \n");
		MrBayesPrint ("   can be switched on or off using the options 'Plot', 'Marglike', and 'Table'.  \n");
		MrBayesPrint ("   By default, all three items are output but only to the screen. If output to   \n");
		MrBayesPrint ("   a file is also desired, set 'Printtofile' to 'Yes'. The name of the output    \n");
		MrBayesPrint ("   file is specified by setting the 'Outputname' option. When a new matrix is    \n");
		MrBayesPrint ("   read in or when the 'Mcmc' output filename or 'Sump' input filename is        \n");
		MrBayesPrint ("   changed, the 'Sump' outputname is changed as well. If you want to output to   \n");
		MrBayesPrint ("   another file than the default, make sure you specify the outputname every     \n");
		MrBayesPrint ("   time you invoke 'Sump'. If the specified outputfile already exists, you will  \n");
		MrBayesPrint ("   be prompted about whether you like to overwrite it or append to it. This      \n");
		MrBayesPrint ("   behavior can be altered using 'Set nowarn=yes'; see the help for the 'Set'    \n");
		MrBayesPrint ("   command. When running 'Sump' you typically want to discard a specified        \n");
		MrBayesPrint ("   number of samples from the beginning of the chain as the burn in. Note that   \n");
		MrBayesPrint ("   the 'Burnin' value of the 'Sump' command is set separately from the 'Burnin'  \n");
		MrBayesPrint ("   values of the 'Sumt' and 'Mcmc' commands. That is, if you issue               \n");
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("   sump burnin = 4000                                                            \n");
		MrBayesPrint ("   sumt burnin = 2000                                                            \n");
		MrBayesPrint ("   sump                                                                          \n");
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("   the burnin of the last 'Sump' command is 4000 and not 2000. The burnin        \n");
		MrBayesPrint ("   values are reset to 0 every time a new matrix is read in. Similarly, 'Plot',  \n");
		MrBayesPrint ("   'Marglike' and 'Table' are all set to 'Yes' and 'Printtofile' to 'No' (the    \n");
		MrBayesPrint ("   default values) when a new matrix is processed. If you have run several       \n");
		MrBayesPrint ("   independent MCMC analyses, you may want to summarize and compare the samples  \n");
		MrBayesPrint ("   from each of these runs. To do this, set 'Nruns' to the number of runs you    \n");
		MrBayesPrint ("   want to compare and make sure that the '.p' files are named using the MrBayes \n");
		MrBayesPrint ("   convention (<filename>.run1.p, <filename>.run2.p, etc). When you run several  \n");
		MrBayesPrint ("   independent analyses simultaneously in MrBayes, the 'Nruns' and 'Filename'    \n");
		MrBayesPrint ("   options are automatically set such that 'Sump' will summarize all the resul-  \n");
		MrBayesPrint ("   ting output files.                                                            \n");
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("   Options:                                                                      \n");
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("   Burnin       -- Determines the number of samples that will be discarded from  \n");
		MrBayesPrint ("                   the input file before calculating summary statistics. If there\n");
		MrBayesPrint ("                   are several input files, the same number of samples will be   \n");
		MrBayesPrint ("                   discarded from each. Note that the burnin is set separately   \n");
		MrBayesPrint ("                   for the 'sump', 'sumt', and 'mcmc' commands.                  \n");
		MrBayesPrint ("   Nruns        -- Determines how many '.p' files from independent analyses that \n");
		MrBayesPrint ("                   will be summarized. If Nruns > 1 then the names of the files  \n");
		MrBayesPrint ("                   are derived from 'Filename' by adding '.run1.p', '.run2.p',   \n");
		MrBayesPrint ("                   etc. If Nruns=1, then the single file name is obtained by     \n");
		MrBayesPrint ("                   adding '.p' to 'Filename'.                                    \n");
		MrBayesPrint ("   Filename     -- The name of the file to be summarized. This is the base of the\n");
		MrBayesPrint ("                   file name to which endings are added according to the current \n");
		MrBayesPrint ("                   setting of the 'Nruns' parameter. If 'Nruns' is 1, then only  \n");
		MrBayesPrint ("                   '.p' is added to the file name. Otherwise, the endings will   \n");
		MrBayesPrint ("                   be '.run1.p', '.run2.p', etc.                                 \n");
		MrBayesPrint ("   Printtofile  -- Determines whether results will be printed to file.           \n");
		MrBayesPrint ("   Outputname   -- Name of the file to which 'sump' results will be printed if   \n");
		MrBayesPrint ("                   'Printtofile' is set to YES.                                  \n");
		MrBayesPrint ("   Plot         -- Determines whether a likelihood plot should be output.        \n");
		MrBayesPrint ("   Marglike     -- Determines whether estimates of marginal model likelihoods    \n");
		MrBayesPrint ("                   should be calculated. The marginal model likelihoods are use- \n");
		MrBayesPrint ("                   ful in Bayesian model testing.                                \n");
		MrBayesPrint ("   Table        -- Determines whether the table summarizing the parameter value  \n");
		MrBayesPrint ("                   samples should be output.                                     \n");
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("   Current settings:                                                             \n");
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("   Parameter       Options                  Current Setting                      \n");
		MrBayesPrint ("   --------------------------------------------------------                      \n");
		MrBayesPrint ("   Burnin          <number>                 %d                                   \n", sumpParams.sumpBurnIn);
		MrBayesPrint ("   Nruns           <number>                 %d                                   \n", sumpParams.numRuns);
		if (sumpParams.numRuns == 1)
			MrBayesPrint ("   Filename        <name>                   %s<.p>\n", sumpParams.sumpFileName);
		else
			MrBayesPrint ("   Filename        <name>                   %s<.run<i>.p>\n", sumpParams.sumpFileName);
		MrBayesPrint ("   Printtofile     Yes/No                   %s                                   \n", sumpParams.printToFile == YES ? "Yes" : "No");
		MrBayesPrint ("   Outputname      <name>                   %s\n", sumpParams.sumpOutfile);
		MrBayesPrint ("   Plot            Yes/No                   %s                                   \n", sumpParams.plot == YES ? "Yes" : "No");
		MrBayesPrint ("   Marglike        Yes/No                   %s                                   \n", sumpParams.margLike == YES ? "Yes" : "No");
		MrBayesPrint ("   Table           Yes/No                   %s                                   \n", sumpParams.table == YES ? "Yes" : "No");
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		}
	else if (!strcmp(helpTkn, "Comparetree"))
		{
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		MrBayesPrint ("   Comparetree                                                                   \n");
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("   This command compares the trees in two files, called \"filename1\" and        \n");
		MrBayesPrint ("   \"filename2\". It will output a bivariate plot of the split frequencies       \n");
		MrBayesPrint ("   as well as plots of the tree distance as a function of the generation. The    \n");
		MrBayesPrint ("   plots can be used to get a quick indication of whether two runs have con-     \n");
		MrBayesPrint ("   verged onto the same set of trees. The \"Comparetree\" command will also      \n");
		MrBayesPrint ("   produce a \".parts\" file and a \".dists\" file (these file endings are added \n");
		MrBayesPrint ("   to the end of the \"Outputname\"). The \".parts\" file contains the paired    \n");
		MrBayesPrint ("   split frequencies from the two tree samples; the \".dists\" file contains the \n");
		MrBayesPrint ("   tree distance values. Note that the \"Sumt\" command provides a different     \n");
		MrBayesPrint ("   set of convergence diagnostics tools that you may also want to explore. Un-   \n");
		MrBayesPrint ("   like \"Comparetree\", \"Sumt\" can compare more than two tree samples and     \n");
		MrBayesPrint ("   will calculate consensus trees and split frequencies from the pooled samples. \n");
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("   Parameter       Options                  Current Setting                      \n");
		MrBayesPrint ("   --------------------------------------------------------                      \n");
		MrBayesPrint ("   Filename1       <name>                   %s                                   \n", comptreeParams.comptFileName1);
		MrBayesPrint ("   Filename2       <name>                   %s                                   \n", comptreeParams.comptFileName2);
		MrBayesPrint ("   Outputname      <name>                   %s                                   \n", comptreeParams.comptOutfile);
		MrBayesPrint ("   Burnin          <number>                 %d                                   \n", comptreeParams.comptBurnIn);
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		}
	else if (!strcmp(helpTkn, "Sumt"))
		{
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		MrBayesPrint ("   Sumt                                                                          \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   This command is used to produce summary statistics for trees sampled during   \n");
	    MrBayesPrint ("   a Bayesian MCMC analysis. You can either summarize trees from one individual  \n");
	    MrBayesPrint ("   analysis, or trees coming from several independent analyses. In either case,  \n");
	    MrBayesPrint ("   all the sampled trees are read in and the proportion of the time any single   \n");
	    MrBayesPrint ("   taxon bipartition is found is counted. The proportion of the time that the    \n");
	    MrBayesPrint ("   bipartition is found is an approximation of the posterior probability of the  \n");
	    MrBayesPrint ("   bipartition. (Remember that a taxon bipartition is defined by removing a      \n");
	    MrBayesPrint ("   branch on the tree, dividing the tree into those taxa to the left and right   \n");
	    MrBayesPrint ("   of the removed branch. This set is called a taxon bipartition.) The branch    \n");
	    MrBayesPrint ("   length of the bipartition is also recorded, if branch lengths have been saved \n");
	    MrBayesPrint ("   to file. The result is a list of the taxon bipartitions found, the frequency  \n");
	    MrBayesPrint ("   with which they were found, the posterior probability of the bipartition      \n");
	    MrBayesPrint ("   and, if the branch lengths were recorded, the mean and variance of the length \n");
	    MrBayesPrint ("   of the branch.                                                                \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   The partition information is output to a file with the suffix \".parts\".     \n");
	    MrBayesPrint ("   A consensus tree is also printed to a file with the suffix \".con\" and       \n");
	    MrBayesPrint ("   printed to the screen as a cladogram, and as a phylogram if branch lengths    \n");
	    MrBayesPrint ("   have been saved. The consensus tree is either a 50 percent majority rule tree \n");
	    MrBayesPrint ("   or a majority rule tree showing all compatible partitions. If branch lengths  \n");
	    MrBayesPrint ("   have been recorded during the run, the \".con\" file will contain a consensus \n");
	    MrBayesPrint ("   tree with branch lengths and interior nodes labelled with support values.     \n");
	    MrBayesPrint ("   This tree can be viewed in a program such as TreeView.                        \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   Finally, MrBayes produces a file with the ending \".trprobs\" that contains a \n");
	    MrBayesPrint ("   list of all the trees that were found during the MCMC analysis, sorted by     \n");
	    MrBayesPrint ("   their probabilities. This list of trees can be used to construct a credible   \n");
	    MrBayesPrint ("   set of trees. For example, if you want to construct a 95 percent credible set \n");
	    MrBayesPrint ("   of trees, you include all of those trees whose cumulated probability is less  \n");
	    MrBayesPrint ("   than or equal to 0.95. You have the option of displaying the trees to the     \n");
	    MrBayesPrint ("   screen using the \"Showtreeprobs\" option. The default is to not display the  \n");
	    MrBayesPrint ("   trees to the screen; the number of different trees sampled by the chain can   \n");
	    MrBayesPrint ("   be quite large. If you are analyzing a large set of taxa, you may actually    \n");
	    MrBayesPrint ("   want to skip the calculation of tree probabilities entirely by setting        \n");
	    MrBayesPrint ("   \"Calctreeprobs\" to NO.                                                      \n");
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("   When calculating summary statistics you probably want to skip those trees that\n");
	    MrBayesPrint ("   were sampled in the initial part of the run, the so-called burn-in period. The\n");
	    MrBayesPrint ("   number of skipped samples is controlled by the \"burnin\" setting. The default\n");
	    MrBayesPrint ("   is 0 but you typically want to override this setting.                         \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   If you are summarizing the trees sampled in several independent analyses,     \n");
	    MrBayesPrint ("   such as those resulting from setting the \"Nruns\" option of the \"Mcmc\" com-\n");
	    MrBayesPrint ("   mand to a value larger than 1, MrBayes will also calculate convergence diag-  \n");
	    MrBayesPrint ("   nostics for the sampled topologies and branch lengths. These values can help  \n");
	    MrBayesPrint ("   you determine whether it is likely that your chains have converged.           \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   The \"Sumt\" command expands the \"Filename\" according to the current values \n");
	    MrBayesPrint ("   of the \"Nruns\" and \"Ntrees\" options. For instance, if both \"Nruns\" and  \n");
	    MrBayesPrint ("   \"Ntrees\" are set to 1, \"Sumt\" will try to open a file named               \n");
	    MrBayesPrint ("   \"<Filename>.t\". If \"Nruns\" is set to 2 and \"Ntrees\" to 1, then \"Sumt\" \n");
	    MrBayesPrint ("   will open two files, \"<Filename>.run1.t\" and \"<Filename>.run2.t\", etc.    \n");
	    MrBayesPrint ("   By default, the \"Filename\" option will be set such that \"Sumt\" auto-      \n");
	    MrBayesPrint ("   matically summarizes all the results from your immediately preceding \"Mcmc\" \n");
	    MrBayesPrint ("   command. You can also use the \"Sumt\" command to summarize tree samples in   \n");
	    MrBayesPrint ("   older analyses. If you want to do that, remember to first read in a matrix    \n");
	    MrBayesPrint ("   so that MrBayes knows what taxon names to expect in the trees. Then set the   \n");
	    MrBayesPrint ("   \"Nruns\", \"Ntrees\" and \"Filename\" options appropriately.                 \n");
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("   Options:                                                                      \n");
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("   Burnin        -- Determines the number of samples that will be discarded from \n");
		MrBayesPrint ("                    the input file before calculating summary statistics. If     \n");
		MrBayesPrint ("                    there are several input files, the same number of samples    \n");
		MrBayesPrint ("                    will be discarded from each. Note that the burnin is set     \n");
		MrBayesPrint ("                    separately for the 'sumt', 'sump', and 'mcmc' commands.      \n");
		MrBayesPrint ("   Nruns         -- Determines how many '.t' files from independent analyses that\n");
		MrBayesPrint ("                    will be summarized. If Nruns > 1 then the names of the files \n");
		MrBayesPrint ("                    are derived from 'Filename' by adding '.run1.t', '.run2.t',  \n");
		MrBayesPrint ("                    etc. If Nruns=1 and Ntrees=1 (see below), then only '.t' is  \n");
		MrBayesPrint ("                    added to 'Filename'.                                         \n");
		MrBayesPrint ("   Ntrees        -- Determines how many trees there are in the sampled model. If \n");
		MrBayesPrint ("                    'Ntrees' > 1 then the names of the files are derived from    \n");
		MrBayesPrint ("                    'Filename' by adding '.tree1.t', '.tree2.t', etc. If there   \n");
		MrBayesPrint ("                    are both multiple trees and multiple runs, the filenames will\n");
		MrBayesPrint ("                    be '<Filename>.tree1.run1.t', '<Filename>.tree1.run2.t', etc.\n");
		MrBayesPrint ("   Filename      -- The name of the file(s) to be summarized. This is the base of\n");
		MrBayesPrint ("                    the file name, to which endings are added according to the   \n");
		MrBayesPrint ("                    current settings of the 'Nruns' and 'Ntrees' options.        \n");
		MrBayesPrint ("   Displaygeq    -- The minimum probability of partitions to display.            \n");
		MrBayesPrint ("   Contype       -- Type of consensus tree. 'Halfcompat' results in a 50 % major-\n");
		MrBayesPrint ("                    ity rule tree, 'Allcompat' adds all compatible groups to such\n");
		MrBayesPrint ("                    a tree.                                                      \n");
		MrBayesPrint ("   Calctreeprobs -- Determines whether tree probabilities should be calculated.  \n");
		MrBayesPrint ("   Showtreeprobs -- Determines whether tree probabilities should be displayed on \n");
		MrBayesPrint ("                    screen. \n");
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("   Current settings:                                                             \n");
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("   Parameter       Options                  Current Setting                      \n");
		MrBayesPrint ("   --------------------------------------------------------                      \n");
		MrBayesPrint ("   Burnin          <number>                 %d                                   \n", sumtParams.sumtBurnIn);
		MrBayesPrint ("   Nruns           <number>                 %d                                   \n", sumtParams.numRuns);
		MrBayesPrint ("   Ntrees          <number>                 %d                                   \n", sumtParams.numTrees);
		if (sumtParams.numRuns == 1 && sumtParams.numTrees == 1)
			MrBayesPrint ("   Filename        <name>                   %s<.t>\n", sumtParams.sumtFileName);
		else if (sumtParams.numRuns == 1 && sumtParams.numTrees > 1)
			MrBayesPrint ("   Filename        <name>                   %s<.tree<i>.t>\n", sumtParams.sumtFileName);
		else if (sumtParams.numRuns > 1 && sumtParams.numTrees == 1)
			MrBayesPrint ("   Filename        <name>                   %s<.run<i>.t>\n", sumtParams.sumtFileName);
		else if (sumtParams.numRuns > 1 && sumtParams.numTrees > 1)
			MrBayesPrint ("   Filename        <name>                   %s<.tree<i>.run<i>.t>\n", sumtParams.sumtFileName);
		MrBayesPrint ("   Displaygeq      <number>                 %1.2lf                               \n", sumtParams.freqDisplay);
		MrBayesPrint ("   Contype         Halfcompat/Allcompat     %s\n", sumtParams.sumtConType);
		MrBayesPrint ("   Calctreeprobs   Yes/No                   %s                                  \n", sumtParams.calcTrprobs == YES ? "Yes" : "No");
		MrBayesPrint ("   Showtreeprobs   Yes/No                   %s                                  \n", sumtParams.showSumtTrees == YES ? "Yes" : "No");
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		}
	else if (!strcmp(helpTkn, "Tree"))
		{
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		MrBayesPrint ("   Tree                                                                          \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   This command is used by MrBayes to write trees to a nexus tree file. Trees    \n");
	    MrBayesPrint ("   are written in the Newick format. For instance,							    \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("      tree ((1,2),3,4);                                                          \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   describes an unrooted tree with taxa 1 and 2 being more closely related to    \n");
	    MrBayesPrint ("   each other than to taxa 3 and 4. If branch lengths are saved to file, they    \n");
	    MrBayesPrint ("   are given after a colon sign immediately following the terminal taxon or the  \n");
		MrBayesPrint ("   interior node they refer to. An example of an unrooted tree with branch       \n");
		MrBayesPrint ("   lengths is:                                                                   \n");
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("      tree ((1:0.064573,2:0.029042):0.041239,3:0.203988,4:0.187654);             \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("   Trees that are rooted (clock trees) are written with a basal dichotomy        \n");
	    MrBayesPrint ("   instead of a basal trichotomy. If the tree described above had been rooted    \n");
		MrBayesPrint ("   on the branch leading to taxon 4, it would have been represented as:          \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("      tree (((1,2),3),4);                                                        \n");
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		}
	else if (!strcmp(helpTkn, "Lset"))
		{
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		MrBayesPrint ("   Lset                                                                          \n");
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("   This command sets the parameters of the likelihood model. The likelihood      \n");
		MrBayesPrint ("   function is the probability of observing the data conditional on the phylo-   \n");
		MrBayesPrint ("   genetic model. In order to calculate the likelihood, you must assume a        \n");
		MrBayesPrint ("   model of character change. This command lets you tailor the biological        \n");
		MrBayesPrint ("   assumptions made in the phylogenetic model. The correct usage is              \n");
	    MrBayesPrint ("                                                                                 \n");
	    MrBayesPrint ("      lset <parameter>=<option> ... <parameter>=<option>                         \n");
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("   For example, \"lset nst=6 rates=gamma\" would set the model to a general      \n");
		MrBayesPrint ("   model of DNA substition (the GTR) with gamma-distributed rate variation       \n");
		MrBayesPrint ("   across sites.                                                                 \n");
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("   Options:                                                                      \n");
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("   Applyto   -- This option allows you to apply the lset commands to specific    \n");
		MrBayesPrint ("                partitions. This command should be the first in the list of      \n");
		MrBayesPrint ("                commands specified in lset. Moreover, it only makes sense to     \n");
		MrBayesPrint ("                be using this command if the data have been partitioned. A       \n");
		MrBayesPrint ("                default partition is set on execution of a matrix. If the data   \n");
		MrBayesPrint ("                are homogeneous (i.e., all of the same data type), then this     \n");
		MrBayesPrint ("                partition will not subdivide the characters. Up to 30 other      \n");
		MrBayesPrint ("                partitions can be defined, and you can switch among them using   \n");
		MrBayesPrint ("                \"set partition=<partition name>\". Now, you may want to         \n");
		MrBayesPrint ("                specify different models to different partitions of the data.    \n");
		MrBayesPrint ("                Applyto allows you to do this. For example, say you have         \n");
		MrBayesPrint ("                partitioned the data by codon position, and you want to apply    \n");
		MrBayesPrint ("                a nst=2 model to the first two partitions and nst=6 to the       \n");
		MrBayesPrint ("                last. This could be implemented in two uses of lset:             \n");
		MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("                   lset applyto=(1,2) nst=2                                      \n");
		MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("                   lset applyto=(3) nst=6                                        \n");
		MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("                The first applies the parameters after \"applyto\" to the        \n");
		MrBayesPrint ("                first and second partitions. The second lset applies nst=6       \n");
		MrBayesPrint ("                to the third partition. You can also use applyto=(all), which    \n");
		MrBayesPrint ("                attempts to apply the parameter settings to all of the data      \n");
		MrBayesPrint ("                partitions. Importantly, if the option is not consistent with    \n");
		MrBayesPrint ("                the data in the partition, the program will not apply the        \n");
		MrBayesPrint ("                lset option to that partition.                                   \n");
		MrBayesPrint ("   Nucmodel  -- This specifies the general form of the nucleotide substitution   \n");
		MrBayesPrint ("                model. The options are \"4by4\" [the standard model of DNA       \n");
		MrBayesPrint ("                substitution in which there are only four states (A,C,G,T/U)],   \n");
		MrBayesPrint ("                \"doublet\" (a model appropriate for modelling the stem regions  \n");
		MrBayesPrint ("                of ribosomal genes where the state space is the 16 doublets of   \n");
		MrBayesPrint ("                nucleotides), and \"codon\" (the substitution model is expanded  \n");
		MrBayesPrint ("                around triplets of nucleotides--a codon).                        \n");
		MrBayesPrint ("   Nst       -- Sets the number of substitution types: \"1\" constrains all of   \n");
		MrBayesPrint ("                the rates to be the same (e.g., a JC69 or F81 model); \"2\" all- \n");
		MrBayesPrint ("                ows transitions and transversions to have potentially different  \n");
		MrBayesPrint ("                rates (e.g., a K80 or HKY85 model); \"6\" allows all rates to    \n");
		MrBayesPrint ("                be different, subject to the constraint of time-reversibility    \n");
		MrBayesPrint ("                (e.g., a GTR model).                                             \n");
		MrBayesPrint ("   Code      -- Enforces the use of a particular genetic code. The default       \n");
		MrBayesPrint ("                is the universal code. Other options include \"vertmt\" for      \n");
		MrBayesPrint ("                vertebrate mitocondrial DNA, \"mycoplasma\", \"yeast\",          \n");
		MrBayesPrint ("                \"ciliates\", and \"metmt\" (for metazoan mitochondrial DNA      \n");
		MrBayesPrint ("                except vertebrates).                                             \n");
		MrBayesPrint ("   Ploidy    -- Specifies the ploidy of the organism. Options are \"Haploid\"    \n");
		MrBayesPrint ("                or \"Diploid\". This option is used when a coalescence prior     \n");
		MrBayesPrint ("                is used on trees.                                                \n");
		MrBayesPrint ("   Rates     -- Sets the model for among-site rate variation. In general, the    \n");
		MrBayesPrint ("                rate at a site is considered to be an unknown random variable.   \n");
		MrBayesPrint ("                The valid options are:                                           \n");
		MrBayesPrint ("                * equal    -- No rate variation across sites.                    \n");
		MrBayesPrint ("                * gamma    -- Gamma-distributed rates across sites. The rate     \n");
		MrBayesPrint ("                              at a site is drawn from a gamma distribution.      \n");
		MrBayesPrint ("                              The gamma distribution has a single parameter      \n");
		MrBayesPrint ("                              that describes how much rates vary.                \n");
		MrBayesPrint ("                * adgamma  -- Autocorrelated rates across sites. The marg-       \n");
		MrBayesPrint ("                              inal rate distribution is gamma, but adjacent      \n");
		MrBayesPrint ("                              sites have correlated rates.                       \n");
		MrBayesPrint ("                * propinv  -- A proportion of the sites are invariable.          \n");
		MrBayesPrint ("                * invgamma -- A proportion of the sites are invariable while     \n");
		MrBayesPrint ("                              the rate for the remaining sites are drawn from    \n");
		MrBayesPrint ("                              a gamma distribution.                              \n");
		MrBayesPrint ("                Note that MrBayes versions 2.0 and earlier supported options     \n");
		MrBayesPrint ("                that allowed site specific rates (e.g., ssgamma). In versions    \n");
		MrBayesPrint ("                3.0 and later, site specific rates are allowed, but set using    \n");
		MrBayesPrint ("                the 'prset ratepr' command for each partition.                   \n");
		MrBayesPrint ("   Ngammacat -- Sets the number of rate categories for the gamma distribution.   \n");
		MrBayesPrint ("                The gamma distribution is continuous. However, it is virtually   \n");
		MrBayesPrint ("                impossible to calculate likelihoods under the continuous gamma   \n");
		MrBayesPrint ("                distribution. Hence, an approximation to the continuous gamma    \n");
		MrBayesPrint ("                is used; the gamma distribution is broken into ncat categories   \n");
		MrBayesPrint ("                of equal weight (1/ncat). The mean rate for each category rep-   \n");
		MrBayesPrint ("                resents the rate for the entire cateogry. This option allows     \n");
		MrBayesPrint ("                you to specify how many rate categories to use when approx-      \n");
		MrBayesPrint ("                imating the gamma. The approximation is better as ncat is inc-   \n");
		MrBayesPrint ("                reased. In practice, \"ncat=4\" does a reasonable job of         \n");
		MrBayesPrint ("                approximating the continuous gamma.                              \n");
		MrBayesPrint ("   Nbetacat  -- Sets the number of rate categories for the beta distribution.    \n");
		MrBayesPrint ("                A symmetric beta distribution is used to model the station-      \n");
		MrBayesPrint ("                ary frequencies when morphological data are used. This option    \n");
		MrBayesPrint ("                specifies how well the beta distribution will be approx-         \n");
		MrBayesPrint ("                imated.                                                          \n");
		MrBayesPrint ("   Omegavar  -- Allows the nonsynonymous/synonymous rate ratio (omega) to vary   \n");
		MrBayesPrint ("                across codons. Ny98 assumes that there are three classes, with   \n");
		MrBayesPrint ("                potentially different omega values (omega1, omega2, omega3):     \n");
		MrBayesPrint ("                omega2 = 1; 0 < omega1 <1; and omega3 > 1. Like the Ny98 model,  \n");
		MrBayesPrint ("                the M3 model has three omega classes. However, their values are  \n");
		MrBayesPrint ("                less constrained, with omega1 < omega2 < omega3. The default     \n");
		MrBayesPrint ("                (omegavar = equal) has no variation on omega across sites.       \n");
		MrBayesPrint ("   Covarion  -- This forces the use of a covarion-like model of substitution     \n");
		MrBayesPrint ("                for nucleotide or amino acid data. The valid options are \"yes\" \n");
		MrBayesPrint ("                and \"no\". The covarion model allows the rate at a site to      \n");
		MrBayesPrint ("                change over its evolutionary history. Specifically, the site     \n");
		MrBayesPrint ("                is either on or off. When it is off, no substitutions are poss-  \n");
		MrBayesPrint ("                ible. When the process is on, substitutions occur according to   \n");
		MrBayesPrint ("                a specified substitution model (specified using the other        \n");
		MrBayesPrint ("                lset options).                                                   \n");
		MrBayesPrint ("   Coding    -- This specifies how characters were sampled. If all site pat-     \n");
		MrBayesPrint ("                terns had the possibility of being sampled, then \"all\" should  \n");
		MrBayesPrint ("                be specified (the default). Otherwise \"variable\" (only var-    \n");
		MrBayesPrint ("                iable characters had the possibility of being sampled),          \n");
		MrBayesPrint ("                \"noabsence\" (characters for which all taxa were coded as       \n");
		MrBayesPrint ("                absent were not sampled), and \"nopresence\" (characters for     \n");
		MrBayesPrint ("                which all taxa were coded as present were not sampled. \"All\"   \n");
		MrBayesPrint ("                works for all data types. However, the others only work for      \n");
		MrBayesPrint ("                morphological (all/variable) or restriction site (all/variable/  \n");
		MrBayesPrint ("                noabsence/nopresence) data.                                      \n");
		MrBayesPrint ("   Parsmodel -- This forces calculation under the so-called parsimony model      \n");
		MrBayesPrint ("                described by Tuffley and Steel (1998). The options are \"yes\"   \n");
		MrBayesPrint ("                or \"no\". Note that the biological assumptions of this model    \n");
		MrBayesPrint ("                are anything but parsimonious. In fact, this model assumes many  \n");
		MrBayesPrint ("                more parameters than the next most complicated model imple-      \n");
		MrBayesPrint ("                mented in this program. If you really believe that the pars-     \n");
		MrBayesPrint ("                imony model makes the biological assumptions described by        \n");
		MrBayesPrint ("                Tuffley and Steel, then the parsimony method is miss-named.      \n");
		/*MrBayesPrint ("   Augment   -- This allows the chain to consider the missing entries of         \n");
		MrBayesPrint ("                the data matrix as random variables. A Gibbs sampler is          \n");
		MrBayesPrint ("                used to sample states.                                           \n");*/
	    MrBayesPrint ("                                                                                 \n");
	    if (numCurrentDivisions == 0)
	    	tempInt = 1;
	    else
	    	tempInt = numCurrentDivisions;
	    for (i=0; i<tempInt; i++)
	    	{
		    if (numCurrentDivisions == 0)
				MrBayesPrint ("   Default model settings:                                                       \n");
			else
				MrBayesPrint ("   Model settings for partition %d:                                              \n", i+1);
	    	MrBayesPrint ("                                                                                 \n");
			MrBayesPrint ("   Parameter    Options                               Current Setting            \n");
			MrBayesPrint ("   ------------------------------------------------------------------            \n");		
			MrBayesPrint ("   Nucmodel     4by4/Doublet/Codon                    %s                         \n", modelParams[i].nucModel);
			MrBayesPrint ("   Nst          1/2/6                                 %s                         \n", modelParams[i].nst);
			MrBayesPrint ("   Code         Universal/Vertmt/Mycoplasma/                                     \n");
			MrBayesPrint ("                Yeast/Ciliates/Metmt                  %s                         \n", modelParams[i].geneticCode);
			MrBayesPrint ("   Ploidy       Haploid/Diploid                       %s                         \n", modelParams[i].ploidy);
			MrBayesPrint ("   Rates        Equal/Gamma/Propinv/Invgamma/Adgamma  %s                         \n", modelParams[i].ratesModel);
			MrBayesPrint ("   Ngammacat    <number>                              %d                         \n", modelParams[i].numGammaCats);
			MrBayesPrint ("   Nbetacat     <number>                              %d                         \n", modelParams[i].numBetaCats);
			MrBayesPrint ("   Omegavar     Equal/Ny98/M3                         %s                         \n", modelParams[i].omegaVar);
			MrBayesPrint ("   Covarion     No/Yes                                %s                         \n", modelParams[i].covarionModel);
			MrBayesPrint ("   Coding       All/Variable/Noabsencesites/                                     \n");
			MrBayesPrint ("                Nopresencesites                       %s                         \n", modelParams[i].coding);
			MrBayesPrint ("   Parsmodel    No/Yes                                %s                         \n", modelParams[i].parsModel);
			/*MrBayesPrint ("   Augment      No/Yes                                %s                         \n", modelParams[i].augmentData);*/
			MrBayesPrint ("   ------------------------------------------------------------------            \n");		
	    	MrBayesPrint ("                                                                                 \n");
			}
		}
	else if (!strcmp(helpTkn, "Report"))
		{
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		MrBayesPrint ("   Report                                                                        \n");
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("   This command allows you to control how the posterior distribution is          \n");
		MrBayesPrint ("   reported. For rate parameters, it allows you to choose among several popular  \n");
		MrBayesPrint ("   parameterizations. The report command also allows you to request printing of  \n");
		MrBayesPrint ("   some model aspects that are usually not reported. For instance, if a node is  \n");
		MrBayesPrint ("   constrained in the analysis, MrBayes can print the probabilities of the       \n");
		MrBayesPrint ("   ancestral states at that node. Similarly, if there is rate variation in the   \n");
		MrBayesPrint ("   model, MrBayes can print the inferred site rates, and if there is omega varia-\n");
		MrBayesPrint ("   tion, MrBayes can print the inferred omega (positive selection) values for    \n");
		MrBayesPrint ("   each codon. In a complex model with several partitions, each partition is     \n");
		MrBayesPrint ("   controlled separately using the same 'Applyto' mechanism as in the 'Lset' and \n");
		MrBayesPrint ("   'Prset' commands.                                                             \n");
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("   Options:                                                                      \n");
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("   Applyto   -- This option allows you to apply the report commands to specific  \n");
		MrBayesPrint ("                partitions. This command should be the first in the list of      \n");
		MrBayesPrint ("                commands specified in 'report'.                                  \n");
		MrBayesPrint ("                For example,                                                     \n");
		MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("                   report applyto=(1,2) tratio=ratio                             \n");
		MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("                   report applyto=(3) tratio=dirichlet                           \n");
		MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("                would result in the transition and transversion rates of the     \n");
		MrBayesPrint ("                first and second partitions in the model being reported as a     \n");
		MrBayesPrint ("                ratio and the transition and transversion rates of the third     \n");
		MrBayesPrint ("                partition being reported as proportions of the rate sum (the     \n");
		MrBayesPrint ("                Dirichlet parameterization).                                     \n");
		MrBayesPrint ("   Tratio    -- This specifies the report format for the transition and trans-   \n");
		MrBayesPrint ("                version rates of a nucleotide substituion model with nst=2.      \n");
		MrBayesPrint ("                If 'ratio' is selected, the rates will be reported as a ratio    \n");
		MrBayesPrint ("                (transition rate/transversion rate). If 'dirichlet' is selected, \n");
		MrBayesPrint ("                the transition and transversion rates will instead be reported   \n");
		MrBayesPrint ("                as proportions of the rate sum. For example, if the transition   \n");
		MrBayesPrint ("                rate is three times the transversion rate and 'ratio' is selec-  \n");
		MrBayesPrint ("                ted, this will reported as a single value, '3.0'. If 'dirichlet' \n");
		MrBayesPrint ("                is selected instead, the same rates will be reported using two   \n");
		MrBayesPrint ("                values, '0.75 0.25'. The sum of the Dirichlet values is always 1.\n");
		MrBayesPrint ("                Although the Dirichlet format may be unfamiliar to some users,   \n");
		MrBayesPrint ("                it is more convenient for specifying priors than the ratio       \n");
		MrBayesPrint ("                format.                                                          \n");
		MrBayesPrint ("   Revmat    -- This specifies the report format for the substitution rates of   \n");
		MrBayesPrint ("                a GTR substitution model for nucleotide or amino acid data. If   \n");
		MrBayesPrint ("                'ratio' is selected, the rates will be reported scaled to the    \n");
		MrBayesPrint ("                G-T rate (for nucleotides) or the Y-V rate (for amino acids). If \n");
		MrBayesPrint ("                'dirichlet' is specified instead, the rates are reported as pro- \n");
		MrBayesPrint ("                portions of the rate sum. For instance, assume that the C-T rate \n");
		MrBayesPrint ("                is twice the A-G rate and four times the transversion rates,     \n");
		MrBayesPrint ("                which are equal. If the report format is set to 'ratio', this    \n");
		MrBayesPrint ("                would be reported as '1.0 2.0 1.0 1.0 4.0 1.0' since the rates   \n");
		MrBayesPrint ("                are reported in the order rAC, rAG, rAT, rCG, rCT, rGT and scaled\n");
		MrBayesPrint ("                relative to the last rate, the G-T rate. If 'dirichlet' is selec-\n");
		MrBayesPrint ("                ted instead, the same rates would have been reported as '0.1 0.2 \n");
		MrBayesPrint ("                0.1 0.1 0.4 0.1' since the rates are now scaled so that they sum \n");
		MrBayesPrint ("                to 1.0. The Dirichlet format is the parameterization used for    \n");
		MrBayesPrint ("                formulating priors on the rates.                                 \n");
		MrBayesPrint ("   Ratemult  -- This specifies the report format used for the rate multiplier of \n");
		MrBayesPrint ("                different model partitions. Three formats are available. If      \n");
		MrBayesPrint ("                'scaled' is selected, then rates are scaled such that the mean   \n");
		MrBayesPrint ("                rate per site across partitions is 1.0. If 'ratio' is chosen,    \n");
		MrBayesPrint ("                the rates are scaled relative to the rate of the first parti-    \n");
		MrBayesPrint ("                tion. Finally, if 'dirichlet' is chosen, the rates are given as  \n");
		MrBayesPrint ("                proportions of the rate sum. The latter is the format used       \n");
		MrBayesPrint ("                when formulating priors on the rate multiplier.                  \n");
		MrBayesPrint ("   Ancstates -- If this option is set to 'yes', MrBayes will print the pro-      \n");
		MrBayesPrint ("                bability of the ancestral states at all constrained nodes. Typ-  \n");
		MrBayesPrint ("                ically, you are interested in the ancestral states of only a few \n");
		MrBayesPrint ("                characters and only at one node in the tree. To perform such     \n");
		MrBayesPrint ("                an analysis, first define and enforce a topology constraint      \n");
		MrBayesPrint ("                using 'constraint' and 'prset topologypr = constraints (...)'.   \n");
		MrBayesPrint ("                Then put the character(s) of interest in a separate partition and\n");
		MrBayesPrint ("                set MrBayes to report the ancestral states for that partition.   \n");
		MrBayesPrint ("                For instance, if the characters of interest are in partition 2,  \n");
		MrBayesPrint ("                use 'report applyto=(2) ancstates=yes' to force MrBayes to print \n");
		MrBayesPrint ("                the probability of the ancestral states of those characters at   \n");
		MrBayesPrint ("                the constrained node to the '.p' file.                           \n");
		MrBayesPrint ("   Siterates -- If this option is set to 'yes' and the relevant model has rate   \n");
		MrBayesPrint ("                variation across sites, the mean site rate in the posterior will \n");
		MrBayesPrint ("                be reported for each site to the '.p' file.                      \n");
		MrBayesPrint ("   Possel    -- If this option is set to 'yes' and the relevant model has omega  \n");
		MrBayesPrint ("                variation across sites, the mean omega value for each model site \n");
		MrBayesPrint ("                (codon in this case) will be written to the '.p' file.           \n");
		MrBayesPrint ("                                                                                 \n");
	    if (numCurrentDivisions == 0)
	    	tempInt = 1;
	    else
	    	tempInt = numCurrentDivisions;
	    for (i=0; i<tempInt; i++)
	    	{
		    if (numCurrentDivisions == 0)
				MrBayesPrint ("   Current settings:                                                         \n");
			else
				MrBayesPrint ("   Current settings for partition %d:                                              \n", i+1);
	    	MrBayesPrint ("                                                                                 \n");
			MrBayesPrint ("   Parameter       Options                  Current Setting                      \n");
			MrBayesPrint ("   --------------------------------------------------------                      \n");
			MrBayesPrint ("   Tratio          Ratio/Dir                %s                                   \n", modelParams[i].tratioFormat);
			MrBayesPrint ("   Revmat          Ratio/Dir                %s                                   \n", modelParams[i].revmatFormat);
			MrBayesPrint ("   Ratemult        Scaled/Ratio/Dir         %s                                   \n", modelParams[i].ratemultFormat);
			MrBayesPrint ("   Ancstates       Yes/No                   %s                                   \n", modelParams[i].inferAncStates);
			MrBayesPrint ("   Siterates       Yes/No                   %s                                   \n", modelParams[i].inferSiteRates);
			MrBayesPrint ("   Possel          Yes/No                   %s                                   \n", modelParams[i].inferPosSel);
			MrBayesPrint ("                                                                                 \n");
			MrBayesPrint ("   ------------------------------------------------------------------            \n");		
	    	MrBayesPrint ("                                                                                 \n");
			}
		}
	else if (!strcmp(helpTkn, "Manual"))
		{
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		MrBayesPrint ("   Manual                                                                          \n");
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("   This command allows you to generate a text file containing help information   \n");
		MrBayesPrint ("   on all the available commands. This text file can be used as an up-to-date    \n");
		MrBayesPrint ("   command reference. You can set the name of the text file using the            \n");
		MrBayesPrint ("   \"filename\" option; the default is \"commref_mb<version>.txt\".               \n");
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("   Parameter       Options                  Current Setting                      \n");
		MrBayesPrint ("   --------------------------------------------------------                      \n");
		MrBayesPrint ("   Filename        <name>                   %s                                   \n", manFileName);
	    MrBayesPrint ("                                                                                 \n");
		MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
		}
	else
		{
		return (ERROR);
		}
		
	return (NO_ERROR);
	
}





/* IsAmbig: This function returns YES if character is set as ambiguous
   either by using parenthetic notation or by ambiguity codes. It returns
   NO if character is unambiguous, missing or gapped */ 
int IsAmbig (int charCode, int dType)

{

	if (dType == DNA || dType == RNA || dType == STANDARD
		|| dType == RESTRICTION || dType == PROTEIN)
		{
		if (charCode != MISSING && charCode != GAP)
			if (NBits(charCode) > 1)
				return (YES);
		}
	else if (dType == CONTINUOUS)
		{
		/* do nothing, these cannot be partly ambiguous */

		}
	else
		{
		MrBayesPrint ("Unknown datatype in \"IsAmbig\"\n", spacer);
		}

	return (NO);

}





int IsArgValid (char *tk, char *validArg)

{

	int			i, j, k, tkLen, targetLen, numDiff, numStrMatches;
	char		tempStr[100];
	ParmInfoPtr	p;

	p = paramPtr;
	tkLen = (int) strlen(tk);

	numStrMatches = i = j = 0;
	do
		{
		if (p->valueList[i] == '|' || p->valueList[i] == '\0')
			{
			tempStr[j++] = '\0';
			targetLen = (int) strlen(tempStr);
			if (tkLen <= targetLen)
				{
				numDiff = 0;
				for (k=0; k<tkLen; k++)
					if (ChangeCase(tk[k]) != ChangeCase(tempStr[k]))
						numDiff++;
				if (numDiff == 0)
					{
					numStrMatches++;
					strcpy (validArg, tempStr);
					}
				}
			j = 0;
			}
		else
			tempStr[j++] = p->valueList[i];
		i++;
		}
	while (p->valueList[i] != '\0');
		
	if (numStrMatches == 0)
		{
		MrBayesPrint ("%s   No valid match for argument \"%s\"\n", spacer, tk);
		return (ERROR);
		}
	else if (numStrMatches == 1)
		{
		return (NO_ERROR);
		}
	else
		{
		MrBayesPrint ("%s   Argument \"%s\" is ambiguous\n", spacer, tk);
		return (ERROR);
		}
		
}





int IsIn (char ch, char *s)

{

	while (*s)
		{
		if (*s++ == ch)
			return 1;
		}
	return 0;

}





int IsMissing (int charCode, int dType)

{

	if (dType == DNA || dType == RNA)
		{
		if (charCode == 15 || charCode == 16)
			return (YES);
		}
	else if (dType == STANDARD || dType == PROTEIN)
		{
		if (charCode == MISSING || charCode == GAP)
			return (YES);
		}
	else if (dType == RESTRICTION)
		{
		if (charCode == 3 || charCode == 4)
			return (YES);
		}
	else if (dType == CONTINUOUS)
		{

		}
	else
		{
		MrBayesPrint ("Unknown datatype in \"IsMissing\"\n", spacer);
		}
	return (NO);

}





int IsSame (char *s1, char *s2)

{

	int			i, nDiff, isIdentical, len;
	
	isIdentical = YES;
	if (strlen(s1) != strlen(s2))
		isIdentical = NO; /* strings cannot be identical because they are different lengths */
	
	/* now, we go through both strings, one character at a time, to see if
	   any are different */
	if (strlen(s1) > strlen(s2))
		len = (int) strlen(s2);
	else
		len = (int) strlen(s1);
	i = nDiff = 0;
	while (i < len)
		{
		if (tolower(s1[i]) != tolower(s2[i]))
			nDiff++;
		i++;
		}
	if (nDiff == 0 && isIdentical == YES)
		return (SAME);
	else if (nDiff == 0 && isIdentical == NO)
		return (CONSISTENT_WITH);
	else
		return (DIFFERENT);

}





int IsWhite (char c)

{

	if (c == ' ' || c == '\t' || c == '\n' || c == '\r')
		{
		if (c == '\n' || c == '\r')
			return 2;
		return 1;
		}
	return 0;
	
}





int LineTermType (FILE *fp)

{

	int			ch, nextCh, term;

	term = LINETERM_UNIX;	/* default if no line endings are found */
	while ((ch = getc(fp)) != EOF)
		{
		if ((ch == '\n') || (ch == '\r'))
			{
			if (ch == '\n')
				term = LINETERM_UNIX;
			else /* ch = '\r' */
				{
				/* First test below handles one-line MAC file */
				if (((nextCh = getc(fp)) == EOF) || (nextCh != '\n'))
					term = LINETERM_MAC;
				else
					term = LINETERM_DOS;
				}
			break;
			}
		}
	(void)fseek(fp, 0L, 0);		/* rewind */
	
	return (term);

}





int LongestLine (FILE *fp)

{

	int			ch, lineLength, longest;
	
	longest = 0;
	lineLength = 0;
	while ((ch = fgetc(fp)) != EOF)
		{
		if ((ch == '\n') || (ch == '\r'))
			{
			if (lineLength > longest)
				longest = lineLength;
			lineLength = 0;
			}
		else
			lineLength++;
		}
	(void)fseek(fp, 0L, 0);		/* rewind */
	
	return (longest);

}





/* NBits: count bits in an int */
int NBits (int x)

{

	int n=0;

	for (n=0; x != 0; n++)
		x &= (x-1);
	
	return n;

}





int NucID (char nuc)

{

	char		n;
	
	if (nuc == 'U' || nuc == 'u')
		n = 'T';
	else
		n = nuc;

	if (n == 'A' || n == 'a')
		{
		return 1;
		}
	else if (n == 'C' || n == 'c')
		{
		return 2;
		}
	else if (n == 'G' || n == 'g')
		{
		return 4;
		}
	else if (n == 'T' || n == 't')
		{
		return 8;
		}
	else if (n == 'R' || n == 'r')
		{
		return 5;
		}
	else if (n == 'Y' || n == 'y')
		{
		return 10;
		}
	else if (n == 'M' || n == 'm')
		{
		return 3;
		}
	else if (n == 'K' || n == 'k')
		{
		return 12;
		}
	else if (n == 'S' || n == 's')
		{
		return 6;
		}
	else if (n == 'W' || n == 'w')
		{
		return 9;
		}
	else if (n == 'H' || n == 'h')
		{
		return 11;
		}
	else if (n == 'B' || n == 'b')
		{
		return 14;
		}
	else if (n == 'V' || n == 'v')
		{
		return 7;
		}
	else if (n == 'D' || n == 'd')
		{
		return 13;
		}
	else if (n == 'N' || n == 'n')
		{
		return 15;
		}
	else if (n == gapId)
		{
		return GAP;
		}
	else if (n == missingId)
		{
		return MISSING;
		}
	else
		return -1;
		
}





FILE *OpenBinaryFileR (char *name)

{

	FILE		*fp;

	if ((fp = fopen (name, "rb")) == NULL)  
		{   
		MrBayesPrint ("%s   Could not open file \"%s\"\n", spacer, name);
		return (NULL);
		}
	else
		return (fp);
	
}





FILE *OpenTextFileR (char *name)

{

	FILE		*fp;

	if ((fp = fopen (name, "r")) == NULL)  
		{   
		MrBayesPrint ("%s   Could not open file \"%s\"\n", spacer, name);
		return (NULL);
		}
	else
		return (fp);
	
}





FILE *OpenTextFileA (char *name)

{

	FILE		*fp;

	if ((fp = fopen (name, "a+")) == NULL)  
		{   
		MrBayesPrint ("%s   Could not open file \"%s\"\n", spacer, name);
		return (NULL);
		}
	else
		return (fp);
	
}





FILE *OpenTextFileW (char *name)

{

	FILE		*fp;

	if ((fp = fopen (name, "w+")) == NULL)  
		{   
		MrBayesPrint ("%s   Could not open file \"%s\"\n", spacer, name);
		return (NULL);
		}
	else
		return (fp);
	
}




/*-------| ParseCommand |------------------------------------------------
|
|   This function is used to parse a file. The expected format is:
|   
|      command parameter=value parameter=value ... ;
|
|   For example, the following is a valid line for this parser:
|
|      lset nst=2;
|
|   In some cases, however, the format is:
|
|      command stuff more stuff ... ;
|
|   For example, when reading a data file, the matrix command might be:
|
|      matrix
|         taxon_1 data
|         taxon_2 data
|         taxon_3 data
|         ;
|
|   Like before, the command and all of the stuff for that command are
|   terminated by a semicolon.
|
*/
int ParseCommand (char *s)

{

	int				rc, tokenType, inError, numMatches;
	char				errStr[100];

	tokenP = &s[0];
	
	inError = NO;
	do
		{
		/* Get the next token. A token is a valid word in a line. Token type is defined in "globals.h". */
		GetToken (&tokenType);
		if (strlen(token) > 0)
			{
#			if defined (SHOW_TOKENS)
			MrBayesPrint ("%s\n", token);
#			endif
			if (tokenType == LEFTCOMMENT)
				{
				/* If the token is a left comment "[", then we don't want to */
				/* actually process commands until we find a right comment.  */
				inComment = YES;
				numComments++;
				}
			if (inComment == NO && inForeignBlock == NO)
				{
				if (tokenType != SEMICOLON)
					{
					/* If the token is not a semicolon, then we will be processing 
					   either a command or a parameter. */
					if (expecting == Expecting(COMMAND))
						{
						/* We are expecting to find a command (defined above in "commands[]"). Find the 
						   correct command and set a pointer to that command. */
						commandPtr = NULL;
						if (FindValidCommand (token, &numMatches) == ERROR)
							{
							/* We couldn't find the command or the user did not specify enough letters
							   to unambiguously determine the command. The command pointer (commandPtr)
							   is NULL. */
							if (numMatches == 0)    
								MrBayesPrint ("%s   Could not find command \"%s\"\n", spacer, token);
							else 
								MrBayesPrint ("%s   Ambiguous command \"%s\"\n", spacer, token);
							inError = YES;
							}
						else
							{
							/* We did find a valid command. Set what we are expecting to see next. */
							expecting = commandPtr->expect;
							
							/* Check to see if we have one of the so-called special cases in which a 
							   command is not necessarily followed by a parameter (e.g., matrix). If we
							   do have a special case, then we want to set the parameter pointer (paramPtr)
							   appropriately. In this case, simply go to the first parameter in the parmList. */
							if (commandPtr->specialCmd == YES)
								{
								isFirstMatrixRead = YES;
								foundFirst = NO;
								paramPtr = paramTable + commandPtr->parmList[0];
								}
							if (strcmp(commandPtr->string, "Execute")==0)
								{
								readWord = YES;
								}
							}
						}
					else 
						{
						/* We are expecting to find a parameter or a value for the parameter, not a command. */
						if ((expecting & Expecting(PARAMETER)) == Expecting(PARAMETER) && 
						    (expecting & Expecting(tokenType)) != Expecting(tokenType))
							{
							/* Specifically, if we are here, we need to go through the parameter list,
							   checking to see if the token is a valid parameter. */
							expecting = (expecting & Expecting(PARAMETER));
							if (FindValidParam (token, &numMatches) == ERROR)
								{
								/* The token is not a valid parameter. */
								if (numMatches == 0)
									MrBayesPrint ("%s   Could not find parameter \"%s\"\n", spacer, token);
								else 
									MrBayesPrint ("%s   Ambiguous parameter \"%s\"\n", spacer, token);
								inError = YES;
								}
							else
								{
								/* The token is a valid parameter. Call the appropriate function ("DoXxxxParm"). */
								if ((paramPtr->fp)(paramPtr->string, token) == ERROR)
									{
									if (strcmp("Xxxxxxxxxx", paramPtr->string))
										MrBayesPrint ("%s   Error when setting parameter \"%s\" (1)\n", spacer, paramPtr->string);
									inError = YES;
									}
								}
							}
						else
							{
							/* Otherwise, we are expecting a value for the parameter. Call the appropriate function ("DoXxxxParm"). */
							if ((expecting & Expecting(tokenType)) != 0)
								expecting = (expecting & Expecting(tokenType));
							if ((Expecting(tokenType) & expecting) == Expecting(tokenType))
								{
								if ((paramPtr->fp)(paramPtr->string, token) == ERROR)
									{
									if (strcmp("Xxxxxxxxxx", paramPtr->string))
										MrBayesPrint ("%s   Error when setting parameter \"%s\" (2)\n", spacer, paramPtr->string);
									inError = YES;
									}
								}
							else
								{
								inError = YES;
								WhatVariableExp (expecting, errStr);
								MrBayesPrint ("%s   Expecting %s\n", spacer, errStr+1);  /* there will be an initial space in errStr so print from pos 1 */
								MrBayesPrint ("%s   Instead found '%s' in command '%s'\n", spacer, token, commandPtr->string);
								}
							}
						}
					}
				else
					{
					/* The token is a semicolon. This means that we are at the end of processing one command. We
					   need to clean things up. We do this by calling the finishing function ("DoXxxx"). */
					if ((Expecting(SEMICOLON) & expecting) == Expecting(SEMICOLON))
						{
						if (commandPtr->cmdFxnPtr != NULL)
							{
							/* Finish up the command here. */
							rc = (commandPtr->cmdFxnPtr) ();
							if (rc  == ERROR || rc == ABORT)
								{
								if (rc == ABORT)
									{
									MrBayesPrint ("   Mcmc run aborted\n");
									}
								else
									{
									MrBayesPrint ("%s   Error in command \"%s\"\n", spacer, commandPtr->string);
									inError = YES;
									}
								}
							}
						/* if the user typed "quit", then we want to bail out of this loop, with a NO_ERROR_QUIT */
						if (!strcmp(commandPtr->string, "Quit"))
							return (NO_ERROR_QUIT);
						expecting = Expecting(COMMAND);
						}
					else
						{
						inError = YES;
						WhatVariableExp (expecting, errStr);
						MrBayesPrint ("%s   Expecting %s\n", spacer, errStr);
						}
					}
				}
			/* Check to see if a comment is terminated. A comment can either be a right comment "]" or, if we were in a foreign nexus block
			   (e.g., a "paup" block) the terminating comment will be "end". */
			if (tokenType == RIGHTCOMMENT)
				{
				if (inComment == NO)
					{
					MrBayesPrint ("%s   Found \"]\", without having previously found \"[\"\n", spacer);
					inError = YES; 
					}
				numComments--;
				if (numComments == 0)
					inComment = NO;
				}
			if ((IsSame(token, "end") == SAME || IsSame(token, "endblock") == SAME) && inForeignBlock == YES)
				{
				strcpy (spacer, "");
				inForeignBlock = NO;
				}
			}
		
		} while (*token && inError == NO);
		
	if (inError == YES)
		return (ERROR);
	else
		return (NO_ERROR);
		
}





void PrintYesNo (int yn, char s[4])

{

	if (yn == YES)
		strcpy (s, "Yes");
	else
		strcpy (s, "No");
		
}





int ProtID (char aa)

{

	
	if (aa == 'A' || aa == 'a')      /* Ala */
		{
		return 1;
		}
	else if (aa == 'R' || aa == 'r') /* Arg */
		{
		return 2;
		}
	else if (aa == 'N' || aa == 'n') /* Asn */
		{
		return 4;
		}
	else if (aa == 'D' || aa == 'd') /* Asp */
		{
		return 8;
		}
	else if (aa == 'C' || aa == 'c') /* Cys */
		{
		return 16;
		}
	else if (aa == 'Q' || aa == 'q') /* Gln */
		{
		return 32;
		}
	else if (aa == 'E' || aa == 'e') /* Glu */
		{
		return 64;
		}
	else if (aa == 'G' || aa == 'g') /* Gly */
		{
		return 128;
		}
	else if (aa == 'H' || aa == 'h') /* His */
		{
		return 256;
		}
	else if (aa == 'I' || aa == 'i') /* Ile */
		{
		return 512;
		}
	else if (aa == 'L' || aa == 'l') /* Leu */
		{
		return 1024;
		}
	else if (aa == 'K' || aa == 'k') /* Lys */
		{
		return 2048;
		}
	else if (aa == 'M' || aa == 'm') /* Met */
		{
		return 4096;
		}
	else if (aa == 'F' || aa == 'f') /* Phe */
		{
		return 8192;
		}
	else if (aa == 'P' || aa == 'p') /* Pro */
		{
		return 16384;
		}
	else if (aa == 'S' || aa == 's') /* Ser */
		{
		return 32768;
		}
	else if (aa == 'T' || aa == 't') /* Thr */
		{
		return 65536;
		}
	else if (aa == 'W' || aa == 'w') /* Trp */
		{
		return 131072;
		}
	else if (aa == 'Y' || aa == 'y') /* Tyr */
		{
		return 262144;
		}
	else if (aa == 'V' || aa == 'v') /* Val */
		{
		return 524288;
		}
	else if (aa == 'X' || aa == 'x') /* Nonidentified */
		{
		return MISSING;
		}
	else if (aa == gapId)
		{
		return GAP;
		}
	else if (aa == missingId)
		{
		return MISSING;
		}
	else
		return -1;
		
}





int RemoveLastFromString (char *s1)

{

	int		i, j, numPrev, numRemoved;
	
	/* We remove the last name from the string simply by deleting the last "|". */
	   
	i = numPrev = 0;
	while (s1[i] != '\0')
		{
		if (s1[i] == '|')
			numPrev++;
		i++;
		}
		
	i = j = numRemoved = 0;
	while (s1[i] != '\0')
		{
		if (s1[i] == '|')
			j++;
		if (numPrev == j)
			{
			s1[i] = ' ';
			numRemoved++;
			break;
			}
		i++;
		}

	if (numRemoved != 1)
		{
		MrBayesPrint ("%s   Could not find name to remove\n", spacer);
		return (ERROR);
		}

	return (NO_ERROR);
	
}





int MBResID (char nuc)

{

	char		n;
	
	n = nuc;

	if (n == '0' || n == 'a' || n == 'A')
		{
		return 1;
		}
	else if (n == '1' || n == 'b' || n == 'B')
		{
		return 2;
		}
	else if (n == gapId)
		{
		return GAP;
		}
	else if (n == missingId)
		{
		return MISSING;
		}
	else
		return -1;
		
}





int SetPartitionInfo (int n)

{

	int			i, j, numDivisions, partDataType[MAX_NUM_DIVS];
	
	numDivisions = 0;
	for (i=0; i<numChar; i++)
		{
		j = charInfo[i].partitionId[n];
		
		if (j > numDivisions)
			numDivisions = j;

		partDataType[j - 1] = charInfo[i].charType;
		}

	for (i=0; i<numDivisions; i++)
		{
		modelParams[i].dataType = partDataType[i];
		modelParams[i].nStates = NumStates (i);
		if (modelParams[i].dataType == STANDARD)
			strcpy(modelParams[i].coding, "Variable"); 
		else if (modelParams[i].dataType == RESTRICTION)
			strcpy(modelParams[i].coding, "Noabsencesites"); 
		}
	
	return (numDivisions);
	
}





void SetUpParms (void)

{

	ParmInfoPtr p = paramTable;

	PARAM   (  0, "NEXUS",          DoNexusParm,       "NEXUS|\0");
	PARAM   (  1, "Data",           DoBeginParm,       "\0");
	PARAM   (  2, "Mrbayes",        DoBeginParm,       "\0");
	PARAM   (  3, "Xxxxxxxxxx",     DoBeginParm,       "\0");
	PARAM   (  4, "Ntax",           DoDimensionsParm,  "\0");
	PARAM   (  5, "Nchar",          DoDimensionsParm,  "\0");
	PARAM   (  6, "Interleave",     DoFormatParm,      "Yes|No|\0");
	PARAM   (  7, "Datatype",       DoFormatParm,      "Dna|Rna|Protein|Restriction|Standard|Continuous|Mixed|\0");
	PARAM   (  8, "Gap",            DoFormatParm,      "\0");
	PARAM   (  9, "Missing",        DoFormatParm,      "\0");
	PARAM   ( 10, "Matchchar",      DoFormatParm,      "\0");
	PARAM   ( 11, "MatrixInfo",     DoMatrixParm,      "\0");
	PARAM   ( 12, "Filename",       DoExecuteParm,     "\0");
	PARAM   ( 13, "Autoclose",      DoSetParm,         "Yes|No|\0");
	PARAM   ( 14, "Partition",      DoSetParm,         "\0");
	PARAM   ( 15, "Xxxxxxxxxx",     DoCharsetParm,     "\0");
	PARAM   ( 16, "Xxxxxxxxxx",     DoPartitionParm,   "\0");
	PARAM   ( 17, "Seed",           DoMcmcParm,        "\0");
	PARAM   ( 18, "Ngen",           DoMcmcParm,        "\0");
	PARAM   ( 19, "Samplefreq",     DoMcmcParm,        "\0");
	PARAM   ( 20, "Printfreq",      DoMcmcParm,        "\0");
	PARAM   ( 21, "Nchains",        DoMcmcParm,        "\0");
	PARAM   ( 22, "Temp",           DoMcmcParm,        "\0");
	PARAM   ( 23, "Filename",       DoMcmcParm,        "\0");
	PARAM   ( 24, "Burnin",         DoMcmcParm,        "\0");
	PARAM   ( 25, "Startingtree",   DoMcmcParm,        "Random|User|\0");
	PARAM   ( 26, "Nperts",         DoMcmcParm,        "\0");
	PARAM   ( 27, "Savebrlens",     DoMcmcParm,        "Yes|No|\0");
	PARAM   ( 28, "Nucmodel",       DoLsetParm,        "4by4|Doublet|Codon|\0");
	PARAM   ( 29, "Nst",            DoLsetParm,        "1|2|6|Mixed|\0");
	PARAM   ( 30, "Aamodel",        DoLsetParm,        "Poisson|Equalin|Jones|Dayhoff|Mtrev|Mtmam|Wag|Rtrev|Cprev|Vt|Blosum|Blossum|\0");
	PARAM   ( 31, "Parsmodel",      DoLsetParm,        "Yes|No|\0");
	PARAM   ( 32, "Omegavar",       DoLsetParm,        "Equal|Ny98|M3|M10|\0");
	PARAM   ( 33, "Code",           DoLsetParm,        "Universal|Vertmt|Mycoplasma|Yeast|Ciliates|Metmt|\0");
	PARAM   ( 34, "Coding",         DoLsetParm,        "All|Variable|Noabsencesites|Nopresencesites|Informative|\0");
	PARAM   ( 35, "Seqerror",       DoPrsetParm,       "\0");
	PARAM   ( 36, "Tratiopr",       DoPrsetParm,       "Beta|Fixed|\0");
	PARAM   ( 37, "Revmatpr",       DoPrsetParm,       "Dirichlet|Fixed|\0");
	PARAM   ( 38, "Omegapr",        DoPrsetParm,       "Uniform|Exponential|Fixed|\0");
	PARAM   ( 39, "Statefreqpr",    DoPrsetParm,       "Dirichlet|Fixed|\0");
	PARAM   ( 40, "Ngammacat",      DoLsetParm,        "\0");
	PARAM   ( 41, "Shapepr",        DoPrsetParm,       "Uniform|Exponential|Fixed|\0");
	PARAM   ( 42, "Ratecorrpr",     DoPrsetParm,       "Uniform|Fixed|\0");
	PARAM   ( 43, "Pinvarpr",       DoPrsetParm,       "Uniform|Fixed|\0");
	PARAM   ( 44, "Covswitchpr",    DoPrsetParm,       "Uniform|Exponential|Fixed|\0");
	PARAM   ( 45, "Xxxxxxxxxx",     DoExcludeParm,     "\0");
	PARAM   ( 46, "Xxxxxxxxxx",     DoIncludeParm,     "\0");
	PARAM   ( 47, "Xxxxxxxxxx",     DoDeleteParm,      "\0");
	PARAM   ( 48, "Xxxxxxxxxx",     DoRestoreParm,     "\0");
	PARAM   ( 49, "Xxxxxxxxxx",     DoTaxasetParm,     "\0");
	PARAM   ( 50, "Xxxxxxxxxx",     DoHelpParm,        "\0");
	PARAM   ( 51, "Applyto",        DoLsetParm,        "\0");
	PARAM   ( 52, "Rates",          DoLsetParm,        "Equal|Gamma|Propinv|Invgamma|Adgamma|\0");
	PARAM   ( 53, "Covarion",       DoLsetParm,        "Yes|No|\0");
	PARAM   ( 54, "Applyto",        DoPrsetParm,       "\0");
	PARAM   ( 55, "Tratio",         DoLinkParm,        "\0");
	PARAM   ( 56, "Revmat",         DoLinkParm,        "\0");
	PARAM   ( 57, "Omega",          DoLinkParm,        "\0");
	PARAM   ( 58, "Statefreq",      DoLinkParm,        "\0");
	PARAM   ( 59, "Shape",          DoLinkParm,        "\0");
	PARAM   ( 60, "Pinvar",         DoLinkParm,        "\0");
	PARAM   ( 61, "Correlation",    DoLinkParm,        "\0");
	PARAM   ( 62, "Ratemultiplier", DoLinkParm,        "\0");
	PARAM   ( 63, "Switchrates",    DoLinkParm,        "\0");
	PARAM   ( 64, "Symdirihyperpr", DoPrsetParm,       "Uniform|Exponential|Fixed|\0");
	PARAM   ( 65, "Xxxxxxxxxx",     DoCtypeParm,       "\0");
	PARAM   ( 66, "Xxxxxxxxxx",     DoConstraintsParm, "\0");
	PARAM   ( 67, "Topologypr",     DoPrsetParm,       "Uniform|Constraints|\0");
	PARAM   ( 68, "Brlenspr",       DoPrsetParm,       "Unconstrained|Clock|\0");
	PARAM   ( 69, "Speciationpr",   DoPrsetParm,       "Uniform|Exponential|Fixed|\0");
	PARAM   ( 70, "Extinctionpr",   DoPrsetParm,       "Uniform|Exponential|Fixed|\0");
	PARAM   ( 71, "Thetapr",        DoPrsetParm,       "Uniform|Exponential|Fixed|\0");
	PARAM   ( 72, "Topology",       DoLinkParm,        "\0");
	PARAM   ( 73, "Brlens",         DoLinkParm,        "\0");
	PARAM   ( 74, "Speciationrate", DoLinkParm,        "\0");
	PARAM   ( 75, "Extinctionrate", DoLinkParm,        "\0");
	PARAM   ( 76, "Theta",          DoLinkParm,        "\0");
	PARAM   ( 77, "Ratepr",         DoPrsetParm,       "Variable|Dirichlet|Fixed|\0");
	PARAM   ( 78, "Xxxxxxxxxx",     DoUserTreeParm,    "\0");
	PARAM   ( 79, "Xxxxxxxxxx",     DoOutgroupParm,    "\0");
	PARAM   ( 80, "Xxxxxxxxxx",     DoTreeParm,        "\0");
	PARAM   ( 81, "Filename",       DoSumtParm,        "\0");
	PARAM   ( 82, "Burnin",         DoSumtParm,        "\0");
	PARAM   ( 83, "Contype",        DoSumtParm,        "Halfcompat|Allcompat|\0");
	PARAM   ( 84, "Xxxxxxxxxx",     DoTranslateParm,   "\0");
	PARAM   ( 85, "Swapfreq",       DoMcmcParm,        "\0");
	PARAM   ( 86, "Start",          DoLogParm,         "\0");
	PARAM   ( 87, "Stop",           DoLogParm,         "\0");
	PARAM   ( 88, "Filename",       DoLogParm,         "\0");
	PARAM   ( 89, "Append",         DoLogParm,         "\0");
	PARAM   ( 90, "Replace",        DoLogParm,         "\0");
	PARAM   ( 91, "Nbetacat",       DoLsetParm,        "\0");
	PARAM   ( 92, "Augment",		DoLsetParm,        "Yes|No|\0");
	PARAM   ( 93, "Xxxxxxxxxx",     DoPairsParm,       "\0");
	PARAM   ( 94, "Xxxxxxxxxx",     DoBreaksParm,      "\0");
	PARAM   ( 95, "Nowarnings",     DoSetParm,         "Yes|No|\0");
	PARAM   ( 96, "Showtreeprobs",  DoSumtParm,        "Yes|No|\0");
	PARAM   ( 97, "Filename",       DoSumpParm,        "\0");
	PARAM   ( 98, "Burnin",         DoSumpParm,        "\0");
	PARAM   ( 99, "Reweight",       DoMcmcParm,        "\0");
	PARAM   (100, "Displaygeq",     DoSumtParm,        "\0");
	PARAM   (101, "Ny98omega1pr",   DoPrsetParm,       "Beta|Fixed|\0");
	PARAM   (102, "Ny98omega3pr",   DoPrsetParm,       "Uniform|Exponential|Fixed|\0");
	PARAM   (103, "Codoncatfreqs",  DoPrsetParm,       "Dirichlet|Fixed|\0");
	PARAM   (104, "Sampleprob",     DoPrsetParm,       "\0");
	PARAM   (105, "Aamodelpr",      DoPrsetParm,       "Fixed|Mixed|\0");
	PARAM   (106, "Aamodel",        DoLinkParm,        "\0");
	PARAM   (107, "Filename",       DoPlotParm,        "\0");
	PARAM   (108, "Parameter",      DoPlotParm,        "\0");
	PARAM   (109, "Match",          DoPlotParm,        "Perfect|Consistentwith|All|\0");
	PARAM   (110, "Burnin",         DoPlotParm,        "\0");
	PARAM   (111, "Brownscalepr",   DoPrsetParm,       "Uniform|Gamma|Gammamean|Fixed|\0");
	PARAM   (112, "Browncorrpr",    DoPrsetParm,       "Uniform|Fixed|\0");
	PARAM   (113, "Pbf",            DoMcmcParm,        "Yes|No|\0");
	PARAM   (114, "Pbfinitburnin",  DoMcmcParm,        "\0");
	PARAM   (115, "Pbfsamplefreq",  DoMcmcParm,        "\0");
	PARAM   (116, "Pbfsampletime",  DoMcmcParm,        "\0");
	PARAM   (117, "Pbfsampleburnin",DoMcmcParm,        "\0");
	PARAM   (118, "Growthpr",       DoPrsetParm,       "Uniform|Exponential|Fixed|Normal|\0");
	PARAM   (119, "Growthrate",     DoLinkParm,        "\0");
	PARAM   (120, "Xxxxxxxxxx",     DoCalibrateParm,   "Uniform|Offsetexponential|Fixed|Unconstrained|\0");
	PARAM   (121, "Calwaitpr",      DoPrsetParm,       "Uniform|Exponential|Fixed|\0");
	PARAM   (122, "M3omegapr",      DoPrsetParm,       "Exponential|Fixed|\0");
	PARAM   (123, "Applyto",        DoReportParm,      "\0");
	PARAM   (124, "Tratio",         DoReportParm,      "Dirichlet|Ratio|\0");
	PARAM   (125, "Revmat",         DoReportParm,      "Dirichlet|Ratio|\0");
	PARAM   (126, "Ratemult",       DoReportParm,      "Dirichlet|Scaled|Ratio|\0");
	PARAM   (127, "Filename",       DoManualParm,      "\0");
	PARAM   (128, "Filename1",      DoCompareTreeParm, "\0");
	PARAM   (129, "Filename2",      DoCompareTreeParm, "\0");
	PARAM   (130, "Outputname",     DoCompareTreeParm, "\0");
	PARAM   (131, "Burnin",         DoCompareTreeParm, "\0");
	PARAM   (132, "Ploidy",         DoLsetParm,        "Haploid|Diploid|\0");
	PARAM   (133, "Swapadjacent",   DoMcmcParm,        "Yes|No|\0");
	PARAM   (134, "Treeheightpr",   DoPrsetParm,       "Gamma|Exponential|\0");
	PARAM   (135, "Ancstates",      DoReportParm,      "Yes|No|\0");
	PARAM   (136, "Siterates",      DoReportParm,      "Yes|No|\0");
	PARAM   (137, "Possel",         DoReportParm,      "Yes|No|\0");
	PARAM   (138, "Plot",           DoSumpParm,        "Yes|No|\0");
	PARAM   (139, "Table",          DoSumpParm,        "Yes|No|\0");
	PARAM   (140, "Marglike",       DoSumpParm,        "Yes|No|\0");
	PARAM   (141, "Printtofile",    DoSumpParm,        "Yes|No|\0");
	PARAM   (142, "Outputname",     DoSumpParm,        "\0");
	PARAM   (143, "Redirect",       DoMcmcParm,        "Yes|No|\0");
	PARAM   (144, "Swapseed",       DoMcmcParm,        "\0");
	PARAM   (145, "Runidseed",      DoMcmcParm,        "\0");
	PARAM   (146, "Quitonerror",    DoSetParm,         "Yes|No|\0");
	PARAM   (147, "Printbrlens",    DoSumtParm,        "Yes|No|\0");
	PARAM   (148, "Brlensgeq",		DoSumtParm,        "\0");
	PARAM   (149, "Minpartfreq",	DoMcmcParm,        "\0");
	PARAM   (150, "Allchains",	 	DoMcmcParm,        "Yes|No|\0");
	PARAM   (151, "Mcmcdiagn",		DoMcmcParm,        "Yes|No|\0");
	PARAM   (152, "Diagnfreq",		DoMcmcParm,        "\0");
	PARAM   (153, "Nruns",			DoMcmcParm,        "\0");
	PARAM   (154, "Stoprule",	    DoMcmcParm,        "Yes|No|\0");
	PARAM   (155, "Stopval",		DoMcmcParm,        "\0");
	PARAM   (156, "Relburnin",		DoMcmcParm,        "Yes|No|\0");
	PARAM   (157, "Burninfrac",		DoMcmcParm,        "\0");
	PARAM   (158, "Allcomps",		DoMcmcParm,        "Yes|No|\0");
	PARAM   (159, "Printall",		DoMcmcParm,        "Yes|No|\0");
	PARAM   (160, "Printmax",		DoMcmcParm,        "\0");
	PARAM   (161, "Data",	        DoMcmcParm,        "Yes|No|\0");
	PARAM   (162, "Nruns",	        DoSumpParm,        "\0");
	PARAM   (163, "Allruns",	    DoSumpParm,        "Yes|No|\0");
	PARAM   (164, "Nruns",	        DoSumtParm,        "\0");
	PARAM   (165, "Ntrees",	        DoSumtParm,        "\0");
	PARAM   (166, "Calctrprobs",	DoSumtParm,        "Yes|No|\0");
	PARAM   (167, "Ordertaxa",	    DoMcmcParm,        "Yes|No|\0");
	PARAM   (168, "Ordertaxa",	    DoSumtParm,        "Yes|No|\0");
	PARAM   (169, "Aarevmatpr",     DoPrsetParm,       "Dirichlet|Fixed|\0");
	PARAM   (170, "Nswaps",         DoMcmcParm,       "\0");
	PARAM   (171, "Autooverwrite",  DoSetParm,         "Yes|No|\0");
	/* NOTE: If a change is made to the parameter table, make certain you
	         change the number of elements (now 200) in paramTable[] (global.h: may not be necessary 
	         and at the top of this file). */

}





void ShowNodes (TreeNode *p, int indent, int isThisTreeRooted)

{

	if (p != NULL)
		{
		MrBayesPrint ("   ");
		if (p->left == NULL && p->right == NULL && p->anc != NULL)
			{
			MrBayesPrint("%*cN %d (l=%d r=%d a=%d) %1.15lf (%s) ", 
			indent, ' ', Dex(p), Dex(p->left), Dex(p->right), Dex(p->anc), p->length, p->label);
			}
		else if (p->left != NULL && p->right == NULL && p->anc == NULL)
			{
			if (isThisTreeRooted == NO)
				{
				if (p->label[0] == '\0' || p->label[0] == '\n' || p->label[0] == ' ')
					MrBayesPrint("%*cN %d (l=%d r=%d a=%d) (---) ", 
					indent, ' ', Dex(p), Dex(p->left), Dex(p->right), Dex(p->anc));
				else
					MrBayesPrint("%*cN %d (l=%d r=%d a=%d) (%s) ", 
					indent, ' ', Dex(p), Dex(p->left), Dex(p->right), Dex(p->anc), p->label);
				}
			else
				{
				MrBayesPrint("%*cN %d (l=%d r=%d a=%d) X.XXXXXX ", 
				indent, ' ', Dex(p), Dex(p->left), Dex(p->right), Dex(p->anc));
				}
			}
		else
			{
			if (p->anc != NULL)
				{
				if (p->anc->anc == NULL && isThisTreeRooted == YES)
					MrBayesPrint("%*cN %d (l=%d r=%d a=%d) X.XXXXXX ", 
					indent, ' ', Dex(p), Dex(p->left), Dex(p->right), Dex(p->anc));
				else	
					MrBayesPrint("%*cN %d (l=%d r=%d a=%d) %1.15lf ", 
					indent, ' ', Dex(p), Dex(p->left), Dex(p->right), Dex(p->anc), p->length);
				}
			}
		if (isThisTreeRooted == YES)
			MrBayesPrint ("depth=%1.15lf\n", p->nodeDepth);
		else
			MrBayesPrint ("\n");
		ShowNodes (p->left,  indent + 2, isThisTreeRooted);
		ShowNodes (p->right, indent + 2, isThisTreeRooted);
		}
   
}





int ShowTree (TreeNode *r, int isThisTreeRooted, int n)

{

	int 			i, j, k, nNodes, x, nLines, nLevels, levelDepth, from, to;
	char			treeLine[SCREENWIDTH2], labelLine[100];
	TreeNode		**downPass, *p;

	/* get down pass sequence */
	if (isThisTreeRooted == YES)
		nNodes = 2 * n;
	else
		nNodes = 2 * n - 2;
	downPass = (TreeNode **)SafeMalloc((size_t) (2 * numTaxa * sizeof(TreeNode *)));
	if (!downPass)
		{
		MrBayesPrint ("   ERROR: Could not allocate downPass\n");
		return (ERROR);
		}
	i = 0;
	GetUserDownPass (r, downPass, &i);
	
	/* get coordinates */
	x = 0;
	nLines = 0;
	for (i=0; i<nNodes; i++)
		{
		p = downPass[i];
		if (p->left == NULL && p->right == NULL)
			{
			p->x = x;
			x += 2;
			p->y = 0;
			nLines += 2;
			}
		else if (p->left != NULL && p->right != NULL && p->anc != NULL)
			{
			p->x = p->left->x + (p->right->x - p->left->x) / 2;
			if (p->left->y > p->right->y)
				p->y = p->left->y + 1;
			else
				p->y = p->right->y + 1;
			}
		else
			{
			p->x = x;
			x += 2;
			p->y = 0;
			}
		} 
	/* print tree out, line-by-line */
	levelDepth = SCREENWIDTH / r->left->y;
	nLevels = r->left->y;
	for (j=0; j<=nLines-2; j++)
		{
		if (j % 2 == 0)
			{
			for (i=0; i<nNodes; i++)
				{
				p = downPass[i];
				if (p->left == NULL && p->x == j)
					{
					strcpy (labelLine, p->label);
					}
				}
			}
			
		for (i=0; i<SCREENWIDTH-1; i++)
			treeLine[i] = ' ';
		treeLine[SCREENWIDTH-1] = '\0';
		for (i=0; i<nNodes; i++)
			{
			p = downPass[i];
			if (p->anc != NULL)
				{
				if (p->anc->anc != NULL)
					{
					if (p->x == j)
						{
						from = (nLevels - p->anc->y) * levelDepth;
						to   = (nLevels - p->y) * levelDepth;
						if (p->y == 0)
							to = SCREENWIDTH-1;
						if (to >= SCREENWIDTH)
							to = SCREENWIDTH-1;
							
						for (k=from; k<to; k++)
							treeLine[k] = '-';
						if (p->anc->left == p)
							treeLine[from] = '/';
						else
							treeLine[from] = '\\';
						if (p->left != NULL)
							{
							treeLine[to] = '+';
							}
						if (p->anc->anc == r && p->anc->right == p)
							{
							if (isThisTreeRooted == NO)
								treeLine[to] = '+';
							else
								treeLine[from] = '\\';
							}
						}
					else
						{
						if (p->left != NULL && p->right != NULL)
							{
							if (j < p->x && j > p->left->x)
								{
								from = (nLevels - p->y) * levelDepth;
								treeLine[from] = '|';
								}
							else if (j > p->x && j < p->right->x && p->left != NULL)
								{
								from = (nLevels - p->y) * levelDepth;
								treeLine[from] = '|';
								}
							}
						}
					}
				else
					{
					if (p->x == j)
						{
						treeLine[0] = '|'; /* temp */
						}
					else if (j < p->x && j > p->left->x)
						{
						treeLine[0] = '|';
						}
					else if (j > p->x && j < p->right->x)
						{
						treeLine[0] = '|';
						}
					if (isThisTreeRooted == NO)
						{
						if (j > p->x && j <= nLines-2)
							treeLine[0] = '|';
						if (j == p->right->x)
							treeLine[0] = '+';
						}
					else
						{
						if (j == p->x)
							treeLine[0] = '+';
						}
					}
				}
			}
		treeLine[SCREENWIDTH-1] = '\0';
		if (j % 2 == 0)
			MrBayesPrint ("   %s %s\n", treeLine, labelLine);
		else
			MrBayesPrint ("   %s \n", treeLine);
		}

	if (isThisTreeRooted == NO)
		{
		for (i=0; i<SCREENWIDTH; i++)
			treeLine[i] = ' ';
		treeLine[SCREENWIDTH-1] = '\0';
		MrBayesPrint ("   |\n");
		for (k=0; k<SCREENWIDTH; k++)
			treeLine[k] = '-';
		treeLine[SCREENWIDTH-1] = '\0';
		treeLine[0] = '\\';
		strcpy (labelLine, r->label);
		labelLine[19] = '\0';
		MrBayesPrint ("   %s %s\n", treeLine, labelLine);
		}
	
	free (downPass);
	
	return (NO_ERROR);
	   
}





int StandID (char nuc)

{

	char		n;
	
	/* Note that if you change how many states are recognized, you need 
	   to look at IsMissing */
	n = nuc;

	if (n == '0')
		{
		return 1;
		}
	else if (n == '1')
		{
		return 2;
		}
	else if (n == '2')
		{
		return 4;
		}
	else if (n == '3')
		{
		return 8;
		}
	else if (n == '4')
		{
		return 16;
		}
	else if (n == '5')
		{
		return 32;
		}
	else if (n == '6')
		{
		return 64;
		}
	else if (n == '7')
		{
		return 128;
		}
	else if (n == '8')
		{
		return 256;
		}
	else if (n == '9')
		{
		return 512;
		}
	else if (n == missingId)
		{
		return MISSING;
		}
	else if (n == gapId)
		{
		return GAP;
		}
	else
		return -1;
		
}




int StateCode_AA (int n)

{
	if (n == 0)
		return 'A';      /* Ala */
	else if (n == 1)
		return 'R';		 /* Arg */
	else if (n == 2)
		return 'N';		 /* Asn */
	else if (n == 3)
		return 'D';		 /* Asp */
	else if (n == 4)
		return 'C';		 /* Cys */
	else if (n == 5)
		return 'Q';		 /* Gln */
	else if (n == 6)
		return 'E';		 /* Glu */
	else if (n == 7)
		return 'G';		 /* Gly */
	else if (n == 8)
		return 'H';		 /* His */
	else if (n == 9)
		return 'I';		 /* Ile */
	else if (n == 10)
		return 'L';		 /* Leu */
	else if (n == 11)
		return 'K';		 /* Lys */
	else if (n == 12)
		return 'M';		 /* Met */
	else if (n == 13)
		return 'F';		 /* Phe */
	else if (n == 14)
		return 'P';		 /* Pro */
	else if (n == 15)
		return 'S';		 /* Ser */
	else if (n == 16)
		return 'T';		 /* Thr */
	else if (n == 17)
		return 'W';		 /* Trp */
	else if (n == 18)
		return 'Y';		 /* Tyr */
	else if (n == 19)
		return 'V';		 /* Val */
	else
		return '?';
}



int StateCode_NUC4 (int n)

{
	if (n == 0)
		return 'A';
	else if (n == 1)
		return 'C';
	else if (n == 2)
		return 'G';
	else if (n == 3)
		return 'T';
	else return '?';
}




int StateCode_Std (int n)

{
	if (n <= 9 && n >= 0)
		return '0' + n;
	else return '?';
}




void WhatVariableExp (unsigned long int exp, char *st)

{

	int			n;
	
	strcpy (st, "");
	n = 0;
	if (exp == 0)
		strcat(st, " nothing");
	else
		{
		if ((exp & Expecting(COMMAND)) == Expecting(COMMAND))
			{
			strcat(st, " command");
			n++;
			}
		if ((exp & Expecting(PARAMETER)) == Expecting(PARAMETER))
			{
			if (n > 0)
				strcat(st, " or");
			strcat(st, " parameter");
			n++;
			}
		if ((exp & Expecting(EQUALSIGN)) == Expecting(EQUALSIGN))
			{
			if (n > 0)
				strcat(st, " or");
			strcat(st, " =");
			n++;
			}
		if ((exp & Expecting(COLON)) == Expecting(COLON))
			{
			if (n > 0)
				strcat(st, " or");
			strcat(st, " :");
			n++;
			}
		if ((exp & Expecting(SEMICOLON)) == Expecting(SEMICOLON))
			{
			if (n > 0)
				strcat(st, " or");
			strcat(st, " ;");
			n++;
			}
		if ((exp & Expecting(COMMA)) == Expecting(COMMA))
			{
			if (n > 0)
				strcat(st, " or");
			strcat(st, " ,");
			n++;
			}
		if ((exp & Expecting(POUNDSIGN)) == Expecting(POUNDSIGN))
			{
			if (n > 0)
				strcat(st, " or");
			strcat(st, " #");
			n++;
			}
		if ((exp & Expecting(QUESTIONMARK)) == Expecting(QUESTIONMARK))
			{
			if (n > 0)
				strcat(st, " or");
			strcat(st, " ?");
			n++;
			}
		if ((exp & Expecting(DASH)) == Expecting(DASH))
			{
			if (n > 0)
				strcat(st, " or");
			strcat(st, " -");
			n++;
			}
		if ((exp & Expecting(LEFTPAR)) == Expecting(LEFTPAR))
			{
			if (n > 0)
				strcat(st, " or");
			strcat(st, " (");
			n++;
			}
		if ((exp & Expecting(RIGHTPAR)) == Expecting(RIGHTPAR))
			{
			if (n > 0)
				strcat(st, " or");
			strcat(st, " )");
			n++;
			}
		if ((exp & Expecting(LEFTCOMMENT)) == Expecting(LEFTCOMMENT))
			{
			if (n > 0)
				strcat(st, " or");
			strcat(st, " [");
			n++;
			}
		if ((exp & Expecting(RIGHTCOMMENT)) == Expecting(RIGHTCOMMENT))
			{
			if (n > 0)
				strcat(st, " or");
			strcat(st, " ]");
			n++;
			}
		if ((exp & Expecting(ALPHA)) == Expecting(ALPHA))
			{
			if (n > 0)
				strcat(st, " or");
			strcat(st, " <name>");
			n++;
			}
		if ((exp & Expecting(NUMBER)) == Expecting(NUMBER))
			{
			if (n > 0)
				strcat(st, " or");
			strcat(st, " <number>");
			n++;
			}
		if ((exp & Expecting(RETURNSYMBOL)) == Expecting(RETURNSYMBOL))
			{
			if (n > 0)
				strcat(st, " or");
			strcat(st, " return");
			n++;
			}
		if ((exp & Expecting(ASTERISK)) == Expecting(ASTERISK))
			{
			if (n > 0)
				strcat(st, " or");
			strcat(st, " *");
			n++;
			}
		if ((exp & Expecting(BACKSLASH)) == Expecting(BACKSLASH))
			{
			if (n > 0)
				strcat(st, " or");
			strcat(st, " /");
			n++;
			}
		if ((exp & Expecting(BACKSLASH)) == Expecting(BACKSLASH))
			{
			if (n > 0)
				strcat(st, " or");
			strcat(st, " \\");
			n++;
			}
		if ((exp & Expecting(EXCLAMATIONMARK)) == Expecting(EXCLAMATIONMARK))
			{
			if (n > 0)
				strcat(st, " or");
			strcat(st, " !");
			n++;
			}
		if ((exp & Expecting(PERCENT)) == Expecting(PERCENT))
			{
			if (n > 0)
				strcat(st, " or");
			strcat(st, " %");
			n++;
			}
		if ((exp & Expecting(WEIRD)) == Expecting(WEIRD))
			{
			if (n > 0)
				strcat(st, " or");
			strcat(st, " <whatever>");
			n++;
			}
		if ((exp & Expecting(UNKNOWN_TOKEN_TYPE)) == Expecting(UNKNOWN_TOKEN_TYPE))
			{
			if (n > 0)
				strcat(st, " or");
			strcat(st, " no clue");
			n++;
			}
		}

}





char WhichAA (int x)

{

	if (x == 1)
		return ('A');
	else if (x == 2)
		return ('R');
	else if (x == 4)
		return ('N');
	else if (x == 8)
		return ('D');
	else if (x == 16)
		return ('C');
	else if (x == 32)
		return ('Q');
	else if (x == 64)
		return ('E');
	else if (x == 128)
		return ('G');
	else if (x == 256)
		return ('H');
	else if (x == 512)
		return ('I');
	else if (x == 1024)
		return ('L');
	else if (x == 2048)
		return ('K');
	else if (x == 4096)
		return ('M');
	else if (x == 8192)
		return ('F');
	else if (x == 16384)
		return ('P');
	else if (x == 32768)
		return ('S');
	else if (x == 65536)
		return ('T');
	else if (x == 131072)
		return ('W');
	else if (x == 262144)
		return ('Y');
	else if (x == 524288)
		return ('V');
	else if (x > 0 && x < 524288)
		return ('*');
	else if (x == MISSING)
		return ('?');
	else if (x == GAP)
		return ('-');
	else 
		return (' ');
		
}





MrBFlt WhichCont (int x)

{

	return ((MrBFlt)(x / 1000.0));
	
}





char WhichNuc (int x)

{

	if (x == 1)
		return ('A');
	else if (x == 2)
		return ('C');
	else if (x == 3)
		return ('M');
	else if (x == 4)
		return ('G');
	else if (x == 5)
		return ('R');
	else if (x == 6)
		return ('S');
	else if (x == 7)
		return ('V');
	else if (x == 8)
		return ('T');
	else if (x == 9)
		return ('W');
	else if (x == 10)
		return ('Y');
	else if (x == 11)
		return ('H');
	else if (x == 12)
		return ('K');
	else if (x == 13)
		return ('D');
	else if (x == 14)
		return ('B');
	else if (x == 15)
		return ('N');
	else if (x == MISSING)
		return ('?');
	else if (x == GAP)
		return ('-');
	else 
		return (' ');
		
}





char WhichRes (int x)

{

	if (x == 1)
		return ('0');
	else if (x == 2)
		return ('1');
	else if (x == 3)
		return ('*');
	else if (x == MISSING)
		return ('N');
	else if (x == GAP)
		return ('-');
	else 
		return (' ');
		
}





char WhichStand (int x)

{

	if (x == 1)
		return ('0');
	else if (x == 2)
		return ('1');
	else if (x == 4)
		return ('2');
	else if (x == 8)
		return ('3');
	else if (x == 16)
		return ('4');
	else if (x == 32)
		return ('5');
	else if (x == 64)
		return ('6');
	else if (x == 128)
		return ('7');
	else if (x == 256)
		return ('8');
	else if (x == 512)
		return ('9');
	else if (x > 0 && x < 512)
		return ('*');
	else if (x == MISSING)
		return ('N');
	else if (x == GAP)
		return ('-');
	else 
		return (' ');
		
}


