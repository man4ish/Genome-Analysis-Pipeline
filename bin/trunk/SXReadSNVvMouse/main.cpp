/* 
 * File:   main.cpp
 * Author: Manish
 *
 * Created on April 15, 2010, 8:46 AM
 */

#include <stdlib.h>

/*
 * 
 */
#include "CSXReadSNV.h"


static void banner(char *argv[])
{
    printf("Synamatix ReadSNV Rev 4.0.0 Copyright 2010 Synamatix Sdn Bhd\n");
    printf("Built %s %s\n", __DATE__, __TIME__);
    printf("Usage: %s\n", argv[0]);
    printf("\t-f <Input file name>\n"
           "\t-fP <Query Pos1> <Query Pos2> <Filtered filename>\n"
           "\t   Data with Query Pos not within <Query Pos1> to <Query Pos2> inclusive will be filtered.\n"
           "\t-fS <Quality Score1> <Quality Score2> <Filtered filename>\n"
           "\t   Data with Quality Score not within <Quality Score1> to <Quality Score2> inclusive will be filtered.\n"
           "\t-fPS <Query Pos1> <Query Pos2> <Quality Score1> <Quality Score2> <SNV Type -A,S,I,D> <Filtered filename> \n"
           "\t   Data with Query Pos not within <Query Pos1> to <Query Pos2> inclusive OR \n"
           "\t   Quality Score not within <Quality Score1> to <Quality Score2> inclusive will be filtered.\n"
           "\t-fPSx <Query Pos1> <Query Pos2> <Quality Score1> <Quality Score2> <SNV Type -A,S,I,D> <Filtered filename>\n"
           "\t   Data with Query Pos within <Query Pos1> to <Query Pos2> inclusive OR \n"
           "\t   Quality Score within <Quality Score1> to <Quality Score2> inclusive will be filtered.\n"
           "\t-gCF <Chromosome Freq report name>\n"
           "\t-gCO <Chromosome Offset Freq report name>\n"
           "\t-gCS <Chromosome Strand Freq report name>\n"
           //"\t-gDs <SNP master-list file name> <Sample-ID>\n"
           //"\t-gDi <INDEL master-list file name> <Sample-ID>\n"
           "\t-gD <master-list file name> <SNV Type eg. S,I,D> <Sample-ID>\n"
           "\t-gDr <Density Files' Path> <Sample-ID>\n"
           "\t-gF <Freq report name>\n"

           /*
           "\t-gL <Gene database> <Exon database> <Gene Name database> <dbSNP database> <dbINDEL database> "
           " <CNV database> <Read Density Files' Path> "
           " <Min Supporting Reads> <Max Supporting Reads> <Min Avg. Quality Score> <Min Read Stregth> <Max Read Density> "
           " <Min Bin Deviation> <SNV Type eg. A,S,I,D> <Sample-ID>\n"
           */
           "\t-gL <Gene database> <GeneKeyword Database> <Exon database> <CNV database> <Min Supporting Reads to be printed>"
           " <Max Supporting Reads to be printed> <SNV Type eg. S,I,D> <Sample-ID>\n"
           "\t-gLa <dbSNP database> <dbINDEL database> <Density Files' Path> <Min Supporting Reads>"
           " <Max Supporting Reads> <Min Read Stregth> <Min Quality Score> <Sample-ID>\n"
           "\t-gS <Sample-ID> [Generates SNV Statistic report]\n"
           "\t-gSs <Sample-ID> [Generates SNP Statistic report]\n"
           "\t-gSi <Sample-ID> [Generates INS Statistic report]\n"
           "\t-gSd <Sample-ID> [Generates DEL Statistic report]\n"
           "\t-tSR <Number of Supporting Reads> <SNV Type eg. A,S,I,D>\n"
           "\t-o <Output ascii file name>\n"
           "\t-v [Output to Screen]\n\n");
}


static void parseCmdLine(CSXReadSNV *pObj, int argc, char *argv[])
{
    for (int i = 1; i < argc; i++)
    {
        if (strcmp(argv[i], "-f") == 0)
        {
            if (i >= argc - 1)
            {
            	fprintf(stderr, "Parameter -f requires an \"input file name\".\n");
		banner(argv); exit(9);
            }

            pObj->m_pcInFile = argv[++i]; continue;
        }

        if (strcmp(argv[i], "-fP") == 0) {
            if (i >= argc - 3) {
		fprintf(stderr, "Parameter -fP requires \"Query Pos1 Query Pos2 filename\" \n");
		banner(argv); exit(9);
            }

            int nQP1 = atoi(argv[++i]); int nQP2 = atoi(argv[++i]);
            pObj->setQryPos(nQP1,nQP2,argv[++i]); continue;
	}

        if (strcmp(argv[i], "-fPS") == 0) {
            if (i >= argc - 6) {
		fprintf(stderr, "Invalid Parameters.\n");
		banner(argv); exit(9);
            }

            int nQP1 = atoi(argv[++i]); int nQP2 = atoi(argv[++i]);
            int nQS1 = atoi(argv[++i]); int nQS2 = atoi(argv[++i]);
            char cType = argv[++i][0];
            const char* pcFileName = argv[++i];
            pObj->setFilterParams(nQP1,nQP2,nQS1,nQS2,cType,pcFileName);
            pObj->m_bFilter = true;
            continue;
	}

        if (strcmp(argv[i], "-fPSx") == 0) {
            if (i >= argc - 6) {
		fprintf(stderr, "Invalid Parameters.\n");
		banner(argv); exit(9);
            }

            int nQP1 = atoi(argv[++i]); int nQP2 = atoi(argv[++i]);
            int nQS1 = atoi(argv[++i]); int nQS2 = atoi(argv[++i]);
            char cType = argv[++i][0];
            const char* pcFileName = argv[++i];
            pObj->setFilterParams(nQP1,nQP2,nQS1,nQS2,cType,pcFileName);
            pObj->m_bFilterx = true;
            continue;
	}

        if (strcmp(argv[i], "-fS") == 0) {
            if (i >= argc - 3) {
		fprintf(stderr, "Parameter -fS requires a \"Quality Score1 Query Score2 filename\" \n");
		banner(argv); exit(9);
            }

            int nQS1 = atoi(argv[++i]); int nQS2 = atoi(argv[++i]);
            pObj->setQltyScore(nQS1,nQS2,argv[++i]); continue;
	}

        if (strcmp(argv[i], "-gCF") == 0) {
            if (i >= argc - 1) {
		fprintf(stderr, "Parameter -gC requires a \"file name\".\n");
		banner(argv); exit(9);
            }
            pObj->m_pcChroFreqFile = argv[++i]; continue;
	}

        if (strcmp(argv[i], "-gCO") == 0) {
            if (i >= argc - 1) {
		fprintf(stderr, "Parameter -gCO requires a \"file name\".\n");
		banner(argv); exit(9);
            }
            pObj->m_pcChromOffsetFile = argv[++i]; continue;
	}

        if (strcmp(argv[i], "-gCS") == 0) {
            if (i >= argc - 1) {
		fprintf(stderr, "Parameter -gCS requires a \"file name\".\n");
		banner(argv); exit(9);
            }
            pObj->m_pcChroStrandFile = argv[++i]; continue;
	}

        if (strcmp(argv[i], "-gD") == 0) {
            if (i >= argc - 3) {
                fprintf(stderr, "Invalid Parameter\n");
                banner(argv); exit(9);
            }

            const char *pcFileName = argv[++i];
            char cType = argv[++i][0];
            const char *pcSampleID = argv[++i];
            pObj->setDensityInput(pcFileName, cType, pcSampleID);
            continue;
        }

        /*  
        if (strcmp(argv[i], "-gDi") == 0) {
            if (i >= argc - 2) {
		fprintf(stderr, "Invalid Parameter\n");
		banner(argv); exit(9);
            }

            const char *pcFileName = argv[++i];
            const char *pcSampleID = argv[++i];
            pObj->setINDELDensityInput(pcFileName, pcSampleID);
            continue;
	}

        if (strcmp(argv[i], "-gDs") == 0) {
            if (i >= argc - 2) {
		fprintf(stderr, "Invalid Parameter\n");
		banner(argv); exit(9);
            }

            const char *pcFileName = argv[++i];
            const char *pcSampleID = argv[++i];
            pObj->setSNPDensityInput(pcFileName, pcSampleID);
            continue;
	}
        */

        if (strcmp(argv[i], "-gF") == 0) {
            if (i >= argc - 1) {
		fprintf(stderr, "Parameter -gF requires a \"file name\".\n");
		banner(argv); exit(9);
            }
            pObj->m_pcFreqFile = argv[++i]; continue;
	}

        if (strcmp(argv[i], "-gL") == 0) {
           if (i >= argc - 8) {
		fprintf(stderr, "Invalid Parameters.\n");
		banner(argv); exit(9);
            }
            const char *pcGene = argv[++i];
            const char *pcGeneKeyword = argv[++i];
            const char *pcExon = argv[++i];
            const char *pcCNV = argv[++i];
            unsigned int unMinSupportingReads = abs(atoi(argv[++i]));
            unsigned int unMaxSupportingReads = abs(atoi(argv[++i]));
            char cType = argv[++i][0];
            const char *pcSampleID = argv[++i];

            pObj->setSNVLstInputs(pcGene,pcGeneKeyword,pcExon,pcCNV,unMinSupportingReads,
                                  unMaxSupportingReads,cType,pcSampleID);
            continue;
	}

        if (strcmp(argv[i], "-gLa") == 0) {
           if (i >= argc - 8) {
		fprintf(stderr, "Invalid Parameters.\n");
		banner(argv); exit(9);
            }

            const char *pcdbSnp = argv[++i];
            const char *pcdbIndel = argv[++i];
            const char *pcDensFilePath = argv[++i];
            unsigned int unMinSupportingReads = abs(atoi(argv[++i]));
            unsigned int unMaxSupportingReads = abs(atoi(argv[++i]));
            float fMinReadStrength = atof(argv[++i]);
            int nMinQScore = atoi(argv[++i]);
            const char *pcSampleID = argv[++i];

            pObj->setAvgQSLstInputs(pcdbSnp,pcdbIndel,pcDensFilePath,
                                    unMinSupportingReads,unMaxSupportingReads,
                                    fMinReadStrength,nMinQScore,pcSampleID);
            continue;
	}

        if (strcmp(argv[i], "-gS") == 0) {
            if (i >= argc - 1) {
		//fprintf(stderr, "Parameter -gS requires a \"file name\".\n");
                fprintf(stderr, "Invalid Parameters.\n");
		banner(argv); exit(9);
            }

            //pObj->m_pcStatsFile = argv[++i]; continue;

            struct tm *ptm; time_t ttime;
            (void) time(&ttime); ptm = localtime(&ttime);

            pObj->m_pcSampleID = argv[++i];
            char *pcFileName = new char[24+strlen(pObj->m_pcSampleID)];

            sprintf(pcFileName,"%s_snv_stat_qs_%02d%02d%02d.mtx",pObj->m_pcSampleID,
                       ptm->tm_year-100,ptm->tm_mon+1,ptm->tm_mday);

            pObj->m_pcStatsFile = pcFileName; continue;
	}

        if (strcmp(argv[i], "-gSd") == 0) {
            if (i >= argc - 1) {
		//fprintf(stderr, "Parameter -gSd requires a \"file name\".\n");
                fprintf(stderr, "Invalid Parameters.\n");
		banner(argv); exit(9);
            }

            //pObj->m_pcDELStatsFile = argv[++i]; continue;

            struct tm *ptm; time_t ttime;
            (void) time(&ttime); ptm = localtime(&ttime);

            pObj->m_pcSampleID = argv[++i];
            char *pcFileName = new char[24+strlen(pObj->m_pcSampleID)];

            sprintf(pcFileName,"%s_del_stat_qs_%02d%02d%02d.mtx",pObj->m_pcSampleID,
                       ptm->tm_year-100,ptm->tm_mon+1,ptm->tm_mday);

            pObj->m_pcDELStatsFile = pcFileName; continue;
	}

        if (strcmp(argv[i], "-gSi") == 0) {
            if (i >= argc - 1) {
		//fprintf(stderr, "Parameter -gSi requires a \"file name\".\n");
                fprintf(stderr, "Invalid Parameters.\n");
		banner(argv); exit(9);
            }

            //pObj->m_pcINSStatsFile = argv[++i]; continue;

            struct tm *ptm; time_t ttime;
            (void) time(&ttime); ptm = localtime(&ttime);

            pObj->m_pcSampleID = argv[++i];
            char *pcFileName = new char[24+strlen(pObj->m_pcSampleID)];

            sprintf(pcFileName,"%s_ins_stat_qs_%02d%02d%02d.mtx",pObj->m_pcSampleID,
                       ptm->tm_year-100,ptm->tm_mon+1,ptm->tm_mday);

            pObj->m_pcINSStatsFile = pcFileName; continue;
	}

        if (strcmp(argv[i], "-gSs") == 0) {
            if (i >= argc - 1) {
		//fprintf(stderr, "Parameter -gSs requires a \"file name\".\n");
                fprintf(stderr, "Invalid Parameters.\n");
		banner(argv); exit(9);
            }

            //pObj->m_pcSNPStatsFile = argv[++i]; continue;

            struct tm *ptm; time_t ttime;
            (void) time(&ttime); ptm = localtime(&ttime);

            pObj->m_pcSampleID = argv[++i];
            char *pcFileName = new char[24+strlen(pObj->m_pcSampleID)];

            sprintf(pcFileName,"%s_snp_stat_qs_%02d%02d%02d.mtx",pObj->m_pcSampleID,
                       ptm->tm_year-100,ptm->tm_mon+1,ptm->tm_mday);

            pObj->m_pcSNPStatsFile = pcFileName; continue;
	}

        if (strcmp(argv[i], "-o") == 0) {
            if (i >= argc - 1) {
		fprintf(stderr, "Parameter -o requires an \"output file name\".\n");
		banner(argv); exit(9);
            }

            pObj->m_pcOutFile = argv[++i]; continue;
	}

        if (strcmp(argv[i],"-tSR")==0){
            if (i >= argc - 2){
                fprintf(stderr,"Invalid Parameters. \n");
                banner(argv); exit(9);
            }

            unsigned int unSupportingReads = abs(atoi(argv[++i]));
            char cType = argv[++i][0];
            pObj->setSNVTraceInputs(unSupportingReads,cType);
            pObj->m_bSNVTrace = true;
            continue;
        }

        if (strcmp(argv[i], "-v") == 0) {
            pObj->m_bOut2Screen = true; continue;
	}

	fprintf(stderr, "Invalid parameter: %s\n", argv[i]);
	banner(argv);
	exit(9);
    }
}




int main(int argc, char** argv) {

    if (argc < 3) {banner(argv); exit(9);}

    CSXReadSNV *pSXReadSNV = new CSXReadSNV();

    if (strcmp(argv[1], "-gDr") == 0) {
        if ( argc != 4) {
            fprintf(stderr, "Invalid Parameters.\n");
            banner(argv); exit(9);
        }
        pSXReadSNV->m_pcDensFilePath = argv[2];
        pSXReadSNV->m_pcSampleID = argv[3];
        pSXReadSNV->OutputREADDensity1K();
    }
    else{
        parseCmdLine(pSXReadSNV,argc,argv);
        pSXReadSNV->run();
    }

    if (pSXReadSNV) delete pSXReadSNV;

    return (EXIT_SUCCESS);
}

