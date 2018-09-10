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
    printf("Synamatix ReadSNV Rev 9.0.0 Copyright 2011 Synamatix Sdn Bhd\n");
    //printf("Built %s %s\n", __DATE__, __TIME__);
    printf("Usage: %s\n", argv[0]);
    printf("\t-f <Input file name>\n"
           "\t-fP <Query Pos1> <Query Pos2> <Filtered filename>\n"         
           "\t    Data with Query Pos not within <Query Pos1> to <Query Pos2> inclusive will be filtered.\n"           
           "\t-fS <Quality Score Offset> <Quality Score1> <Quality Score2> <Filtered filename>\n"
           "\t    Data with Quality Score not within <Quality Score1> to <Quality Score2> inclusive will be filtered.\n"
           "\t-fPS <Query Pos1> <Query Pos2> <Quality Score Offset> <Quality Score1> <Quality Score2> <SNV Type -A,S,I,D>\n"
           "\t     <Filtered filename> \n"
           "\t     Data with Query Pos not within <Query Pos1> to <Query Pos2> inclusive OR \n"
           "\t     Quality Score not within <Quality Score1> to <Quality Score2> inclusive will be filtered.\n"
           "\t-fPSx <Query Pos1> <Query Pos2> <Quality Score Offset> <Quality Score1> <Quality Score2> <SNV Type -A,S,I,D>\n"
           "\t      <Filtered filename>\n"
           "\t      Data with Query Pos within <Query Pos1> to <Query Pos2> inclusive OR \n"
           "\t      Quality Score within <Quality Score1> to <Quality Score2> inclusive will be filtered.\n"
           "\t-gCF <Chromosome Freq report name>\n"
           "\t-gCO <Chromosome Offset Freq report name>\n"
           "\t-gCS <Chromosome Strand Freq report name>\n"
           "\t-gD <master-list file name> <SNV Type eg. S,I,D> <Sample-ID>\n"
           "\t-gDr <Density Files' Path> <Sample-ID>\n"
           "\t-gF <Freq report name>\n"
           "\t-gL <Binned Annotatation filename> <Binned Exception Annonatation filename> <Annotation filename>\n"
           "\t    <Min Supporting Reads to be printed> <Max Supporting Reads to be printed> <SNV Type eg. S,I,D>\n"
           "\t    <Sample-ID> <Genome version eg. 37.1> <Output Dir>\n"
           "\t-gLa <dbSNP database> <dbINS database> <dbDEL database> <Density Files' Path> <Min Supporting Reads>\n"
           "\t     <Max Supporting Reads> <Min Read Stregth> <Min Quality Score> <Sample-ID> <Genome version eg. 37.1>\n"
           "\t     <Quality Score Offset> <Human Ref. files path> <Output Dir> <Project SampleID>\n"
           "\t-gLaO <dbSNP database> <dbINDEL database> <Density Files' Path> <Min Supporting Reads>\n"
           "\t      <Max Supporting Reads> <Min Read Stregth> <Min Quality Score> <Sample-ID> <Genome version eg. 37.1>\n"
           "\t      <Quality Score Offset> <Human Ref. files path> <Output Dir> <Project SampleID>\n" 
           "\t-gS <Quality Score Offset> <Sample-ID> [Generates SNV Statistic report]\n"
           "\t-gSs <Quality Score Offset> <Sample-ID> [Generates SNP Statistic report]\n"
           "\t-gSi <Quality Score Offset> <Sample-ID> [Generates INS Statistic report]\n"
           "\t-gSd <Quality Score Offset> <Sample-ID> [Generates DEL Statistic report]\n"
           "\t-a <Ambiguity filename, used only together with gLa and gLaO>\n"
           "\t-c <Out clustered ascii file, used only together with gLa and gLaO>\n"   
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
            	fprintf(stderr, "Parameter -f requires an \"<input file name>\".\n");
		banner(argv); exit(9);
            }

            pObj->m_pcInFile = argv[++i]; continue;
        }

        if (strcmp(argv[i], "-fP") == 0) {
            if (i >= argc - 3) {
		fprintf(stderr, "Parameter -fP requires \"<Query Pos1> <Query Pos2> <filename>\" \n");
		banner(argv); exit(9);
            }

            int nQP1 = atoi(argv[++i]); int nQP2 = atoi(argv[++i]);
            const char *pcFileName = argv[++i];           
 
            pObj->setQryPos(nQP1,nQP2,pcFileName); continue;
	}

        if (strcmp(argv[i], "-fPS") == 0) {
            if (i >= argc - 7) {
		fprintf(stderr, "Invalid Parameters.\n");
		banner(argv); exit(9);
            }

            int nQP1 = atoi(argv[++i]); int nQP2 = atoi(argv[++i]);            
            int nQScoreOfs = atoi(argv[++i]); 
            int nQS1 = atoi(argv[++i]); int nQS2 = atoi(argv[++i]);
            char cType = argv[++i][0];
            const char* pcFileName = argv[++i];
            pObj->setFilterParams(nQP1,nQP2,nQScoreOfs,nQS1,nQS2,cType,pcFileName);
            pObj->m_bFilter = true;
            continue;
	}

        if (strcmp(argv[i], "-fPSx") == 0) {
            if (i >= argc - 7) {
		fprintf(stderr, "Invalid Parameters.\n");
		banner(argv); exit(9);
            }

            int nQP1 = atoi(argv[++i]); int nQP2 = atoi(argv[++i]);
            int nQScoreOfs = atoi(argv[++i]);
            int nQS1 = atoi(argv[++i]); int nQS2 = atoi(argv[++i]);
            char cType = argv[++i][0];
            const char* pcFileName = argv[++i];
            pObj->setFilterParams(nQP1,nQP2,nQScoreOfs,nQS1,nQS2,cType,pcFileName);
            pObj->m_bFilterx = true;
            continue;
	}

        if (strcmp(argv[i], "-fS") == 0) {
            if (i >= argc - 4) {
		fprintf(stderr, "Parameter -fS requires a \"<Quality Score Offset> <Quality Score1>  <Quality Score2> <filename>\" \n");
		banner(argv); exit(9);
            }

            int nQScoreOfs = atoi(argv[++i]); int nQS1 = atoi(argv[++i]); int nQS2 = atoi(argv[++i]);
            pObj->setQltyScore(nQScoreOfs, nQS1,nQS2,argv[++i]); continue;
	}

        if (strcmp(argv[i], "-gCF") == 0) {
            if (i >= argc - 1) {
		fprintf(stderr, "Parameter -gC requires a \"<file name>\".\n");
		banner(argv); exit(9);
            }
            pObj->m_pcChroFreqFile = argv[++i]; continue;
	}

        if (strcmp(argv[i], "-gCO") == 0) {
            if (i >= argc - 1) {
		fprintf(stderr, "Parameter -gCO requires a \"<file name>\".\n");
		banner(argv); exit(9);
            }
            pObj->m_pcChromOffsetFile = argv[++i]; continue;
	}

        if (strcmp(argv[i], "-gCS") == 0) {
            if (i >= argc - 1) {
		fprintf(stderr, "Parameter -gCS requires a \"<file name>\".\n");
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

        if (strcmp(argv[i], "-gF") == 0) {
            if (i >= argc - 1) {
		fprintf(stderr, "Parameter -gF requires a \"<file name>\".\n");
		banner(argv); exit(9);
            }
            pObj->m_pcFreqFile = argv[++i]; continue;
	}

        if (strcmp(argv[i], "-gL") == 0) {
           if (i >= argc - 9) {
		fprintf(stderr, "Invalid Parameters.\n");
		banner(argv); exit(9);
            }
            const char *pcBinnedFile = argv[++i];
            const char *pcExceptionFile = argv[++i];
            const char *pcAnnotateFile = argv[++i];
            unsigned int unMinSupportingReads = abs(atoi(argv[++i]));
            unsigned int unMaxSupportingReads = abs(atoi(argv[++i]));
            char cType = argv[++i][0];
            const char *pcSampleID = argv[++i];
            const char *pcGenomeVer = argv[++i];
            const char *pcOutputDir = argv[++i];

            pObj->setSNVLstInputs(pcBinnedFile,pcExceptionFile,pcAnnotateFile,unMinSupportingReads,
                                  unMaxSupportingReads,cType,pcSampleID,pcGenomeVer,pcOutputDir);
            continue;
	}

        if (strcmp(argv[i], "-gLa") == 0) {

           if (i >= argc - 14) {
		fprintf(stderr, "Invalid Parameters.\n");
		banner(argv); exit(9);
            }

            const char *pcdbSnp = argv[++i];
            const char *pcdbINS = argv[++i];
            const char *pcdbDEL = argv[++i];
            const char *pcDensFilePath = argv[++i];
            unsigned int unMinSupportingReads = abs(atoi(argv[++i]));
            unsigned int unMaxSupportingReads = abs(atoi(argv[++i]));
            float fMinReadStrength = atof(argv[++i]);
            int nMinQScore = atoi(argv[++i]);
            const char *pcSampleID = argv[++i];
            const char *pcGenomeVer = argv[++i];
            int nQScoreOfs = atoi(argv[++i]); 
            const char *pcHsRefFilePath = argv[++i];
            const char *pcOutputDir = argv[++i];
            const char *pcProjSampleID = argv[++i];

            pObj->setAvgQSLstInputs(pcdbSnp,pcdbINS,pcdbDEL,pcDensFilePath,unMinSupportingReads,
                                    unMaxSupportingReads,fMinReadStrength,nMinQScore,pcSampleID,
                                    pcGenomeVer,nQScoreOfs,pcHsRefFilePath,pcOutputDir,pcProjSampleID);
            continue;
	}

        if (strcmp(argv[i], "-gLaO") == 0) {

           if (i >= argc - 13) {
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
            const char *pcGenomeVer = argv[++i];
            int nQScoreOfs = atoi(argv[++i]);
            const char *pcHsRefFilePath = argv[++i];
            const char *pcOutputDir = argv[++i]; 
            const char *pcProjSampleID = argv[++i];

            pObj->setInputs_gLaO(pcdbSnp,pcdbIndel,pcDensFilePath,unMinSupportingReads,
                                 unMaxSupportingReads,fMinReadStrength,nMinQScore,pcSampleID,
                                 pcGenomeVer,nQScoreOfs,pcHsRefFilePath,pcOutputDir,
                                 pcProjSampleID);
            continue;
        }

        if (strcmp(argv[i], "-gS") == 0) {
            if (i >= argc - 2) {
                fprintf(stderr, "Invalid Parameters.\n");
		banner(argv); exit(9);
            }

            struct tm *ptm; time_t ttime;
            (void) time(&ttime); ptm = localtime(&ttime);

            pObj->m_nQScoreOfs = atoi(argv[++i]);  
            pObj->m_pcSampleID = argv[++i];
            char *pcFileName = new char[24+strlen(pObj->m_pcSampleID)];

            sprintf(pcFileName,"%s_snv_stat_qs_%02d%02d%02d.mtx",pObj->m_pcSampleID,
                       ptm->tm_year-100,ptm->tm_mon+1,ptm->tm_mday);

            pObj->m_pcStatsFile = pcFileName; continue;
	}

        if (strcmp(argv[i], "-gSd") == 0) {
            if (i >= argc - 2) {
		//fprintf(stderr, "Parameter -gSd requires a \"file name\".\n");
                fprintf(stderr, "Invalid Parameters.\n");
		banner(argv); exit(9);
            }

            //pObj->m_pcDELStatsFile = argv[++i]; continue;

            struct tm *ptm; time_t ttime;
            (void) time(&ttime); ptm = localtime(&ttime);

            pObj->m_nQScoreOfs = atoi(argv[++i]); 
            pObj->m_pcSampleID = argv[++i];
            char *pcFileName = new char[24+strlen(pObj->m_pcSampleID)];

            sprintf(pcFileName,"%s_del_stat_qs_%02d%02d%02d.mtx",pObj->m_pcSampleID,
                       ptm->tm_year-100,ptm->tm_mon+1,ptm->tm_mday);

            pObj->m_pcDELStatsFile = pcFileName; continue;
	}

        if (strcmp(argv[i], "-gSi") == 0) {
            if (i >= argc - 2) {
		//fprintf(stderr, "Parameter -gSi requires a \"file name\".\n");
                fprintf(stderr, "Invalid Parameters.\n");
		banner(argv); exit(9);
            }

            //pObj->m_pcINSStatsFile = argv[++i]; continue;

            struct tm *ptm; time_t ttime;
            (void) time(&ttime); ptm = localtime(&ttime);

            pObj->m_nQScoreOfs = atoi(argv[++i]);
            pObj->m_pcSampleID = argv[++i];
            char *pcFileName = new char[24+strlen(pObj->m_pcSampleID)];

            sprintf(pcFileName,"%s_ins_stat_qs_%02d%02d%02d.mtx",pObj->m_pcSampleID,
                       ptm->tm_year-100,ptm->tm_mon+1,ptm->tm_mday);

            pObj->m_pcINSStatsFile = pcFileName; continue;
	}

        if (strcmp(argv[i], "-gSs") == 0) {
            if (i >= argc - 2) {
		//fprintf(stderr, "Parameter -gSs requires a \"file name\".\n");
                fprintf(stderr, "Invalid Parameters.\n");
		banner(argv); exit(9);
            }

            //pObj->m_pcSNPStatsFile = argv[++i]; continue;

            struct tm *ptm; time_t ttime;
            (void) time(&ttime); ptm = localtime(&ttime);

            pObj->m_nQScoreOfs = atoi(argv[++i]); 
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

        if (strcmp(argv[i], "-v") == 0) {pObj->m_bOut2Screen = true; continue;}

        if (strcmp(argv[i], "-a") == 0) {
           if (i >= argc - 1) {
                fprintf(stderr, "Parameter -o requires an \"output file name\".\n");
                banner(argv); exit(9);
           }
           pObj->m_pcAmbiguityFile = argv[++i]; 
           pObj->m_bWAmbiguityFiles = true; continue;
        }

        if (strcmp(argv[i], "-c") == 0) {pObj->m_bOutClusteredFile = true; continue;}

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

