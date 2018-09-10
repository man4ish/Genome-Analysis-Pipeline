#include "CSXReadMultiBpSNV.h"

void banner(char *argv[])
{
    printf("Synamatix SXReadMultiSNV Copyright 2010 Synamatix Sdn Bhd\n");
    printf("Built %s %s\n", __DATE__, __TIME__);
    printf("Usage: %s\n", argv[0]);
    printf("\t-f <Input filename>\n"
           "\t-prn <Output filename>\n"
           "\t-gLa <dbSNP database> <dbINDEL database> <Density Files' Path> <Min Supporting Reads>"
           " <Max Supporting Reads> <Min Read Stregth> <Min Quality Score> <Sample-ID> <Genome version eg. 37.1>\n");

}


void parseCmdLine(CSXReadMultiBpSNV *pObj, int argc, char *argv[])
{

    for (int i=1; i<argc; i++) 
    {
        if (strcmp(argv[i], "-f") == 0)
        {
            if (i >= argc - 1)
            {
                printf("Parameter -f requires an \"input file name\".\n");
                banner(argv); exit(9);
            }
           
            pObj->m_pcInFile = argv[++i]; continue; 
        }

        if (strcmp(argv[i],"-prn") == 0)  
        {
           if (i >= argc - 1)
           {
               printf("Parameter -f requires an \"output file name\".\n");
           }

           pObj->PrintSNVFileContents(argv[++i]);break;    
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
            const char *pcGenomeVer = argv[++i];
            pObj->setMbpSNVListInput(pcdbSnp,pcdbIndel,pcDensFilePath,
                                     unMinSupportingReads,unMaxSupportingReads,
                                     fMinReadStrength,nMinQScore,pcSampleID,pcGenomeVer);
            pObj->Cluster();
            continue;
        }

    }     
}


int main(int argc, char *argv[])
{
    if (argc < 5) {banner(argv); exit(9);}
    
    CSXReadMultiBpSNV *pSXReadMultiBpSNV = new CSXReadMultiBpSNV();

    parseCmdLine(pSXReadMultiBpSNV,argc,argv); 

    if (pSXReadMultiBpSNV) delete pSXReadMultiBpSNV;      
    return (EXIT_SUCCESS);
}
