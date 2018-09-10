#include "SXAppendDrugInfo.h"

FILE *g_pfDrugInfo=NULL, *g_pfIn=NULL, *g_pfOut=NULL;
DrugInfo_List g_DrugInfoList;
    
static void banner(char *argv[])
{   
    printf("Synamatix SXAppendDrugInfo Copyright 2010 Synamatix Sdn Bhd\n");
    printf("Built %s %s\n", __DATE__, __TIME__);
    printf("Usage: %s\n", argv[0]);
    printf("\t<DrugInfo filename> <Input filename> <Output filename>\n\n");
}


void GenerateDrugInfoList()
{
    fprintf(stdout,"Generating DrugInList...\n");

    char acbuf[81920],*pChr=NULL; char *pcKey; stDrugInfo *pstDrugInfo; DrugInfo_List::iterator itr;

    unsigned long ulTest1,ulTest2; ulTest1=ulTest2=0;

    while(!feof(g_pfDrugInfo))
    {
        if (++ulTest1 == 1000000)
        {
            fprintf(stdout,"Total Recs Processed = %luM\n", ++ulTest2);
            ulTest1=0;
        }

        if (!fgets(acbuf,sizeof(acbuf),g_pfDrugInfo)){break;}

        pChr = strchr(acbuf,'\r'); if (pChr) {*pChr='\0';}
        pChr = strchr(acbuf,'\n'); if (pChr) {*pChr='\0';}

        if (acbuf[0] == '\0'||acbuf[0] == '#') continue;

        pChr = strtok(acbuf,"\t"); //Chromosome        
        pChr = strtok(NULL,"\t");  //Start
        pChr = strtok(NULL,"\t");  //Stop
        pChr = strtok(NULL,"\t");  //Strand
        pChr = strtok(NULL,"\t");  //Observed Variant Alleles
        pChr = strtok(NULL,"\t");  //Reference base (NCBI)
        pChr = strtok(NULL,"\t");  //Class
        pChr = strtok(NULL,"\t");  //Functional Localization
        pChr = strtok(NULL,"\t");  //JSNP ID
        pChr = strtok(NULL,"\t");  //HGVBaseG2P ID
        pChr = strtok(NULL,"\t");  //Reference / Other Allele (HapMap)
        pChr = strtok(NULL,"\t");  //African Ancestry in SW USA [ASW]
        pChr = strtok(NULL,"\t");  //CEPH Collection [CEU]
        pChr = strtok(NULL,"\t");  //Han Chinese in Beijing, China[CHB]
        pChr = strtok(NULL,"\t");  //Chinese in Metropolitan Denver,CO,USA [CHD]
        pChr = strtok(NULL,"\t");  //Gujarati Indians in Houston, TX, USA[GIH]
        pChr = strtok(NULL,"\t");  //Japanese in Tokyo, Japan [JPT]
        pChr = strtok(NULL,"\t");  //Luhya in Webuye, Kenya[LWK]
        pChr = strtok(NULL,"\t");  //Mexican Ancestry in LA, CA, USA[MEX]        
        pChr = strtok(NULL,"\t");  //Masaai in Kinyawa, Kenya[MKK]
        pChr = strtok(NULL,"\t");  //Toscani in Italia [TSI]
        pChr = strtok(NULL,"\t");  //Yoruba in Ibadan, Nigeria[YRI]

        pChr = strtok(NULL,"\t");  //dbSNP rsID
        pcKey = new char[strlen(pChr)+1]; strcpy(pcKey,pChr);

        itr = g_DrugInfoList.find(pcKey);

        if (itr == g_DrugInfoList.end()){  
           pChr = strtok(NULL,"\t");  //(Phenotype)GWAS

           pstDrugInfo = new stDrugInfo;

           pChr = strtok(NULL,"\t");  //Evidence
           pstDrugInfo->pcEvidence = new char[strlen(pChr)+1]; strcpy(pstDrugInfo->pcEvidence,pChr);

           pChr = strtok(NULL,"\t");  //Annotation
           pstDrugInfo->pcAnnotation = new char[strlen(pChr)+1]; strcpy(pstDrugInfo->pcAnnotation,pChr);

           pChr = strtok(NULL,"\t"); //Drugs

           if (pChr==NULL) {
               pstDrugInfo->pcDrugs = new char[2]; strcpy(pstDrugInfo->pcDrugs,"-");
               pstDrugInfo->pcDrugClasses = new char[2]; strcpy(pstDrugInfo->pcDrugClasses,"-");
               pstDrugInfo->pcDisease = new char[2]; strcpy(pstDrugInfo->pcDisease,"-");                 
           }
           else{
	       pstDrugInfo->pcDrugs = new char[strlen(pChr)+1]; strcpy(pstDrugInfo->pcDrugs,pChr);
 
               pChr = strtok(NULL,"\t");  //Drug Classes
               pstDrugInfo->pcDrugClasses = new char[strlen(pChr)+1]; strcpy(pstDrugInfo->pcDrugClasses,pChr);

               pChr = strtok(NULL,"\t");  //Disease
               pstDrugInfo->pcDisease = new char[strlen(pChr)+1]; strcpy(pstDrugInfo->pcDisease,pChr);                    
           }
           g_DrugInfoList[pcKey]=pstDrugInfo;  
        }
        else {
           delete [] pcKey;  
        }
    }
 
    // Testing ...
    /*fprintf(stdout,"size = %d\n",g_DrugInfoList.size());    
    itr = g_DrugInfoList.begin();
    while (itr!=g_DrugInfoList.end()){
          if ((*itr->second).pcEvidence[0] != '-') {
              fprintf(stdout,"%s\t%s\t%s\t%s\t%s\n",(*itr->second).pcEvidence,
                     (*itr->second).pcAnnotation,(*itr->second).pcDrugs,(*itr->second).pcDrugClasses,
                     (*itr->second).pcDisease);
          }
          itr++;   
    }

   fprintf(stdout,"End...\n");
   */
}


void ProcessRecs()
{
    fprintf(stdout,"ProcessRecs...\n");

    DrugInfo_List::iterator itr;
    char acbuf[81920], acOut[81920]; char *pChr=NULL, *pcRsid=NULL;   
    unsigned long ulTest1,ulTest2; ulTest1=ulTest2=0; bool bFound;

    while (!feof(g_pfIn))
    {
        ulTest1++;
        if (ulTest1 == 1000000)
        {
            ulTest2++;
            fprintf(stderr,"Total Recs Processed = %luM\n", ulTest2);
            ulTest1=0;
        }

        if (!fgets(acbuf,sizeof(acbuf),g_pfIn)){break;}

        pChr = strchr(acbuf,'\r'); if (pChr) {*pChr='\0';}
        pChr = strchr(acbuf,'\n'); if (pChr) {*pChr='\0';}

        if (acbuf[0] == '\0'||acbuf[0] == '#'){
           fprintf(g_pfOut,"%s\n",acbuf); continue;}

        if (acbuf[0] == 'C'){
           fprintf(g_pfOut,"%s\tEvidence\tAnnotation\tDrugs\tDrug Classes\tDisease\n",
                            acbuf); continue;
        }

        strcpy(acOut,acbuf);
        pChr = strtok(acbuf,"\t"); //Chromosome  
        pChr = strtok(NULL,"\t");  //GiNumber
        pChr = strtok(NULL,"\t");  //Offset
        pChr = strtok(NULL,"\t");  //Deleted_Base
        pChr = strtok(NULL,"\t");  //Fwd_DEL_Reads
        pChr = strtok(NULL,"\t");  //Rvs_DEL_Reads
        pChr = strtok(NULL,"\t");  //Total_DEL_Reads
        pChr = strtok(NULL,"\t");  //Total_Read_Density
        pChr = strtok(NULL,"\t");  //PEnd_Count
        pChr = strtok(NULL,"\t");  //Avg_QryPos
        pChr = strtok(NULL,"\t");  //Avg_QScore
        pChr = strtok(NULL,"\t");  //dbSNP

        fprintf(g_pfOut,"%s\t",acOut); bFound = false;
             
        if (pChr[0] != '-')
        {  
            pcRsid = new char[strlen(pChr)+1]; strcpy(pcRsid, pChr);
            pChr = strtok(pcRsid,",");

            while(pChr)        
            {  
                itr = g_DrugInfoList.find(pChr);    

                if (itr != g_DrugInfoList.end())
                {
                   fprintf(g_pfOut,"%s\t%s\t%s\t%s\t%s\n",(*itr->second).pcEvidence,(*itr->second).pcAnnotation,
                                   (*itr->second).pcDrugs,(*itr->second).pcDrugClasses,(*itr->second).pcDisease);
                   bFound = true;
                   break;
                }
                //else {fprintf(g_pfOut,"-\t-\t-\t-\t-\n");}
                pChr = strtok(NULL,",");
            }            
        }   
        //else {fprintf(g_pfOut,"-\t-\t-\t-\t-\n");}        

       if (!bFound){fprintf(g_pfOut,"-\t-\t-\t-\t-\n");}  
    }
}


void CloseFiles()
{
    if (g_pfIn) fclose(g_pfIn); if (g_pfOut) fclose(g_pfOut);
}


int main(int argc, char** argv) {

    if (argc < 4){
        fprintf(stdout,"Parameters not sufficient...\n\n");   goto Err_Exit;
    }

    g_pfDrugInfo = fopen(argv[1],"r");
    if (!g_pfDrugInfo){
        fprintf(stdout,"Failed to open %s ...\n",argv[1]); goto Err_Exit;
    }

    g_pfIn = fopen(argv[2],"r");
    if (!g_pfIn){
        fprintf(stdout,"Failed to open %s ...\n",argv[2]); goto Err_Exit;
    }

    g_pfOut = fopen(argv[3],"w");
    if (!g_pfOut){
        fprintf(stdout,"Failed to open %s ...\n",argv[3]); goto Err_Exit;
    }

    GenerateDrugInfoList(); ProcessRecs(); CloseFiles();

    return (EXIT_SUCCESS);

Err_Exit:
    CloseFiles(); banner(argv); exit(9);
}

