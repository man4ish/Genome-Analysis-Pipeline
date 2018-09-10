#include <stdio.h>
#include <stdlib.h>
#include <string.h>


enum enmAlleleVar{A2C=0,G2C,T2C,A2G,C2G,T2G,C2A,G2A,T2A,A2T,C2T,G2T};

const short enmAlleleVar_Size = G2T+1;

typedef struct _stAlleleVar
{
    unsigned int undbSnp[2];         
    unsigned int unTsTv[2];
    unsigned int unSynoNonSyno[2];
    unsigned int unGene; 
    unsigned int unExonIntron[2];
    unsigned int unPromoter;
    unsigned int unmiRNA;
    unsigned int unUTR;     

    _stAlleleVar()
    {
        unGene=unPromoter=unmiRNA=unUTR=0;

        for (short i=0; i<2; i++)  
        {  
            undbSnp[i]=unTsTv[i]=unSynoNonSyno[i]=unExonIntron[i]=0; 
        }
    }
}stAlleleVar;


typedef struct _stStatRpt
{
    unsigned int undbSnp[2];
    unsigned int unGene[2]; 
    unsigned int unExon[2];
    unsigned int unPromoter[2];
    unsigned int unmiRNA[2];
    unsigned int unUTR[2];
    unsigned int unSynonymous[2];
    unsigned int unNonSynonymous[2];
    unsigned int unNonSense[2];
    unsigned int unMissense[2];
    unsigned int unHet[2];
    unsigned int unHom[2];   
    unsigned int unGreyArea[2];
    unsigned int unTs[2];
    unsigned int unTv[2];
    stAlleleVar AllelVar[enmAlleleVar_Size];           
    unsigned int unRef[4];
    unsigned int unVar[4];
    unsigned int unScore[10]; 

    _stStatRpt()
    {
         for(int i=0;i<2;i++) {
             undbSnp[i]=unGene[i]=unExon[i]=unPromoter[i]=unmiRNA[i]=unUTR[i]=unSynonymous[i]=unNonSynonymous[i]=
             unNonSense[i]=unMissense[i]=unHet[i]=unHom[i]=unGreyArea[i]=unTs[i]=unTv[i]=0;
         } 

         for (int i=0; i<4; i++){ unRef[i]=unVar[i]=0; }
         for (int i=0; i<10; i++){ unScore[i]=0; }
    }

} stStatRpt;


stStatRpt g_stStatRpt;

enmAlleleVar g_eAlleleVar;
char g_acCodon[28]={'A','R','N','D','C','E','Q','G','H','I','L','K','M','F','P',
                    'S','T','W','Y','V','U','O','B','Z','J','X','-','*'};

const short g_cnSize = sizeof(g_acCodon);
unsigned int g_unNSCnt[g_cnSize][g_cnSize];
int g_nMinSR=1;

FILE *g_pfIn=NULL, *g_pfOut=NULL;

static void banner(char *argv[])
{
    printf("Synamatix SXPrintStatsSNP Copyright 2011 Synamatix Sdn Bhd\n");
    printf("Built %s %s\n", __DATE__, __TIME__);
    printf("Usage: %s\n", argv[0]);
    printf("\t<Input filename> <Min. Supporting Read> <Output filename>\n\n");
}


void CloseFiles()
{    if (g_pfIn) fclose(g_pfIn); if (g_pfOut) fclose(g_pfOut);
}


inline short GetCodonIdx(char cIn)
{
    for (short i=0; i < g_cnSize; i++)
    {
        if (g_acCodon[i]==cIn) return i;
    }   

    return -1;
}

void ProcessInputFile()
{
     fprintf(stdout,"ProcessInputFile...\n");

     unsigned long ulTest1=0, ulTest2=0;
     char acAllele[4],acbuf[40960]; char *pChr=NULL;
     short nRow=0,nCol=0,nIdx=0; 
     bool bNonSynonymous, bdbSnp;  

     while (!feof(g_pfIn)){

         if (++ulTest1 == 100000) {fprintf(stdout,"Total of Recs Processed = %lux0.1M\n",++ulTest2); ulTest1=0;}

         if (!fgets(acbuf,sizeof(acbuf),g_pfIn)) break;               

         pChr = strchr(acbuf,'\r'); if (pChr) *pChr=0;
         pChr = strchr(acbuf,'\n'); if (pChr) *pChr=0;

         if (acbuf[0]=='#') { fprintf(g_pfOut,"%s\n",acbuf); continue; }
         if (acbuf[0]=='C'||acbuf[0]==0) continue;


         pChr = strtok(acbuf,"\t");//Chromosome	
         pChr = strtok(NULL,"\t"); //GiNumber	
         pChr = strtok(NULL,"\t"); //Offset	
         pChr = strtok(NULL,"\t"); //Nucleotide_Variant	
         strcpy(acAllele,pChr);

         pChr = strtok(NULL,"\t"); //Fwd_SNP_Reads      
         pChr = strtok(NULL,"\t"); //Rvs_SNP_Reads      
         pChr = strtok(NULL,"\t"); //Total_SNP_Reads    
         if (atoi(pChr) < g_nMinSR) continue;

         switch (acAllele[0]){
             case 'A': {    
                            g_stStatRpt.unRef[0]++;                                                      
                            if (acAllele[2]=='C') { g_stStatRpt.unVar[1]++; g_eAlleleVar=A2C; }
                            else if (acAllele[2]=='G') { g_stStatRpt.unVar[2]++; g_eAlleleVar=A2G; }
                            else { g_stStatRpt.unVar[3]++; g_eAlleleVar=A2T; }
                            break;
                       } 
             case 'C': {                            
                            g_stStatRpt.unRef[1]++;
                            if (acAllele[2]=='A') { g_stStatRpt.unVar[0]++; g_eAlleleVar=C2A; }
                            else if (acAllele[2]=='G') { g_stStatRpt.unVar[2]++; g_eAlleleVar=C2G; }
                            else { g_stStatRpt.unVar[3]++; g_eAlleleVar=C2T; } 
                            break;
                       }  
             case 'G': {
                            g_stStatRpt.unRef[2]++;
                            if (acAllele[2]=='A') { g_stStatRpt.unVar[0]++; g_eAlleleVar=G2A; }
                            else if (acAllele[2]=='C') { g_stStatRpt.unVar[1]++; g_eAlleleVar=G2C; }
                            else { g_stStatRpt.unVar[3]++; g_eAlleleVar=G2T; }
                            break;
                       }  
             case 'T': {
                            g_stStatRpt.unRef[3]++;
                            if (acAllele[2]=='A') { g_stStatRpt.unVar[0]++; g_eAlleleVar=T2A; }
                            else if (acAllele[2]=='C') { g_stStatRpt.unVar[1]++; g_eAlleleVar=T2C; }
                            else { g_stStatRpt.unVar[2]++; g_eAlleleVar=T2G; }                             
                            break; 
                       }
         } 

         pChr = strtok(NULL,"\t"); //Total_Read_Density	
         pChr = strtok(NULL,"\t"); //PEnd_Count	
         pChr = strtok(NULL,"\t"); //Avg_QryPos	
         pChr = strtok(NULL,"\t"); //Avg_QScore	

         pChr = strtok(NULL,"\t"); //dbSNP	

         bdbSnp = (pChr[0]!='-')? true:false;

         g_stStatRpt.undbSnp[!bdbSnp]++;
         g_stStatRpt.AllelVar[g_eAlleleVar].undbSnp[!bdbSnp]++;
                         
         pChr = strtok(NULL,"\t"); //Gene_Name
         if (pChr[0]!='-') 
         { 
             if (bdbSnp) g_stStatRpt.unGene[0]++; else g_stStatRpt.unGene[1]++; 
             g_stStatRpt.AllelVar[g_eAlleleVar].unGene++;
         }

         pChr = strtok(NULL,"\t"); //Gene_Description	
         pChr = strtok(NULL,"\t"); //Keyword	
         pChr = strtok(NULL,"\t"); //miRNA	

         if (pChr[0]!='-') 
         {             
            g_stStatRpt.unmiRNA[!bdbSnp]++; g_stStatRpt.AllelVar[g_eAlleleVar].unmiRNA++;
         }

         pChr = strtok(NULL,"\t"); //Promoter	
         if (pChr[0]=='Y') 
         { 
            g_stStatRpt.unPromoter[!bdbSnp]++; g_stStatRpt.AllelVar[g_eAlleleVar].unPromoter++;  
         }

         pChr = strtok(NULL,"\t"); //UTR	
         if (pChr[0]=='Y') 
         {
            g_stStatRpt.unUTR[!bdbSnp]++; g_stStatRpt.AllelVar[g_eAlleleVar].unUTR++; 
         }

         pChr = strtok(NULL,"\t"); //Exon	
         if (pChr[0]!='-') 
         {     
            g_stStatRpt.unExon[!bdbSnp]++; g_stStatRpt.AllelVar[g_eAlleleVar].unExonIntron[0]++;
         }
         else
         {
             g_stStatRpt.AllelVar[g_eAlleleVar].unExonIntron[1]++; 
         } 

         pChr = strtok(NULL,"\t"); //Zygosity	
         if (pChr[0]=='-') { g_stStatRpt.unGreyArea[!bdbSnp]++; }
         else if (pChr[1]=='e') { g_stStatRpt.unHet[!bdbSnp]++; }
         else { g_stStatRpt.unHom[!bdbSnp]++; }          
            
         pChr = strtok(NULL,"\t"); //Local_Copy_Number	
         pChr = strtok(NULL,"\t"); //Mappable_Bases	
         pChr = strtok(NULL,"\t"); //Mappable_Read_Densities	
         pChr = strtok(NULL,"\t"); //Mappable_Repeat_Densities	
         pChr = strtok(NULL,"\t"); //Mappable_Bases_Repeat	
         pChr = strtok(NULL,"\t"); //Repeat_Densities_MBP	
         pChr = strtok(NULL,"\t"); //Known_CNVRegion
	
         pChr = strtok(NULL,"\t"); //Transversions(Tv)_Transitions(Ts)	
         if (pChr[1]=='s') 
         {             
            g_stStatRpt.unTs[!bdbSnp]++; g_stStatRpt.AllelVar[g_eAlleleVar].unTsTv[0]++;
         }
         else if (pChr[1]=='v') 
         { 
            g_stStatRpt.unTv[!bdbSnp]++; g_stStatRpt.AllelVar[g_eAlleleVar].unTsTv[1]++;
         }
           
         pChr = strtok(NULL,"\t"); //Synonymous(S)_Non-Synonymous(NS)  

         bNonSynonymous = false;   
         if (pChr[0] == 'S') 
         { 
            g_stStatRpt.unSynonymous[!bdbSnp]++; g_stStatRpt.AllelVar[g_eAlleleVar].unSynoNonSyno[0]++;  
         }    
         else if (pChr[0] == 'N' && pChr[1] == 'S') { 
            g_stStatRpt.unNonSynonymous[!bdbSnp]++; g_stStatRpt.AllelVar[g_eAlleleVar].unSynoNonSyno[1]++; 
            bNonSynonymous = true; 
         }

         pChr = strtok(NULL,"\t"); //Protein_Variant

         if (bNonSynonymous)
         {
             nRow=GetCodonIdx(pChr[0]); nCol=GetCodonIdx(pChr[2]); g_unNSCnt[nRow][nCol]++;             
         }

         pChr = strtok(NULL,"\t"); //Missense
         if (pChr[0] == 'Y') { if (bdbSnp) g_stStatRpt.unMissense[0]++;  else g_stStatRpt.unMissense[1]++; }

         pChr = strtok(NULL,"\t"); //Nonsense
         if (pChr[0] == 'Y') { if (bdbSnp) g_stStatRpt.unNonSense[0]++; else g_stStatRpt.unNonSense[1]++; }

         pChr = strtok(NULL,"\t"); //PEScore
         pChr = strtok(NULL,"\t"); //SNPReadScore
         pChr = strtok(NULL,"\t"); //StrandScore
         pChr = strtok(NULL,"\t"); //ReadDenScore
         pChr = strtok(NULL,"\t"); //ReptDenScore
         pChr = strtok(NULL,"\t"); //AQSScore
         pChr = strtok(NULL,"\t"); //ConfidentScore

         nIdx = ((short)((atof(pChr)*100)+0.5)); 
         nIdx= (short)((((float)nIdx-0.1)/10));

         if (nIdx < 0) nIdx=0;

         g_stStatRpt.unScore[nIdx]++;
     }

     fprintf(stdout,"ProcessInputFile...done!!\n");
}


void OutputProteinVarMatrix()
{     
    fprintf(g_pfOut,"Protein Change\nFrom Protein\\To Protein");
    
    for (short i=0; i < g_cnSize; i++)
    {
        fprintf(g_pfOut,"\t%c",g_acCodon[i]);
    }

    fprintf(g_pfOut,"\n");

    for (short j=0; j < g_cnSize; j++)
    {
        fprintf(g_pfOut,"%c",g_acCodon[j]);
        for (short k=0; k < g_cnSize; k++){
            fprintf(g_pfOut,"\t%u",g_unNSCnt[j][k]);
        }  
        fprintf(g_pfOut,"\n");
    }
}


inline double Division(unsigned long ulVal, unsigned long ulDenominator)
{
    if (ulDenominator == 0) return 0;
    return (ulVal == 0)? 0:(double)ulVal/ulDenominator;
}

void OutputStats()
{
    unsigned int unTotalSNP=g_stStatRpt.undbSnp[0]+g_stStatRpt.undbSnp[1];
    unsigned int unTotalGene=g_stStatRpt.unGene[0]+g_stStatRpt.unGene[1];
    unsigned int unTotalExon=g_stStatRpt.unExon[0]+g_stStatRpt.unExon[1];
    unsigned int unTotalHet=g_stStatRpt.unHet[0]+g_stStatRpt.unHet[1];
    unsigned int unTotalHom=g_stStatRpt.unHom[0]+g_stStatRpt.unHom[1];
    unsigned int unTotalTs=g_stStatRpt.unTs[0]+g_stStatRpt.unTs[1];
    unsigned int unTotalTv=g_stStatRpt.unTv[0]+g_stStatRpt.unTv[1]; 
    
    double dbKnownSNP=(double)Division(g_stStatRpt.undbSnp[0],unTotalSNP)*100;

    unsigned int unA2C=g_stStatRpt.AllelVar[A2C].undbSnp[0]+g_stStatRpt.AllelVar[A2C].undbSnp[1];
    unsigned int unG2C=g_stStatRpt.AllelVar[G2C].undbSnp[0]+g_stStatRpt.AllelVar[G2C].undbSnp[1];
    unsigned int unT2C=g_stStatRpt.AllelVar[T2C].undbSnp[0]+g_stStatRpt.AllelVar[T2C].undbSnp[1];
    unsigned int unA2G=g_stStatRpt.AllelVar[A2G].undbSnp[0]+g_stStatRpt.AllelVar[A2G].undbSnp[1];;
    unsigned int unC2G=g_stStatRpt.AllelVar[C2G].undbSnp[0]+g_stStatRpt.AllelVar[C2G].undbSnp[1];
    unsigned int unT2G=g_stStatRpt.AllelVar[T2G].undbSnp[0]+g_stStatRpt.AllelVar[T2G].undbSnp[1];
    unsigned int unC2A=g_stStatRpt.AllelVar[C2A].undbSnp[0]+g_stStatRpt.AllelVar[C2A].undbSnp[1];
    unsigned int unG2A=g_stStatRpt.AllelVar[G2A].undbSnp[0]+g_stStatRpt.AllelVar[G2A].undbSnp[1];
    unsigned int unT2A=g_stStatRpt.AllelVar[T2A].undbSnp[0]+g_stStatRpt.AllelVar[T2A].undbSnp[1];
    unsigned int unA2T=g_stStatRpt.AllelVar[A2T].undbSnp[0]+g_stStatRpt.AllelVar[A2T].undbSnp[1];
    unsigned int unC2T=g_stStatRpt.AllelVar[C2T].undbSnp[0]+g_stStatRpt.AllelVar[C2T].undbSnp[1];
    unsigned int unG2T=g_stStatRpt.AllelVar[G2T].undbSnp[0]+g_stStatRpt.AllelVar[G2T].undbSnp[1];   

    double dbA2C=(double)Division(unA2C,unTotalSNP)*100;
    double dbG2C=(double)Division(unG2C,unTotalSNP)*100;
    double dbT2C=(double)Division(unT2C,unTotalSNP)*100;
    double dbA2G=(double)Division(unA2G,unTotalSNP)*100;
    double dbC2G=(double)Division(unC2G,unTotalSNP)*100;
    double dbT2G=(double)Division(unT2G,unTotalSNP)*100;  
    double dbC2A=(double)Division(unC2A,unTotalSNP)*100;
    double dbG2A=(double)Division(unG2A,unTotalSNP)*100;
    double dbT2A=(double)Division(unT2A,unTotalSNP)*100;
    double dbA2T=(double)Division(unA2T,unTotalSNP)*100;
    double dbC2T=(double)Division(unC2T,unTotalSNP)*100;

    double dbRefA=(double)Division(g_stStatRpt.unRef[0],unTotalSNP)*100;
    double dbRefC=(double)Division(g_stStatRpt.unRef[1],unTotalSNP)*100;
    double dbRefG=(double)Division(g_stStatRpt.unRef[2],unTotalSNP)*100; 
    double dbVarA=(double)Division(g_stStatRpt.unVar[0],unTotalSNP)*100;
    double dbVarC=(double)Division(g_stStatRpt.unVar[1],unTotalSNP)*100; 
    double dbVarG=(double)Division(g_stStatRpt.unVar[2],unTotalSNP)*100;

    float fPercent=0.0, fCumlPercent=0.0;  

    //fprintf(g_pfOut,"\nMin. Supporting Reads = %d\n",g_nMinSR);
  
    fprintf(g_pfOut,"\n\nSNP Summary\tKnown\tNovel\tTotal\n");
    fprintf(g_pfOut,"Total SNP (Count)\t%u\t%u\t%u\n",g_stStatRpt.undbSnp[0],g_stStatRpt.undbSnp[1],unTotalSNP);                    
    fprintf(g_pfOut,"Total SNP (by %%)\t%.2f%%\t%.2f%%\t100.00%%\n",dbKnownSNP,100.00-dbKnownSNP);

    fprintf(g_pfOut,"SNP in Gene (Count)\t%u\t%u\t%u\n",g_stStatRpt.unGene[0],g_stStatRpt.unGene[1],unTotalGene);
    fprintf(g_pfOut,"SNP in Gene (By %%)\t%.2f%%\t%.2f%%\t%.2f%%\n",Division(g_stStatRpt.unGene[0],unTotalSNP)*100,
                                                             Division(g_stStatRpt.unGene[1],unTotalSNP)*100,
                                                             Division(unTotalGene,unTotalSNP)*100);
                                                               
    fprintf(g_pfOut,"SNP in Exon (Count)\t%u\t%u\t%u\n",g_stStatRpt.unExon[0],g_stStatRpt.unExon[1],unTotalExon);
    fprintf(g_pfOut,"SNP in Exon (By %%)\t%.2f%%\t%.2f%%\t%.2f%%\n",Division(g_stStatRpt.unExon[0],unTotalSNP)*100,
                                                             Division(g_stStatRpt.unExon[1],unTotalSNP)*100,
                                                             Division(unTotalExon,unTotalSNP)*100); 

    fprintf(g_pfOut,"SNP in Promoter\t%u\t%u\t%u\n",g_stStatRpt.unPromoter[0],g_stStatRpt.unPromoter[1],
                                                    g_stStatRpt.unPromoter[0]+g_stStatRpt.unPromoter[1]);
    fprintf(g_pfOut,"SNP in miRNA\t%u\t%u\t%u\n",g_stStatRpt.unmiRNA[0],g_stStatRpt.unmiRNA[1],
                                                 g_stStatRpt.unmiRNA[0]+g_stStatRpt.unmiRNA[1]);
    fprintf(g_pfOut,"SNP in UTR\t%u\t%u\t%u\n",g_stStatRpt.unUTR[0],g_stStatRpt.unUTR[1],
                                               g_stStatRpt.unUTR[0]+g_stStatRpt.unUTR[1]);
    fprintf(g_pfOut,"\n");
    fprintf(g_pfOut,"Synonymous SNP\t%u\t%u\t%u\n",g_stStatRpt.unSynonymous[0],g_stStatRpt.unSynonymous[1],
                                                   g_stStatRpt.unSynonymous[0]+g_stStatRpt.unSynonymous[1]);
    fprintf(g_pfOut,"Non-Synonymous SNP\t%u\t%u\t%u\n",g_stStatRpt.unNonSynonymous[0],g_stStatRpt.unNonSynonymous[1],
                                                   g_stStatRpt.unNonSynonymous[0]+g_stStatRpt.unNonSynonymous[1]); 
    fprintf(g_pfOut,"NonSense SNP\t%u\t%u\t%u\n",g_stStatRpt.unNonSense[0],g_stStatRpt.unNonSense[1],
                                                 g_stStatRpt.unNonSense[0]+g_stStatRpt.unNonSense[1]);
    fprintf(g_pfOut,"Missense SNP\t%u\t%u\t%u\n",g_stStatRpt.unMissense[0],g_stStatRpt.unMissense[1],
                                                 g_stStatRpt.unMissense[0]+g_stStatRpt.unMissense[1]);
    fprintf(g_pfOut,"\n");
    fprintf(g_pfOut,"Heterozygous SNP(het)\t%u\t%u\t%u\n",g_stStatRpt.unHet[0],g_stStatRpt.unHet[1],unTotalHet); 
    fprintf(g_pfOut,"Homozygous SNP(hom)\t%u\t%u\t%u\n",g_stStatRpt.unHom[0],g_stStatRpt.unHom[1],unTotalHom); 
    fprintf(g_pfOut,"Grey Area\t%u\t%u\t%u\n",g_stStatRpt.unGreyArea[0],g_stStatRpt.unGreyArea[1],
                                              g_stStatRpt.unGreyArea[0]+g_stStatRpt.unGreyArea[1]);
    fprintf(g_pfOut,"Het/Hom Ratio\t%.2f\t%.2f\t%.2f\n",(double)Division(g_stStatRpt.unHet[0],g_stStatRpt.unHom[0]),
                                                        (double)Division(g_stStatRpt.unHet[1],g_stStatRpt.unHom[1]),
                                                        (double)Division(unTotalHet,unTotalHom));
    fprintf(g_pfOut,"\n");
    fprintf(g_pfOut,"Transitions SNP(Ts)\t%u\t%u\t%u\n",g_stStatRpt.unTs[0],g_stStatRpt.unTs[1],unTotalTs);  
    fprintf(g_pfOut,"Transversions SNP(Tv)\t%u\t%u\t%u\n",g_stStatRpt.unTv[0],g_stStatRpt.unTv[1],unTotalTv);
    fprintf(g_pfOut,"Ts/Tv Ratio\t%.2f\t%.2f\t%.2f\n",(double)Division(g_stStatRpt.unTs[0],g_stStatRpt.unTv[0]),
                                                      (double)Division(g_stStatRpt.unTs[1],g_stStatRpt.unTv[1]),
                                                      (double)Division(unTotalTs,unTotalTv)); 
    fprintf(g_pfOut,"\n");
    fprintf(g_pfOut,"Polymorphysm Type\tTs/Tv\tTotal SNP-Count\tTotal SNP-Percentage\tNovel SNP\t"
                    "Known SNP\tSynonymous SNP\tNon-Synonymous SNP\tSNP in Gene\tSNP in Exon\t"
                    "SNP in Intron\tSNP in Promoter\tSNP in miRNA\tSNP in UTR\n");

    fprintf(g_pfOut,"A>C\t%s\t%u\t%.2f%%\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\n",                     
                    (g_stStatRpt.AllelVar[A2C].unTsTv[0]!=0)?"Ts":"Tv",unA2C,dbA2C,
                    g_stStatRpt.AllelVar[A2C].undbSnp[1],g_stStatRpt.AllelVar[A2C].undbSnp[0],
                    g_stStatRpt.AllelVar[A2C].unSynoNonSyno[0],g_stStatRpt.AllelVar[A2C].unSynoNonSyno[1],
                    g_stStatRpt.AllelVar[A2C].unGene,g_stStatRpt.AllelVar[A2C].unExonIntron[0],
                    g_stStatRpt.AllelVar[A2C].unExonIntron[1],g_stStatRpt.AllelVar[A2C].unPromoter,
                    g_stStatRpt.AllelVar[A2C].unmiRNA,g_stStatRpt.AllelVar[A2C].unUTR); 

    fprintf(g_pfOut,"G>C\t%s\t%u\t%.2f%%\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\n",
                    (g_stStatRpt.AllelVar[G2C].unTsTv[0]!=0)?"Ts":"Tv",unG2C,dbG2C,
                    g_stStatRpt.AllelVar[G2C].undbSnp[1],g_stStatRpt.AllelVar[G2C].undbSnp[0],
                    g_stStatRpt.AllelVar[G2C].unSynoNonSyno[0],g_stStatRpt.AllelVar[G2C].unSynoNonSyno[1],
                    g_stStatRpt.AllelVar[G2C].unGene,g_stStatRpt.AllelVar[G2C].unExonIntron[0],
                    g_stStatRpt.AllelVar[G2C].unExonIntron[1],g_stStatRpt.AllelVar[G2C].unPromoter,
                    g_stStatRpt.AllelVar[G2C].unmiRNA,g_stStatRpt.AllelVar[G2C].unUTR);

    fprintf(g_pfOut,"T>C\t%s\t%u\t%.2f%%\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\n",
                    (g_stStatRpt.AllelVar[T2C].unTsTv[0]!=0)?"Ts":"Tv",unT2C,dbT2C,
                    g_stStatRpt.AllelVar[T2C].undbSnp[1],g_stStatRpt.AllelVar[T2C].undbSnp[0],
                    g_stStatRpt.AllelVar[T2C].unSynoNonSyno[0],g_stStatRpt.AllelVar[T2C].unSynoNonSyno[1],
                    g_stStatRpt.AllelVar[T2C].unGene,g_stStatRpt.AllelVar[T2C].unExonIntron[0],
                    g_stStatRpt.AllelVar[T2C].unExonIntron[1],g_stStatRpt.AllelVar[T2C].unPromoter,
                    g_stStatRpt.AllelVar[T2C].unmiRNA,g_stStatRpt.AllelVar[T2C].unUTR);  

    fprintf(g_pfOut,"A>G\t%s\t%u\t%.2f%%\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\n",
                    (g_stStatRpt.AllelVar[A2G].unTsTv[0]!=0)?"Ts":"Tv",unA2G,dbA2G,
                    g_stStatRpt.AllelVar[A2G].undbSnp[1],g_stStatRpt.AllelVar[A2G].undbSnp[0],
                    g_stStatRpt.AllelVar[A2G].unSynoNonSyno[0],g_stStatRpt.AllelVar[A2G].unSynoNonSyno[1],
                    g_stStatRpt.AllelVar[A2G].unGene,g_stStatRpt.AllelVar[A2G].unExonIntron[0],
                    g_stStatRpt.AllelVar[A2G].unExonIntron[1],g_stStatRpt.AllelVar[A2G].unPromoter,
                    g_stStatRpt.AllelVar[A2G].unmiRNA,g_stStatRpt.AllelVar[A2G].unUTR);
  
    fprintf(g_pfOut,"C>G\t%s\t%u\t%.2f%%\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\n",
                    (g_stStatRpt.AllelVar[C2G].unTsTv[0]!=0)?"Ts":"Tv",unC2G,dbC2G,
                    g_stStatRpt.AllelVar[C2G].undbSnp[1],g_stStatRpt.AllelVar[C2G].undbSnp[0],
                    g_stStatRpt.AllelVar[C2G].unSynoNonSyno[0],g_stStatRpt.AllelVar[C2G].unSynoNonSyno[1],
                    g_stStatRpt.AllelVar[C2G].unGene,g_stStatRpt.AllelVar[C2G].unExonIntron[0],
                    g_stStatRpt.AllelVar[C2G].unExonIntron[1],g_stStatRpt.AllelVar[C2G].unPromoter,
                    g_stStatRpt.AllelVar[C2G].unmiRNA,g_stStatRpt.AllelVar[C2G].unUTR
                    );

    fprintf(g_pfOut,"T>G\t%s\t%u\t%.2f%%\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\n",
                    (g_stStatRpt.AllelVar[T2G].unTsTv[0]!=0)?"Ts":"Tv",unT2G,dbT2G,
                    g_stStatRpt.AllelVar[T2G].undbSnp[1],g_stStatRpt.AllelVar[T2G].undbSnp[0],
                    g_stStatRpt.AllelVar[T2G].unSynoNonSyno[0],g_stStatRpt.AllelVar[T2G].unSynoNonSyno[1],
                    g_stStatRpt.AllelVar[T2G].unGene,g_stStatRpt.AllelVar[T2G].unExonIntron[0],
                    g_stStatRpt.AllelVar[T2G].unExonIntron[1],g_stStatRpt.AllelVar[T2G].unPromoter,
                    g_stStatRpt.AllelVar[T2G].unmiRNA,g_stStatRpt.AllelVar[T2G].unUTR);

    fprintf(g_pfOut,"C>A\t%s\t%u\t%.2f%%\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\n",
                    (g_stStatRpt.AllelVar[C2A].unTsTv[0]!=0)?"Ts":"Tv",unC2A,dbC2A, 
                    g_stStatRpt.AllelVar[C2A].undbSnp[1],g_stStatRpt.AllelVar[C2A].undbSnp[0],
                    g_stStatRpt.AllelVar[C2A].unSynoNonSyno[0],g_stStatRpt.AllelVar[C2A].unSynoNonSyno[1],
                    g_stStatRpt.AllelVar[C2A].unGene,g_stStatRpt.AllelVar[C2A].unExonIntron[0],
                    g_stStatRpt.AllelVar[C2A].unExonIntron[1],g_stStatRpt.AllelVar[C2A].unPromoter,
                    g_stStatRpt.AllelVar[C2A].unmiRNA,g_stStatRpt.AllelVar[C2A].unUTR 
                    );

    fprintf(g_pfOut,"G>A\t%s\t%u\t%.2f%%\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\n",
                    (g_stStatRpt.AllelVar[G2A].unTsTv[0]!=0)?"Ts":"Tv",unG2A,dbG2A,
                    g_stStatRpt.AllelVar[G2A].undbSnp[1],g_stStatRpt.AllelVar[G2A].undbSnp[0],
                    g_stStatRpt.AllelVar[G2A].unSynoNonSyno[0],g_stStatRpt.AllelVar[G2A].unSynoNonSyno[1],
                    g_stStatRpt.AllelVar[G2A].unGene,g_stStatRpt.AllelVar[G2A].unExonIntron[0],
                    g_stStatRpt.AllelVar[G2A].unExonIntron[1],g_stStatRpt.AllelVar[G2A].unPromoter,
                    g_stStatRpt.AllelVar[G2A].unmiRNA,g_stStatRpt.AllelVar[G2A].unUTR
                    ); 

    fprintf(g_pfOut,"T>A\t%s\t%u\t%.2f%%\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\n",
                    (g_stStatRpt.AllelVar[T2A].unTsTv[0]!=0)?"Ts":"Tv",unT2A,dbT2A, 
                    g_stStatRpt.AllelVar[T2A].undbSnp[1],g_stStatRpt.AllelVar[T2A].undbSnp[0],
                    g_stStatRpt.AllelVar[T2A].unSynoNonSyno[0],g_stStatRpt.AllelVar[T2A].unSynoNonSyno[1],
                    g_stStatRpt.AllelVar[T2A].unGene,g_stStatRpt.AllelVar[T2A].unExonIntron[0],
                    g_stStatRpt.AllelVar[T2A].unExonIntron[1],g_stStatRpt.AllelVar[T2A].unPromoter,
                    g_stStatRpt.AllelVar[T2A].unmiRNA,g_stStatRpt.AllelVar[T2A].unUTR
                    );  

    fprintf(g_pfOut,"A>T\t%s\t%u\t%.2f%%\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\n",
                    (g_stStatRpt.AllelVar[A2T].unTsTv[0]!=0)?"Ts":"Tv",unA2T,dbA2T,
                    g_stStatRpt.AllelVar[A2T].undbSnp[1],g_stStatRpt.AllelVar[A2T].undbSnp[0],
                    g_stStatRpt.AllelVar[A2T].unSynoNonSyno[0],g_stStatRpt.AllelVar[A2T].unSynoNonSyno[1],
                    g_stStatRpt.AllelVar[A2T].unGene,g_stStatRpt.AllelVar[A2T].unExonIntron[0],
                    g_stStatRpt.AllelVar[A2T].unExonIntron[1],g_stStatRpt.AllelVar[A2T].unPromoter,
                    g_stStatRpt.AllelVar[A2T].unmiRNA,g_stStatRpt.AllelVar[A2T].unUTR
                    ); 
    fprintf(g_pfOut,"C>T\t%s\t%u\t%.2f%%\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\n",
                    (g_stStatRpt.AllelVar[C2T].unTsTv[0]!=0)?"Ts":"Tv",unC2T,dbC2T,   
                    g_stStatRpt.AllelVar[C2T].undbSnp[1],g_stStatRpt.AllelVar[C2T].undbSnp[0],
                    g_stStatRpt.AllelVar[C2T].unSynoNonSyno[0],g_stStatRpt.AllelVar[C2T].unSynoNonSyno[1],
                    g_stStatRpt.AllelVar[C2T].unGene,g_stStatRpt.AllelVar[C2T].unExonIntron[0],
                    g_stStatRpt.AllelVar[C2T].unExonIntron[1],g_stStatRpt.AllelVar[C2T].unPromoter,
                    g_stStatRpt.AllelVar[C2T].unmiRNA,g_stStatRpt.AllelVar[C2T].unUTR
                    );
    fprintf(g_pfOut,"G>T\t%s\t%u\t%.2f%%\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\n",
                    (g_stStatRpt.AllelVar[G2T].unTsTv[0]!=0)?"Ts":"Tv",unG2T,
                    100.00-(dbA2C+dbG2C+dbT2C+dbA2G+dbC2G+dbT2G+dbC2A+dbG2A+dbT2A+dbA2T+dbC2T),  
                    g_stStatRpt.AllelVar[G2T].undbSnp[1],g_stStatRpt.AllelVar[G2T].undbSnp[0],
                    g_stStatRpt.AllelVar[G2T].unSynoNonSyno[0],g_stStatRpt.AllelVar[G2T].unSynoNonSyno[1],
                    g_stStatRpt.AllelVar[G2T].unGene,g_stStatRpt.AllelVar[G2T].unExonIntron[0],
                    g_stStatRpt.AllelVar[G2T].unExonIntron[1],g_stStatRpt.AllelVar[G2T].unPromoter,
                    g_stStatRpt.AllelVar[G2T].unmiRNA,g_stStatRpt.AllelVar[G2T].unUTR
                     );
    fprintf(g_pfOut,"\n");


    fprintf(g_pfOut,"Polymorphised Reference Allele\tSNP Count\tPercentage\n");   
    fprintf(g_pfOut,"A\t%u\t%.2f%%\n",g_stStatRpt.unRef[0],dbRefA);
    fprintf(g_pfOut,"C\t%u\t%.2f%%\n",g_stStatRpt.unRef[1],dbRefC);
    fprintf(g_pfOut,"G\t%u\t%.2f%%\n",g_stStatRpt.unRef[2],dbRefG);  
    fprintf(g_pfOut,"T\t%u\t%.2f%%\n",g_stStatRpt.unRef[3],100.00-(dbRefA+dbRefC+dbRefG));  
    fprintf(g_pfOut,"\n");

    fprintf(g_pfOut,"Variant Allele\tSNP Count\tPercentage\n");
    fprintf(g_pfOut,"A\t%u\t%.2f%%\n",g_stStatRpt.unVar[0],dbVarA);
    fprintf(g_pfOut,"C\t%u\t%.2f%%\n",g_stStatRpt.unVar[1],dbVarC);
    fprintf(g_pfOut,"G\t%u\t%.2f%%\n",g_stStatRpt.unVar[2],dbVarG);
    fprintf(g_pfOut,"T\t%u\t%.2f%%\n",g_stStatRpt.unVar[3],100.00-(dbVarA+dbVarC+dbVarG));
    fprintf(g_pfOut,"\n");

    fprintf(g_pfOut,"Confident Score Distribution\tSNP Count\tPercentage\tPercentage (cumulative)\n"); 

    fCumlPercent=fPercent=Division(g_stStatRpt.unScore[0],unTotalSNP)*100;    
    fprintf(g_pfOut,"0.0<Score<=0.1\t%u\t%.2f%%\t%.2f%%\n",g_stStatRpt.unScore[0],fPercent,fCumlPercent);

    for (short i=1;i<10;i++)
    {
        fPercent=Division(g_stStatRpt.unScore[i],unTotalSNP)*100;                  
        fCumlPercent+=fPercent;  
        fprintf(g_pfOut,"%.1f<Score<=%.1f\t%u\t%.2f%%\t%.2f%%\n",(float)i/10,(float)(i+1)/10,g_stStatRpt.unScore[i],
                        fPercent,fCumlPercent);
    }

    fprintf(g_pfOut,"\n");
}


int main(int argc, char** argv)
{
     if (argc < 4){
        fprintf(stdout,"Parameters not sufficient...\n\n");   goto Err_Exit;
     }

     g_pfIn = fopen(argv[1],"r"); if (!g_pfIn){ fprintf(stdout,"Failed to open %s ...\n",argv[1]); goto Err_Exit;}
     g_pfOut = fopen(argv[3],"w"); if (!g_pfOut) { fprintf(stdout,"Failed to open %s ...\n",argv[3]); goto Err_Exit;}
     g_nMinSR = atoi(argv[2]);

     for (short j=0; j<(g_cnSize); j++)
     {
         for (short k=0; k<(g_cnSize); k++) 
              g_unNSCnt[j][k]=0;  
     } 

     ProcessInputFile(); OutputStats(); OutputProteinVarMatrix();

     CloseFiles();
     return (EXIT_SUCCESS);

Err_Exit:
     CloseFiles(); banner(argv); exit(9);
}
