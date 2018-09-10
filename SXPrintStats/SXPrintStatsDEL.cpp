#include <stdio.h>
#include <stdlib.h>
#include <string.h>

enum enmAllele{A=0,C,G,T};

const short enmAllele_Size = T+1;

typedef struct _stAllele
{
    unsigned int undbSnp[2];         
    unsigned int unGene; 
    unsigned int unExonIntron[2];
    unsigned int unPromoter;
    unsigned int unmiRNA;
    unsigned int unUTR;     

    _stAllele()
    {
        unGene=unPromoter=unmiRNA=unUTR=0;

        for (short i=0; i<2; i++)  
        {  
            undbSnp[i]=unExonIntron[i]=0; 
        }
    }
}stAllele;


typedef struct _stStatRpt
{
    unsigned int undbSnp[2];
    unsigned int unGene[2]; 
    unsigned int unExon[2];
    unsigned int unPromoter[2];
    unsigned int unmiRNA[2];
    unsigned int unUTR[2];
    unsigned int unHet[2];
    unsigned int unHom[2];   
    unsigned int unGreyArea[2];
    stAllele Allele[enmAllele_Size];           
    
    unsigned int unScore[10]; 

    _stStatRpt()
    {
         for(int i=0;i<2;i++) {
             undbSnp[i]=unGene[i]=unExon[i]=unPromoter[i]=unmiRNA[i]=unUTR[i]=
             unHet[i]=unHom[i]=unGreyArea[i]=0;
         } 

         for (int i=0; i<10; i++){ unScore[i]=0; }
    }

} stStatRpt;

stStatRpt g_stStatRpt;
enmAllele g_eAllele;
int g_nMinSR=1;
FILE *g_pfIn=NULL, *g_pfOut=NULL;

static void banner(char *argv[])
{
    printf("Synamatix SXPrintStatsDEL Copyright 2011 Synamatix Sdn Bhd\n");
    printf("Built %s %s\n", __DATE__, __TIME__);
    printf("Usage: %s\n", argv[0]);
    printf("\t<Input filename> <Min. Supporting Read> <Output filename>\n\n");
}


void CloseFiles()
{    if (g_pfIn) fclose(g_pfIn); if (g_pfOut) fclose(g_pfOut);
}


void ProcessInputFile()
{
     fprintf(stdout,"ProcessInputFile...\n");

     unsigned long ulTest1=0, ulTest2=0;
     char cAllele,acbuf[40960]; char *pChr=NULL;
     short nIdx=0; 
     bool bdbSnp;  

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
         pChr = strtok(NULL,"\t"); //Deleted_Base	
         cAllele=pChr[0];
         
         pChr = strtok(NULL,"\t"); //Fwd_DEL_Reads	
         pChr = strtok(NULL,"\t"); //Rvs_DEL_Reads	
         pChr = strtok(NULL,"\t"); //Total_DEL_Reads	
         if (atoi(pChr) < g_nMinSR) continue;

         switch (cAllele){
             case 'A': { g_eAllele=A; break; }
             case 'C': { g_eAllele=C; break; }
             case 'G': { g_eAllele=G; break; }
             case 'T': { g_eAllele=T; break; }
         }

         pChr = strtok(NULL,"\t"); //Total_Read_Density	
         pChr = strtok(NULL,"\t"); //PEnd_Count	
         pChr = strtok(NULL,"\t"); //Avg_QryPos	
         pChr = strtok(NULL,"\t"); //Avg_QScore	

         pChr = strtok(NULL,"\t"); //dbSNP	

         bdbSnp = (pChr[0]!='-')? true:false;

         g_stStatRpt.undbSnp[!bdbSnp]++;
         g_stStatRpt.Allele[g_eAllele].undbSnp[!bdbSnp]++;
                         
         pChr = strtok(NULL,"\t"); //Gene_Name
         if (pChr[0]!='-') 
         { 
             if (bdbSnp) g_stStatRpt.unGene[0]++; else g_stStatRpt.unGene[1]++; 
             g_stStatRpt.Allele[g_eAllele].unGene++;
         }

         pChr = strtok(NULL,"\t"); //Gene_Description	
         pChr = strtok(NULL,"\t"); //Keyword	
         pChr = strtok(NULL,"\t"); //miRNA	

         if (pChr[0]!='-') 
         {             
            g_stStatRpt.unmiRNA[!bdbSnp]++; g_stStatRpt.Allele[g_eAllele].unmiRNA++;
         }

         pChr = strtok(NULL,"\t"); //Promoter	
         if (pChr[0]=='Y') 
         { 
            g_stStatRpt.unPromoter[!bdbSnp]++; g_stStatRpt.Allele[g_eAllele].unPromoter++;  
         }

         pChr = strtok(NULL,"\t"); //UTR	
         if (pChr[0]=='Y') 
         {
            g_stStatRpt.unUTR[!bdbSnp]++; g_stStatRpt.Allele[g_eAllele].unUTR++; 
         }

         pChr = strtok(NULL,"\t"); //Exon	
         if (pChr[0]!='-') 
         {     
            g_stStatRpt.unExon[!bdbSnp]++; g_stStatRpt.Allele[g_eAllele].unExonIntron[0]++;
         }
         else
         {
             g_stStatRpt.Allele[g_eAllele].unExonIntron[1]++; 
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
}



inline double Division(unsigned long ulVal, unsigned long ulDenominator)
{
    if (ulDenominator == 0) return 0;
    return (ulVal == 0)? 0:(double)ulVal/ulDenominator;
}

void OutputStats()
{
    unsigned int unTotalINS=g_stStatRpt.undbSnp[0]+g_stStatRpt.undbSnp[1];
    unsigned int unTotalGene=g_stStatRpt.unGene[0]+g_stStatRpt.unGene[1];
    unsigned int unTotalExon=g_stStatRpt.unExon[0]+g_stStatRpt.unExon[1];
    unsigned int unTotalHet=g_stStatRpt.unHet[0]+g_stStatRpt.unHet[1];
    unsigned int unTotalHom=g_stStatRpt.unHom[0]+g_stStatRpt.unHom[1];
    
    double dbKnownSNP=(double)Division(g_stStatRpt.undbSnp[0],unTotalINS)*100;

    unsigned int unA=g_stStatRpt.Allele[A].undbSnp[0]+g_stStatRpt.Allele[A].undbSnp[1];
    unsigned int unC=g_stStatRpt.Allele[C].undbSnp[0]+g_stStatRpt.Allele[C].undbSnp[1];
    unsigned int unG=g_stStatRpt.Allele[G].undbSnp[0]+g_stStatRpt.Allele[G].undbSnp[1];;
    unsigned int unT=g_stStatRpt.Allele[T].undbSnp[0]+g_stStatRpt.Allele[T].undbSnp[1];

    double dbA=(double)Division(unA,unTotalINS)*100;
    double dbC=(double)Division(unC,unTotalINS)*100;
    double dbG=(double)Division(unG,unTotalINS)*100;

    float fPercent=0.0, fCumlPercent=0.0;  

    //fprintf(g_pfOut,"\nMin. Supporting Reads = %d\n",g_nMinSR);

    fprintf(g_pfOut,"\n\nDEL Summary\tKnown\tNovel\tTotal\n");
    fprintf(g_pfOut,"Total DEL (Count)\t%u\t%u\t%u\n",g_stStatRpt.undbSnp[0],g_stStatRpt.undbSnp[1],unTotalINS);                    
    fprintf(g_pfOut,"Total DEL (by %%)\t%.2f%%\t%.2f%%\t100.00%%\n",dbKnownSNP,100.00-dbKnownSNP);

    fprintf(g_pfOut,"DEL in Gene (Count)\t%u\t%u\t%u\n",g_stStatRpt.unGene[0],g_stStatRpt.unGene[1],unTotalGene);
    fprintf(g_pfOut,"DEL in Gene (By %%)\t%.2f%%\t%.2f%%\t%.2f%%\n",Division(g_stStatRpt.unGene[0],unTotalINS)*100,
                                                             Division(g_stStatRpt.unGene[1],unTotalINS)*100,
                                                             Division(unTotalGene,unTotalINS)*100);
                                                               
    fprintf(g_pfOut,"DEL in Exon (Count)\t%u\t%u\t%u\n",g_stStatRpt.unExon[0],g_stStatRpt.unExon[1],unTotalExon);
    fprintf(g_pfOut,"DEL in Exon (By %%)\t%.2f%%\t%.2f%%\t%.2f%%\n",Division(g_stStatRpt.unExon[0],unTotalINS)*100,
                                                             Division(g_stStatRpt.unExon[1],unTotalINS)*100,
                                                             Division(unTotalExon,unTotalINS)*100); 

    fprintf(g_pfOut,"DEL in Promoter\t%u\t%u\t%u\n",g_stStatRpt.unPromoter[0],g_stStatRpt.unPromoter[1],
                                                    g_stStatRpt.unPromoter[0]+g_stStatRpt.unPromoter[1]);
    fprintf(g_pfOut,"DEL in miRNA\t%u\t%u\t%u\n",g_stStatRpt.unmiRNA[0],g_stStatRpt.unmiRNA[1],
                                                 g_stStatRpt.unmiRNA[0]+g_stStatRpt.unmiRNA[1]);
    fprintf(g_pfOut,"DEL in UTR\t%u\t%u\t%u\n",g_stStatRpt.unUTR[0],g_stStatRpt.unUTR[1],
                                               g_stStatRpt.unUTR[0]+g_stStatRpt.unUTR[1]);
    fprintf(g_pfOut,"\n");
    fprintf(g_pfOut,"Heterozygous DEL(het)\t%u\t%u\t%u\n",g_stStatRpt.unHet[0],g_stStatRpt.unHet[1],unTotalHet); 
    fprintf(g_pfOut,"Homozygous DEL(hom)\t%u\t%u\t%u\n",g_stStatRpt.unHom[0],g_stStatRpt.unHom[1],unTotalHom); 
    fprintf(g_pfOut,"Grey Area\t%u\t%u\t%u\n",g_stStatRpt.unGreyArea[0],g_stStatRpt.unGreyArea[1],
                                              g_stStatRpt.unGreyArea[0]+g_stStatRpt.unGreyArea[1]);
    fprintf(g_pfOut,"Het/Hom Ratio\t%.2f\t%.2f\t%.2f\n",(double)Division(g_stStatRpt.unHet[0],g_stStatRpt.unHom[0]),
                                                        (double)Division(g_stStatRpt.unHet[1],g_stStatRpt.unHom[1]),
                                                        (double)Division(unTotalHet,unTotalHom));
    fprintf(g_pfOut,"\n");
    fprintf(g_pfOut,"Deleted Bases\tTotal DEL-Count\tTotal DEL-Percentage\tNovel DEL\t"
                    "Known DEL\tDEL in Gene\tDEL in Exon\t"
                    "DEL in Intron\tDEL in Promoter\tDEL in miRNA\tDEL in UTR\n");

    fprintf(g_pfOut,"A\t%u\t%.2f%%\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\n",                     
                    unA,dbA,g_stStatRpt.Allele[A].undbSnp[1],g_stStatRpt.Allele[A].undbSnp[0],
                    g_stStatRpt.Allele[A].unGene,g_stStatRpt.Allele[A].unExonIntron[0],
                    g_stStatRpt.Allele[A].unExonIntron[1],g_stStatRpt.Allele[A].unPromoter,
                    g_stStatRpt.Allele[A].unmiRNA,g_stStatRpt.Allele[A].unUTR); 

    fprintf(g_pfOut,"C\t%u\t%.2f%%\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\n",
                    unC,dbC,g_stStatRpt.Allele[C].undbSnp[1],g_stStatRpt.Allele[C].undbSnp[0],
                    g_stStatRpt.Allele[C].unGene,g_stStatRpt.Allele[C].unExonIntron[0],
                    g_stStatRpt.Allele[C].unExonIntron[1],g_stStatRpt.Allele[C].unPromoter,
                    g_stStatRpt.Allele[C].unmiRNA,g_stStatRpt.Allele[C].unUTR);

    fprintf(g_pfOut,"G\t%u\t%.2f%%\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\n",
                    unG,dbG,g_stStatRpt.Allele[G].undbSnp[1],g_stStatRpt.Allele[G].undbSnp[0],
                    g_stStatRpt.Allele[G].unGene,g_stStatRpt.Allele[G].unExonIntron[0],
                    g_stStatRpt.Allele[G].unExonIntron[1],g_stStatRpt.Allele[G].unPromoter,
                    g_stStatRpt.Allele[G].unmiRNA,g_stStatRpt.Allele[G].unUTR);
      
    fprintf(g_pfOut,"T\t%u\t%.2f%%\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\n",
                    unT,100.00-(dbA+dbC+dbG),
                    g_stStatRpt.Allele[T].undbSnp[1],g_stStatRpt.Allele[T].undbSnp[0],
                    g_stStatRpt.Allele[T].unGene,g_stStatRpt.Allele[T].unExonIntron[0],
                    g_stStatRpt.Allele[T].unExonIntron[1],g_stStatRpt.Allele[T].unPromoter,
                    g_stStatRpt.Allele[T].unmiRNA,g_stStatRpt.Allele[T].unUTR);
  
    fprintf(g_pfOut,"\n");

    fprintf(g_pfOut,"Confident Score Distribution\tDEL Count\tPercentage\tPercentage (cumulative)\n"); 

    fCumlPercent=fPercent=Division(g_stStatRpt.unScore[0],unTotalINS)*100;    
    fprintf(g_pfOut,"0.0<Score<=0.1\t%u\t%.2f%%\t%.2f%%\n",g_stStatRpt.unScore[0],fPercent,fCumlPercent);

    for (short i=1;i<10;i++)
    {
        fPercent=Division(g_stStatRpt.unScore[i],unTotalINS)*100;                  
        fCumlPercent+=fPercent;  
        fprintf(g_pfOut,"%.1f<Score<=%.1f\t%u\t%.2f%%\t%.2f%%\n",(float)i/10,(float)(i+1)/10,g_stStatRpt.unScore[i],
                        fPercent,fCumlPercent);
    }

    fprintf(g_pfOut,"\n");
}


int main(int argc, char** argv)
{
     if (argc < 3){
        fprintf(stdout,"Parameters not sufficient...\n\n");   goto Err_Exit;
     }

     g_pfIn = fopen(argv[1],"r"); if (!g_pfIn){ fprintf(stdout,"Failed to open %s ...\n",argv[1]); goto Err_Exit;}
     g_pfOut = fopen(argv[3],"w"); if (!g_pfOut) { fprintf(stdout,"Failed to open %s ...\n",argv[3]); goto Err_Exit;}
     g_nMinSR = atoi(argv[2]);

     ProcessInputFile(); OutputStats(); 

     CloseFiles();
     return (EXIT_SUCCESS);

Err_Exit:
     CloseFiles(); banner(argv); exit(9);
}
