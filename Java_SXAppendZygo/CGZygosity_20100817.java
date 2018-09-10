
import java.io.*;
import java.math.BigInteger;
import java.math.BigDecimal;
//import java.lang.Object;
import java.util.Formatter;



public class CGZygosity {

	static int[] InvertedReads = null;;
	static BufferedWriter dbsnplist = null;
	static int DiploidProbabilityParam = 5;

	   
	public static void main(String arg[]) throws EOFException{		
		
		if (arg.length < 2){
			System.out.print("Parameters not sufficient...\n\n");
			System.out.format("SXAppendZygo\n");
			System.out.format("\t<Input filename> <Output filename>\n\n");
		    return;
		}

		File inputfile = null; String strOut="";   
		int delCount, TotReadCount, homoCount = 0, heteroCount = 0, allReadCount = 0;
        double RS, Haploid1 = 0, Haploid2 = 0, Haploid3 = 0, Haploid = 0, Diploid = 0;
        String DeletionData = ""; String[] DeletionDataSplit;
		
		try{
			inputfile = new File(arg[0]);
			FileInputStream file_input = new FileInputStream (inputfile);
			DataInputStream data_in = new DataInputStream(file_input);			
			FileOutputStream file_output = new FileOutputStream(arg[1]);				

			while (true){
				allReadCount++; DeletionData =  data_in.readLine();                                  
                
                if(DeletionData==null) break;
                
                DeletionDataSplit = DeletionData.split("\t");
                                
                delCount = Integer.parseInt(DeletionDataSplit[4]);
                TotReadCount = Integer.parseInt(DeletionDataSplit[7]);
                RS = (double) delCount / (double) TotReadCount;                
                    
                if (delCount < 3){
                	strOut = DeletionData + "\t-\n";
                }
                else if(!(TotReadCount > 275))
                {                
	                BigInteger HeteroNew1 = Factorial4.factorial(delCount);   
	                BigInteger HeteroNew2 = Factorial4.factorial(TotReadCount);   
	                BigInteger HeteroNew3 = Factorial4.factorial(TotReadCount-delCount);  
	                BigDecimal HeteroNew4 = new BigDecimal(DiploidProbabilityParam);	
	                BigDecimal SquaredTerm = new BigDecimal("0");	
	                SquaredTerm = HeteroNew4.pow(TotReadCount);
	
	                BigDecimal a = new BigDecimal(HeteroNew1);
	                BigDecimal b = new BigDecimal(HeteroNew3);
	                BigDecimal c = new BigDecimal(HeteroNew2);
	                BigDecimal d = new BigDecimal("0");
	                BigDecimal e = new BigDecimal("0");
	                BigDecimal f = SquaredTerm;
	                BigDecimal g = new BigDecimal("0");
	                BigDecimal h = new BigDecimal("0");
	                d = a.multiply(b);
	                e = c.divide(d);
	                
	                g = e.divide(f, 200, BigDecimal.ROUND_FLOOR);
	                double hhh = g.doubleValue();                            
	                
	                if (TotReadCount == 2*delCount)
	                	Diploid = (-Math.log10(hhh) -0.301);
	                else
	                	Diploid = -Math.log10(hhh);
	                
	                BigDecimal k = new BigDecimal("0.99");
	                BigDecimal m = new BigDecimal("0.01");
	              
	                double h1 = e.doubleValue();
	                Haploid1 =  (-Math.log10(h1));
	                double h2 = k.doubleValue();
	                Haploid2 =  (-Math.log10(h2)) * (double) delCount;                        
	                double h3 = m.doubleValue();
	                Haploid3 =  (-Math.log10(h3)) * (double) (TotReadCount-delCount);                           
	                Haploid = Haploid1 + Haploid2 + Haploid3;
	                	                	                
	                if(Haploid<Diploid && !Double.isInfinite(Haploid) && !Double.isInfinite(Diploid)){
	                    homoCount++; strOut = DeletionData + "\thom\n"; 	                    	                   
	                }
	                else if(Diploid<Haploid && !Double.isInfinite(Haploid) && !Double.isInfinite(Diploid)){
	                    heteroCount++; strOut = DeletionData + "\thet\n";  	                    
	                }                                
	                else {
	                	strOut = DeletionData + "\t-\n"; 
	                }
                }
                else if(TotReadCount > 275 && delCount>2)
                {                    
                    if(RS<0.7){
                        heteroCount++; strOut = DeletionData + "\thet\n"; 
                    }
                    else{
                        homoCount++; strOut = DeletionData + "\thom\n"; 
                    }
                }			       
                
                //file_output.write(strOut.getBytes());
			}	
			
			file_output.close();
			file_input.close();
		}
		catch (IOException e) {
		      System.out.println(e);
		}
		catch(Exception e)
		{
			e.printStackTrace();
		}
		
		
		System.out.println("Done");
		System.out.println(allReadCount-1);   
        System.out.println("homoCount : " + homoCount);  
        System.out.println("heteroCount : " + heteroCount);  
	}
}
