/*
 TKF91 model and homology testing.
 Paper Implemented: Statistical Alignment: Computational Properties, Homology Testing and Goodness-of-Fit
 J. Hein, C. Wiuf2, B. Knudsen1, M. B. Mùller3 and G. Wibling3   (Hein2000)
 */
 
package TKF91;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Random;
import java.util.Scanner;

public class TKF91 {

    double[] equlibriumProb;
    HashMap<Character,Integer> charPosition = new HashMap<Character,Integer>();
    HashMap<Integer,Character> numberChar = new HashMap<Integer,Character>();
    double[][] rateMatrix;
    double lambda = 0.1, mu=0.11,beta;
    public static void main(String[] args)
    {    
      TKF91 object1 = new TKF91();
      
      //Initializing all the parameters
      object1.equlibriumProb = object1.getEquilibriumProb();
      object1.beta = (1 - Math.exp(object1.lambda - object1.mu))/(object1.mu - object1.lambda*Math.exp(object1.lambda- object1.mu));
      object1.rateMatrix = object1.getRateMatrix();
      object1.charPosition = object1.getCharPosition();
      object1.numberChar = object1.getNumberChar();
      System.out.println("lamba "+ object1.lambda + " mu "+ object1.mu);
    
      // Sequence input from the user  
      Scanner in = new Scanner(System.in);
      String string1,string2;
      System.out.println("Sequence 1");
      string1 = in.nextLine();
      System.out.println("Sequence 2");
      string2 = in.nextLine();
      
      // Calculating the probabilites P(s1), P(s2) and P(s1,s2)
      Double Ps1s2=0.0, Ps1=0.0,Ps2=0.0;
      Ps1 = object1.eqProbString(string1);
      Ps2 = object1.eqProbString(string2);
      Ps1s2 = object1.probCalculation(string1, string2) * Ps1;
      System.out.println("Probability P(s1,s2) is "+ Ps1s2);
      
      //Statistics for homology testing -2*ln(P(s1,s2)/P(s1)P(s2))
      double stat;
      stat =-2 * Math.log(Ps1s2/(Ps1*Ps2));
      System.out.println("Statistic for homology testing U is "+ stat);
      System.out.println();
    }
    
   
    //@getRateMatrix function returns the rate matrix for the alphabets
    public double[][] getRateMatrix()
    {
        double[][] rateMatrix = {{-0.3,0.1,0.1,0.1},{0.1,-0.3,0.1,0.1},{0.1,0.1,-0.3,0.1},{0.1,0.1,0.1,-0.3}};
        return rateMatrix;
    }
    //@getEquilibriumProb function return the equilibrium probability for the alphabets
    public double[] getEquilibriumProb()
    {
        double[] equlibriumProb = {0.25,0.25,0.25,0.25};
        return equlibriumProb;       
    }
    
    public HashMap<Character,Integer> getCharPosition()
    {
        HashMap<Character,Integer> charPosition = new HashMap<Character,Integer>();
        charPosition.put('A', 0);
        charPosition.put('C', 1);
        charPosition.put('G', 2);
        charPosition.put('T', 3);
        return charPosition;
    }
    
    
    //@probCalculation function returns the P(a|b) where a and b are sequences by summing over all probability of all possible statiscal alignment paths
    // Model implemented from the Hein2000 paper
    public double probCalculation(String string1, String string2) 
    {  
        //Base case
       if(string1.length()==0)
       {
           double eqProb= eqProbString(string2);
           return eqProb*immortalDescendant(string2.length());
       }
       
       int i= string1.length()-1;
       int j= string2.length()-1; 
       double term1 = noMortalDescendant(0)*probCalculation(string1.substring(0,i),string2);
       
       double term2  =0;
       double temp,temp1,temp2;
       for(int k=0; k <=j; k++)
       {
           temp1 = mortalDescendant(k+1)*getSubProb(string1.charAt(i),string2.charAt(j-k))*eqProbString(string2.substring(Math.min(j+1,j-k+1),j+1)) ;
           temp2 =noMortalDescendant(k+1)*eqProbString(string2.substring(j-k,j+1));
           temp = temp1 + temp2;
           term2 = term2 + probCalculation(string1.substring(0,i),string2.substring(0,j-k)) * temp;
       }
       return term1 + term2;   
    }
    
    //@mostProbCalculation returns the probabilty of most probable alignment path between the two sequences
     public double mostProbCalculation(String string1, String string2) 
    {  
        //Base case
       if(string1.length()==0)
       {
           double eqProb= eqProbString(string2);
           return eqProb*immortalDescendant(string2.length());
       }
       
       int i= string1.length()-1;
       int j= string2.length()-1; 
       double term1 = noMortalDescendant(0)*probCalculation(string1.substring(0,i),string2);
       
       double term2  =0;
       double temp,temp1,temp2;
       String printString="";
       int storeK=0;
       for(int k=0; k <=j; k++)
       {
           temp1 = mortalDescendant(k+1)*getSubProb(string1.charAt(i),string2.charAt(j-k))*eqProbString(string2.substring(Math.min(j+1,j-k+1),j+1)) ;
           temp2 =noMortalDescendant(k+1)*eqProbString(string2.substring(j-k,j+1));
           temp = Math.max(temp1, temp2);
           if(probCalculation(string1.substring(0,i),string2.substring(0,j-k)) * temp>term2)
           {
               term2 = probCalculation(string1.substring(0,i),string2.substring(0,j-k)) * temp;
               storeK=k;
           }
       return Math.max(term1, term2);   
    }
    
    //@eqProbString function returns the probabilty of a given sequence
    public double eqProbString(String s1)
    {
         double eqProb=1;
         int j = s1.length();
         int pos;
         for(int t=0; t<j; t++)
           {
               pos = charPosition.get(s1.charAt(t));
               eqProb = eqProb*equlibriumProb[pos]*lambda/mu;
           }
         return eqProb*(1-lambda/mu);
    }
    
    //@getsubProb returns the substitution probability of alphabets given the rate matrix
    public double getSubProb(char a, char b)
    {
        int pos1 = charPosition.get(a);
        int pos2 = charPosition.get(b);        
        return Math.exp(rateMatrix[pos1][pos2]);
    }
   //@mortalDescendant function finds the probability of mortal link surviving and leaving k descendant behind 
    public double mortalDescendant(int k)
    {
        if(k==0)
            return 0;
        
        return Math.exp(mu)*(1 - lambda*beta)*Math.pow(lambda*beta, k-1);
    }

//@mortalDescendant function finds the probability of immortal link surviving and leaving k descendant behind
     public double immortalDescendant(int k)
    {   
        if(k==0)
            return 0;
        
        return (1 - lambda*beta)*Math.pow(lambda*beta, k-1);
    }
    
     //@mortalDescendant function finds the probability of mortal link not surviving and leaving k descendant behind 
      public double noMortalDescendant(int k)
    {
        if(k==0)
            return mu*beta;
        
        return (1 - Math.exp(-mu) - mu*beta)*(1 - lambda*beta)*Math.pow(lambda*beta, k-1);
    }
}
