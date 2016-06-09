/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package historicallinguistcs;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Random;
import java.util.Scanner;

public class HistoricalLinguistcs {

    double[] equlibriumProb;
    HashMap<Character,Integer> charPosition = new HashMap<Character,Integer>();
    HashMap<Integer,Character> numberChar = new HashMap<Integer,Character>();
    double[][] rateMatrix;
    double lambda = 0.1, mu=0.11,beta;
    public static void main(String[] args)
    {    
      HistoricalLinguistcs object1 = new HistoricalLinguistcs();
      object1.equlibriumProb = object1.getEquilibriumProb();
      object1.beta = (1 - Math.exp(object1.lambda - object1.mu))/(object1.mu - object1.lambda*Math.exp(object1.lambda- object1.mu));
      object1.rateMatrix = object1.getRateMatrix();
      object1.charPosition = object1.getCharPosition();
      object1.numberChar = object1.getNumberChar();
      System.out.println("lamba "+ object1.lambda + " mu "+ object1.mu);
     // while(true)
      //{
      Scanner in = new Scanner(System.in);
      String string1,string2;
      System.out.println("Sequence 1");
      string1 = in.nextLine();
     // System.out.println("Sequence 2");
     // string2 = in.nextLine();
      Double Ps1s2=0.0, Ps1=0.0,Ps2=0.0;
      Ps1 = object1.eqProbString(string1);
       
      ArrayList<String> simulatedSet = new ArrayList<String>();// = object1.simulateSet(string1);
      simulatedSet.add("AGACGAGTAGC");
      int i=0;
      while(i<simulatedSet.size())
      {
          string2 = simulatedSet.get(i);
          System.out.println(string2);
          Ps2 = object1.eqProbString(string2);
          Ps1s2 = object1.probCalculation(string1, string2) * Ps1;
          System.out.println("Probability P(s1,s2) is "+ Ps1s2);
          double stat;
          stat =-2 * Math.log(Ps1s2/(Ps1*Ps2));
          System.out.println("Statistic for homology testing U is "+ stat);
          System.out.println();
          i++;
      }
      
         
        
//      
//      Ps1s2 = object1.mostProbCalculation(string1, string2) * Ps1;
//      System.out.println("Probability P(s1,s2) for most probable alignment is "+ Ps1s2);
//      stat =-2 * Math.log(Ps1s2/(Ps1*Ps2));
//      System.out.println("Statistic for homology testing U is "+ stat);
     // }
    }
    public List<String> simulateSet(String root)
    {
        ArrayList<String> simulatedString = new ArrayList<String>();
        int i=0;
        String temp = root;
        while(i<20)
        {
            temp = mutateString(root);
            temp = performInsertion(temp);
            temp = performDeletion(temp);
            simulatedString.add(temp);
            i++;
        }
        return simulatedString;
    }
    
    public String mutateString(String root)
    {
        int i=0;
        int j=0;
        char temp;
        Random rand = new Random();
        StringBuffer newString = new StringBuffer(root);
        while(i<root.length())
        {
            j = rand.nextInt(4);
            temp = numberChar.get(j);
            newString.setCharAt(i, temp); 
            i++;
        }
        return newString.toString();
    }
    
    public String performInsertion(String root)
    {
        Random rand = new Random();
        int ins = rand.nextInt(30);
        StringBuffer strng = new StringBuffer(root);
        strng.setLength(ins + root.length());
        int i =0,j=0,k=0;
        char temp;
        while(i<ins)
        {
            j = rand.nextInt(4);
            temp = numberChar.get(j);
          
            k = rand.nextInt(strng.length());
            // System.out.println(i + " "+ ins + "  "+ strng.length() + " " +k);
            strng.insert(k,temp);
            i++;
        }
        return strng.toString();
    }
    
    public String performDeletion(String root)
    {
         Random rand = new Random();
        int del = rand.nextInt(Math.min(30,root.length()));
        StringBuffer strng = new StringBuffer(root);
        int i =0,k=0;
        while(i<del)
        {
            k = rand.nextInt(strng.length());          
            strng.deleteCharAt(k);
            i++;
        }
              
        return strng.toString().trim().replaceAll("[^a-zA-Z]", "");
      
    }
    
    public double[][] getRateMatrix()
    {
        double[][] rateMatrix = {{-0.3,0.1,0.1,0.1},{0.1,-0.3,0.1,0.1},{0.1,0.1,-0.3,0.1},{0.1,0.1,0.1,-0.3}};
        return rateMatrix;
    }
    
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
    
    public HashMap<Integer,Character> getNumberChar()
    {
        HashMap<Integer,Character> charPosition = new HashMap<Integer,Character>();
        charPosition.put(0,'A');
        charPosition.put(1,'C');
        charPosition.put(2,'G');
        charPosition.put(3,'T');
        return charPosition;
    }
    
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
    
     public double mostProbCalculation(String string1, String string2) 
    {  
        //Base case
       if(string1.length()==0)
       {
         //  System.out.println(string1);
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
               if(temp1>temp2)
               {
                   //printString = string1.substring(0,i-1) + string2.charAt(j-k) + string2.substring(Math.min(j+1, j-k+1), j+1);
               }
               else
               {
                //  printString = string1.substring(0,i-1) + string2.substring( j-k, j+1);

               }
           }
       }
       /*
       if(term1>term2)
       {
           System.out.println(string1.substring(0,i-1));
       }
       else
       {
           
           System.out.println(printString);
       }*/
       return Math.max(term1, term2);   
    }
    
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
    
    public double getSubProb(char a, char b)
    {
        int pos1 = charPosition.get(a);
        int pos2 = charPosition.get(b);        
        return Math.exp(rateMatrix[pos1][pos2]);
    }
    
    public double mortalDescendant(int k)
    {
        if(k==0)
            return 0;
        
        return Math.exp(mu)*(1 - lambda*beta)*Math.pow(lambda*beta, k-1);
    }
    
     public double immortalDescendant(int k)
    {   
        if(k==0)
            return 0;
        
        return (1 - lambda*beta)*Math.pow(lambda*beta, k-1);
    }
     
      public double noMortalDescendant(int k)
    {
        if(k==0)
            return mu*beta;
        
        return (1 - Math.exp(-mu) - mu*beta)*(1 - lambda*beta)*Math.pow(lambda*beta, k-1);
    }
}
