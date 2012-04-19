package bioGUI.model;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.lang.StringBuilder;

import java.lang.Math; 
import java.util.ArrayList;
import java.util.Scanner;


public class CGContent {
   public static File inputFile;

	public static void main(String[] args) throws NumberFormatException, IOException {
			int i = 0; 
			
			if(args[0].contains("GFF"))
			{
				GFF[] files = readGFF(args[0]);
				
				for(GFF gene : files)
				{
					System.out.println(gene.filename + " " + gene.start + " " + gene.end);
					
					gcFrameCount(String.valueOf(i), gene.filename + ".txt", gene.start, gene.end, Integer.valueOf(args[1]), Integer.valueOf(args[2]));
				}
			}
			else gcFrameCount(args[0].substring(0, args[0].indexOf('.')), args[0], Integer.valueOf(args[1]), Integer.valueOf(args[2]), Integer.valueOf(args[3]), Integer.valueOf(args[4]));
			
			
			
	}

	public static String readGene(String filename)
	{
		StringBuilder fileContents = new StringBuilder();
		
		File gene;
		
		try
		{
			gene = new File(filename);
			Scanner sc = new Scanner(gene);
			
			sc.nextLine(); //ignore first line
			
			while (sc.hasNextLine()) {
				fileContents.append(sc.nextLine());
			}
			sc.close(); 
		}
		catch(Exception ex)
		{
			System.out.println(ex + " - Invalid filename?");
		}
		
		return fileContents.toString();
	}
	
	public static GFF[] readGFF(String filename)
	{	
		ArrayList<GFF> strandsAL = new ArrayList<GFF>();  
		
		GFF[] strands = new GFF[1]; 
		
		try
		{
			File gene = new File(filename);
			Scanner sc = new Scanner(gene);

			sc.nextLine(); //ignore first line: browser position
			sc.nextLine(); // track name, description, visibility, color
			
			while (sc.hasNextLine()) 
			{
				GFF current = new GFF(); 
				
				current.filename = sc.next();
				sc.next(); // GEP
				sc.next(); //mRNA or CDS
				current.start = sc.nextInt(); 
				current.end = sc.nextInt(); 
				sc.next(); //dot
				current.direction = (sc.next() == "+"); 
				sc.nextLine(); //dot
				
				strandsAL.add(current); 
				
			}
			sc.close(); 
		}
		catch(Exception ex)
		{
			System.out.println(ex + " - Invalid GFF filename?");
		}
		
		return strandsAL.toArray(strands); 
	}
	
	/**
	 * gcFramecount
	 * Outputs a CSV histogram of gc percentages over a genome
	 * @throws IOException 
	 */
	public static double[] gcFrameCount(String CSVName, String filename, int start, int end, int frame, int shift) throws IOException
	{
		long total= 0;
		int sliceCount;

		String sub = readGene(filename); 
				
		end = Math.min(sub.length(), end); 		
		
		sub = sub.substring(start, end);
		
		int length = sub.length(); 
		
		//fine histogram of shift-width percentages
	
		
		if(shift != 0)
			sliceCount = (int)(.5 + (end - start + 1) / shift);
		else
		{
			sliceCount = frame;
			shift = 1; 
		}
		
		//System.out.println("number of slices: " + sliceCount);
		
		double[] sliceHistogram = new double[sliceCount];
		double[] histogram = new double[sliceCount];
		
		String CSVoutput = filename.substring(0, filename.indexOf('.')) + CSVName + ".csv"; 
		
		PrintWriter csv = new PrintWriter(new FileWriter(CSVoutput)); 
	    csv.println("Frame #,GC Percentage"); 
		
		for(int i = 0; i < sliceCount; i++)
		{
			//fill the sliceHistogram
			//if it's the last index, go to the last nucleotide instead of out of bounds
			//System.out.println("Subslice from codon " + i * shift + " to codon " + Math.min((i + 1) * shift, length));

			
			sliceHistogram[i] = gcCount(sub.substring(i * shift, Math.min((i + 1) * shift, length)));
			//System.out.println(sliceHistogram[i] + "%"); 
		}
		

		
		for(int i = 0; i < sliceCount; i++)
		{
			double currentCount = 0; 
			
			if(i * shift > end - frame)
				break; 
			
			//Algorithmically, there's a better way to do this
			//Make the frame by averaging all of the shifts it's made of 
			for(int j = 0; j < (frame / shift); j++)
			{
				if(j + i == sliceHistogram.length)
					break; 
				
				currentCount += sliceHistogram[j + i];
			}
			//System.out.println( i + " " + sliceHistogram[0]);
			
			histogram[i] = (double)(currentCount / (frame/shift));
		
			//System.out.println("Frame #" + i + ": " + histogram[i] + "%");   
			// FILE OUTPUT//
		    csv.println(i + "," + histogram[i] + "%"); 
		}   
		csv.close();
		
		return histogram; 
	}
	
	public static double gcSubCount(String filename, int start, int end) {
		long total= 0;

		String sub = readGene(filename).substring(start, end);
		
		return gcCount(sub);
	}
	
	public static double gcCount(String dna) {
		long total= 0;

		for (char c : dna.toCharArray()) {
			if (c=='C' || c=='G')
				total++;
		}
		return (double) total/dna.length()*100;
	}
}
