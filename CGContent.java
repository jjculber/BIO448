import java.io.File;
import java.io.FileNotFoundException;
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
			int i = 0, frame = 0, shift = 0; 
			boolean valid = false; 
			char isGFF = '\0'; 
			Scanner sc = new Scanner(System.in); 
			String filename = null;
			
			System.out.println("GC Content Analyzer");
			
			System.out.println("Are you reading a GFF File? (y/n)");
			
			//GFF File check - UI input
			while(!valid)
			{
			    isGFF = (char) System.in.read();
				
				if(Character.toLowerCase(isGFF) == 'y' || Character.toLowerCase(isGFF) == 'n')
					valid = true;
				
				else System.out.print("Let's try that again. GFF file - Yes (y) or no (n)? ");
			}
			valid = false; 
			
			//Frame and range check - UI input
			while(!valid)
			{
				System.out.println("How many base pairs should the frame be?");
				frame = sc.nextInt();
				
				System.out.println("How many base pairs should it be shifted by?");
				shift = sc.nextInt();
				
				if(shift > 0 && frame > 0)
					valid = true; 
			}
			valid = false; 
			
			//Read Filename

			System.out.println("What file would you like to analyze?");
			filename = sc.next(); 			

			PrintWriter csv = new PrintWriter(new FileWriter(filename.substring(0, filename.lastIndexOf('.')) + ".csv"));

			csv.println("Start Codon,End Codon,GC Percentage"); 			

			if(Character.toLowerCase(isGFF) == 'y')
			{
				GFF[] files = readGFF(filename);
				
				i = 0;
				for(GFF gene : files)
				{
					i++;
					System.out.println(gene.filename + " " + gene.start + " " + gene.end);
					
					String input = readGene(gene.filename, gene.start, gene.end, !gene.direction);
					
					
					gcPercentCSV(input, Math.min(frame, gene.end - gene.start), shift, csv);
				}
			}
			else 
			{
				int start = 0, end = 0;
				
				while(!valid)
				{
					
					System.out.println("What start position?");
					start = sc.nextInt();
					
					System.out.println("What end position?");
					end = sc.nextInt();
					
					if(shift > 0 && frame > 0)
						valid = true; 
				}

				String input = readGene(filename, start, end, false);
					
				gcPercentCSV(input, frame, shift, csv);
			}
			System.out.println("Done! File " + csv.toString() + " written.");
			System.out.println("Done. Goodbye.");
	}
	
	/**
	 * Returns the reverse complement of a strand of DNA
	 */
	public static String reverseComplement(String gene)
	{
				
		int length = gene.length();
		
		StringBuilder reverse = new StringBuilder();
		
		for(int i = gene.length(); i > 0; --i)
		{
			char result = gene.charAt(i-1);
			
			reverse.append(complement(result));
		}

		return reverse.toString(); 
	}
	
	/**
	 * Given a nucleotide as a char, returns its complement
	 * @param nucleotide 'A' 'C' 'G' or 'T' 
	 * @return 'A' 'C' 'G' or 'T'
	 */
	public static char complement(char nucleotide)
	{
		switch (Character.toUpperCase(nucleotide)) {
		case 'A':
			return 'T';
		case 'T':
			return 'A';
		case 'C':
			return 'G';
		case 'G':
			return 'C';
		default:
			if(nucleotide != ' ')
				System.err.println("Non-nucleotide " + nucleotide + " found");
			return nucleotide;
		}
	}
	
	/**
	 * Reads a contig (?) file and outputs the nucleotide sequence in the specified range
	 * @param filename
	 * @param start
	 * @param end
	 * @param reverseComplement
	 * @return
	 */
	public static String readGene(String filename, int start, int end, boolean reverseComplement)
	{
		int i; 
		StringBuilder fileContents = new StringBuilder();
		
		File gene;
		
		try
		{
			gene = new File(filename);
			Scanner sc = new Scanner(gene);
			
			sc.nextLine(); //ignore first line
			
			for(i = 0; i < start - 50; i += 50 )
			{
				sc.nextLine(); 
			}
			
			//Add first line, from the first codon wanted to the end of the line
			if(end - start > 50)
			{
				fileContents.append(sc.nextLine().substring(start - i));
			}
			else if(end - i < 50)
			{
				fileContents.append(sc.nextLine().substring(start - i, end  - i));
				return (reverseComplement) ? reverseComplement(fileContents.toString()) : fileContents.toString(); 
			}
			//Read the rest of the lines until the line before the end
			for(; i < end - 50; i += 50)
				fileContents.append(sc.nextLine());
			if(sc.hasNextLine())
				fileContents.append(sc.nextLine().substring(0, end - i));
			
			sc.close(); 
		}
		catch(FileNotFoundException ex)
		{
			System.out.println(ex + " - Invalid filename?");
			System.out.println("There's something wrong with the way you typed " + filename + " or something bad happened");
		}
		
		return (reverseComplement) ? reverseComplement(fileContents.toString()) : fileContents.toString(); 
	}
	

	public static GFF[] readGFF(String filename)
	{	
		ArrayList<GFF> strandsAL = new ArrayList<GFF>();  
		
		GFF[] strands = new GFF[0]; 
		
		try
		{
			File gene = new File(filename);
			Scanner sc = new Scanner(gene);

			sc.nextLine(); //ignore first line: browser position
			sc.nextLine(); // track name, description, visibility, color
			
			while (sc.hasNextLine()) 
			{
				GFF current = new GFF(); 
				
				current.filename = sc.next() + ".txt";
				sc.next(); // GEP
				sc.next(); //mRNA or CDS
				current.start = sc.nextInt(); 
				current.end = sc.nextInt(); 
				sc.next(); //dot
				current.direction = (sc.next().equals("+")); 
				sc.nextLine(); //dot
				
				if(current.end < current.start)
				{
					int temp = current.start;
					current.start = current.end;
					current. end = temp; 
				}
				strandsAL.add(current); 
				
			}
			sc.close(); 
		}
		catch(Exception ex)
		{
			System.out.println(ex + " - Invalid GFF filename?");
		}
		
		strands = strandsAL.toArray(strands);
		
		return strands;
	}
	
	/**
	 * Computes the GC percentages of a gene and outputs a CSV file of it
	 * @param gene
	 * @param frame
	 * @param shift
	 * @param outputFile
	 */
	public static void gcPercentCSV(String gene, int frame, int shift, PrintWriter csv)
	{
		
		//String CSVoutput = filename.substring(0, filename.indexOf('.')) + CSVName + ".csv"; 
		try 
		{
			
			double[] percentages = gcFrameCount(gene, frame, shift);
			
			for(int i = 0; i < percentages.length; i++)
			{
				 csv.println(i + 1 + "," + percentages[i] + "%"); 
			}
		
			
			csv.close(); 
		} 
		catch (IOException e) {
			System.out.println("CSV output error");
			e.printStackTrace();
		} 
	}
	
	/**
	 * gcFramecount
	 * Outputs a CSV histogram of gc percentages over a genome
	 * @throws IOException 
	 */
	public static double[][][] gcFrameCount(String gene, int frame, int shift) throws IOException
	{
		long total= 0;
		int sliceCount;
		int length = gene.length(); 
		
		
		//fine histogram of shift-width percentages
	
		
		if(shift != 0)
			sliceCount = (int)(.5 + (length + 1) / shift);
		else
		{
			sliceCount = frame;
			shift = 1; 
		}
		
		//System.out.println("number of slices: " + sliceCount);
		
		double[] sliceHistogram = new double[sliceCount];
		double[] histogram = new double[sliceCount];

		
		for(int i = 0; i < sliceCount; i++)
		{
			//fill the sliceHistogram
			//if it's the last index, go to the last nucleotide instead of out of bounds
			//System.out.println("Subslice from codon " + i * shift + " to codon " + Math.min((i + 1) * shift, length));

			
			sliceHistogram[i] = gcCount(gene.substring(i * shift, Math.min((i + 1) * shift, length)));
			//System.out.println(sliceHistogram[i] + "%"); 
		}

		for(int i = 0; i < sliceCount; i++)
		{
			double currentCount = 0; 
			
			if(i * shift > length - frame)
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
		}   
		
		return histogram; 
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
