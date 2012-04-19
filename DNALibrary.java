import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.lang.StringBuilder;

import java.lang.Math; 
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.Scanner;


public static class DNALibrary {
   public static File inputFile;
	public static Map<String, String[]> aa2codon = new HashMap<String, String[]>();
	public static int[][][] freq = new int[4][4][4];
	
	public static void main(String[] args) throws NumberFormatException, IOException {
			int i = 0, frame = 0, shift = 0; 
			boolean valid = false; 
			Scanner sc = new Scanner(System.in); 
			String filename = null;
			
			System.out.println("GC Content Analyzer, now FASTA only!");
			
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
				else System.out.println("Shift and frame can't be 0");
			}
			valid = false; 
			
			//Read Filename

			System.out.println("What GFF file would you like to analyze?");
			filename = sc.next(); 			

			PrintWriter csv = new PrintWriter(new FileWriter(filename.substring(0, filename.lastIndexOf('.')) + ".csv"));

			csv.println("Start Codon,End Codon,GC Percentage"); 			

			Strand[] files = readStrandsFromFAFSA(filename);
			
			i = 0;
			for(Strand gene : files)
			{
				Strand[] frames;
				i++;
				System.out.println(gene.name + " " + gene.start + " " + gene.end);
				
				String input = readContig(gene.name + ".txt", gene.start, gene.end);
				gene.bases = input;
				frames = gcPercentage(gene, frame, shift);
				
				for(Strand section : frames)
				{
					csv.printf("%d,%d,%.1f%%\n", gene.start + section.start, gene.start + section.end, section.gcPercent);
					System.out.printf("%d,%d,%.1f%%\n", gene.start + section.start, gene.start + section.end, section.gcPercent);
				}
				csv.println();
			}
			csv.close();
			
			System.out.println("Done! File " +filename.substring(0, filename.lastIndexOf('.')) + ".csv" + " written.");
			System.out.println("Done. Goodbye.");
	}
	
	/**
	 * Returns the reverse complement of a strand of DNA
	 */
	public static String reverseComplement(String gene)
	{
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
	 * @return
	 */
	public static String readContig(String filename, int start, int end)
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
				return fileContents.toString(); 
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
		
		return fileContents.toString(); 
	}

	public static Strand[] readStrandsFromGFF(String filename)
	{	
		ArrayList<Strand> strandsAL = new ArrayList<Strand>();  
		
		Strand[] strands = new Strand[0]; 
		
		try
		{
			File gene = new File(filename);
			Scanner sc = new Scanner(gene);

			sc.nextLine(); //ignore first line: browser position
			sc.nextLine(); // track name, description, visibility, color
			
			while (sc.hasNextLine()) 
			{
				Strand current = new Strand(); 
				
				current.name = sc.next();
				sc.next(); // GEP
				sc.next(); //mRNA or CDS
				current.start = sc.nextInt(); 
				current.end = sc.nextInt(); 
				
				//System.out.println(current.start + " " + current.end);
				
				sc.next(); //dot
				current.direction = sc.next().charAt(0); 
				sc.nextLine(); //dot
				current.length = current.end - current.start;
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
			System.out.println(ex + " - Invalid Strand filename?");
		}
		
		strands = strandsAL.toArray(strands);
		
		return strands;
	}
	
	public static Strand[] readStrandsFromFAFSA(String filename)
	{	
		ArrayList<Strand> strandsAL = new ArrayList<Strand>();  
		
		Strand[] strands = new Strand[0]; 
		
		try
		{
			File gene = new File(filename);
			Scanner sc = new Scanner(gene);
			
			while (sc.hasNextLine()) 
			{
				Strand current = new Strand(); 
				
				try{
				current.name = sc.next();
				sc.next(); //dot
				current.type = sc.next(); //mRNA or CDS
				current.start = sc.nextInt(); 
				current.end = sc.nextInt(); 
				sc.next(); //dot
				current.direction = sc.next().charAt(0); 
				sc.next(); //dot
				sc.next(); //"gene_id" hopefully doesn't change
				current.id = sc.next();
				sc.next(); //transcript_id
				current.transcriptId = sc.nextLine();
				
				current.id = current.id.substring(1, current.id.length() - 2);
				current.transcriptId = current.transcriptId.substring(1, current.transcriptId.length() - 2);
				
				System.out.println(current.name + " " + current.type + " " + current.start + " " + current.end + " " + current.direction + " " + current.id + " " + current.transcriptId);
				
				
				current.length = current.end - current.start;
				if(current.end < current.start)
				{
					int temp = current.start;
					current.start = current.end;
					current. end = temp; 
				}
				strandsAL.add(current); 
				}
				catch(Exception ex)
				{
					System.err.println(ex);
					System.err.println("Problem reading FAFSA file");
				}
				
			}
			sc.close(); 
		}
		catch(Exception ex)
		{
			System.out.println(ex + " - Invalid Strand filename?");
		}
		
		strands = strandsAL.toArray(strands);
		
		return strands;
	}
	
	public static Strand[] gcPercentage(Strand gene, int frame, int shift)
	{
		ArrayList<Strand> strands = new ArrayList<Strand>(); 
		
		double[] slices = new double[(int) (1 + ((gene.length/ shift)))]; 
		
		if(shift > frame)
		{
			System.err.println("Shift " + shift + " is larger than frame " + frame);
			return null; 
		}
		
		//Compute percentages of each slice
		int i, j;
		//Short slices are only one frame, don't need recalculating
		if(frame >= gene.length)
		{
			Strand[] single = new Strand[1];
			Strand shorty = new Strand();
			shorty.gcPercent = gcCount(gene.bases);
			shorty.start = 0;
			shorty.end = gene.length;
			single[0] = shorty;
			return single;
		}
		//If the shift divides the frame evenly, we can speed things up
		else if(frame % shift == 0)
		{
			for(i = 0; i < slices.length && (i - 1) * shift < gene.length; i++)
			{
				slices[i] = gcTotal(gene.bases.substring(i * shift, Math.min(gene.length, (i + 1) * shift)));
			}
			for(i = 0; (i - 1) * shift + frame < gene.length; i++)
			{
				Strand current = new Strand();
				current.start = i * shift;
				current.end = Math.min(gene.length, i * shift + frame); 
				current.length = current.end - current.start;
				int numerator = 0;
				for(j = i; j < frame/shift + i; j++)
				{
					numerator += slices[j];
				}
				
				current.gcPercent = (numerator  * 100.0 / current.length);
			//	System.out.println(current.gcPercent);
				strands.add(i, current);
			}
		}
		else
		{
			for(i = 0; i * shift + frame < gene.length; i++)
			{
				Strand current = new Strand();
				current.start = i * shift;
				current.end = i * shift + frame; 

				current.gcPercent = gcCount(gene.bases.substring(current.start, current.end));
				strands.add(i, current);
				
				//last stubby frame
				if((i + 1) * shift + frame > gene.length)
				{
					current.start = i * shift;
					current.end = gene.length;

					current.gcPercent = gcCount(gene.bases.substring(current.start, current.end));
					strands.add(i, current);
					
				}
			}
		}
		Strand[] ret = new Strand[strands.size()];
		return strands.toArray(ret);
	}

	public static double gcCount(String dna) {
		long total= 0;

		for (char c : dna.toCharArray()) {
			if (c=='C' || c=='G')
				total++;
		}
		return (double) total/dna.length()*100;
	}
	
	public static int gcTotal(String dna) {
		int total= 0;

		for (char c : dna.toCharArray()) {
			if (c=='C' || c=='G')
				total++;
		}
		return total;
	}
	
	public static double chiSq(String aminoAcid) {
		String[] codons = aa2codon.get(aminoAcid);
		int degenerate = codons.length;
		int[] freqs = new int[degenerate];
		int total = 0;
		double e;
		double sum = 0.0;

		for (int i = 0; i < degenerate; i++) {
			freqs[i] = getFreq(codons[i]);
			total += freqs[i];
		}
		e = (double) total / degenerate;
		//System.out.print("Expecting " + e + " for " + aminoAcid);

		// Summation of ((# of occurences - E)^2)/E
		for (int i = 0; i < degenerate; i++) {
			sum += Math.pow(freqs[i]-e, 2.0) / e;
		}

		return sum;
	}

	public static int getFreq(String codon) {
		return freq[ntoi(codon.charAt(0))][ntoi(codon.charAt(2))][ntoi(codon.charAt(1))];
	}
	
	public static int ntoi(char neucleotide) {
		if (neucleotide == 'A')
			return 0;
		else if (neucleotide == 'C')
			return 1;
		else if (neucleotide == 'G')
			return 2;
		return 3;
	}
	
	public static void initAA() {
		aa2codon.put("I", new String[]{"ATT", "ATC", "ATA"});
		aa2codon.put("L", new String[]{"CTT", "CTC", "CTA", "CTG", "TTA", "TTG"});
		aa2codon.put("V", new String[]{"GTT", "GTC", "GTA", "GTG"});
		aa2codon.put("F", new String[]{"TTT", "TTC"});
		aa2codon.put("M", new String[]{"ATG"});
		aa2codon.put("C", new String[]{"TGT", "TGC"});
		aa2codon.put("A", new String[]{"GCT", "GCC", "GCA", "GCG"});
		aa2codon.put("G", new String[]{"GGT", "GGC", "GGA", "GGG"});
		aa2codon.put("P", new String[]{"CCT", "CCC", "CCA", "CCG"});
		aa2codon.put("T", new String[]{"ACT", "ACC", "ACA", "ACG"});
		aa2codon.put("S", new String[]{"TCT", "TCC", "TCA", "TCG", "AGT", "AGC"});
		aa2codon.put("Y", new String[]{"TAT", "TAC"});
		aa2codon.put("W", new String[]{"TGG"});
		aa2codon.put("Q", new String[]{"CAA", "CAG"});
		aa2codon.put("N", new String[]{"AAT", "AAC"});
		aa2codon.put("H", new String[]{"CAT", "CAC"});
		aa2codon.put("E", new String[]{"GAA", "GAG"});
		aa2codon.put("D", new String[]{"GAT", "GAC"});
		aa2codon.put("K", new String[]{"AAA", "AAG"});
		aa2codon.put("R", new String[]{"CGT", "CGC", "CGA", "CGG", "AGA", "AGG"});
		aa2codon.put("STOP", new String[]{"TAA", "TAG", "TGA"});
	}

	
	public static class Strand
	{
		String bases; 
		//Name of the file this strand came from
		String name;
		//number of nucleotides long
		int length;
		//start position on gene
		int start;
		//end position on gene
		int end;
		//mRNA or CDS
		String type;
		//+ or -
		char direction;
		//GC percentage if calculated
		double gcPercent;
		//ID of this gene
		String id;
		//Transcript ID
		String transcriptId;
		
		public Strand()
		{
		}
		
		public Strand(int start, int end, char direction, String bases)
		{
			this.start = start;
			this.end = end;
			this.direction = direction;
			this.length = end - start;
			
			if(this.direction != '+' && this.direction != '-')
				System.err.println("Strand made with odd direction " + direction);
		}
		
		public Strand (String filename, int start, int end, boolean reverseComplement)
		{
			int i; 
			StringBuilder fileContents = new StringBuilder();

			this.start = start;
			this.end = end;
			this.length = end - start;
			this.name = filename.substring(0, filename.indexOf("."));
			
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
					
					return;
				}
				//Read the rest of the lines until the line before the end
				else 
				{
					for(; i < end - 50; i += 50)		
						fileContents.append(sc.nextLine());
			
					if(sc.hasNextLine())
						fileContents.append(sc.nextLine().substring(0, end - i));
				}
				sc.close(); 
			}
			catch(FileNotFoundException ex)
			{
				System.out.println(ex + " - Invalid filename?");
				System.out.println("There's something wrong with the way you typed " + filename + " or something bad happened");
			}
			this.bases = fileContents.toString(); 
		}
		
	}
}


