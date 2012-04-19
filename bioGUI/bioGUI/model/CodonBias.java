package bioGUI.model;

import java.io.BufferedWriter;
import java.util.HashMap;
import java.util.Map;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.Scanner;
import java.lang.Math;

public class CodonBias {
	public static Map<String, String[]> aa2codon = new HashMap<String, String[]>();
	public static int[][][] freq = new int[4][4][4];


	public static void main(String[] args) throws IOException {
		File gffFile;
		FileWriter outFile;

		StringBuffer output = new StringBuffer();
		int startPos, stopPos;
		Scanner sc = new Scanner(System.in);
		final String START = "ATG";
		final List<String> STOP = Arrays.asList("TAA", "TAG", "TGA");
		
		// USER INPUT
		System.out.println("Please enter the GFF file name:");
		gffFile = new File(sc.nextLine());

		GFF[] exons = CGContent.readGFF(gffFile.getAbsolutePath());

		String path = gffFile.getAbsolutePath();
		String fileSuffix = "_"+path.substring(path.lastIndexOf("\\")+1, path.length()-4) + ".csv";
		path = path.substring(0, path.lastIndexOf("\\")+1);
		
		System.out.println(exons[0].filename);
		String fastaFileName = path + exons[0].filename + ".txt";
		System.out.println("Looking for FASTA file at " + fastaFileName);
		File fastaFile = new File(fastaFileName);

		System.out.println("Saving files in " + path);

		String outFileName = path + "genes" + fileSuffix;
		outFile = new FileWriter(outFileName);
		BufferedWriter out = new BufferedWriter(outFile);
		System.out.println("Please enter the start position(use 0 for beginning):");
		startPos = sc.nextInt();
		System.out.println("Please enter the stop position(use -1 for end):");
		stopPos = sc.nextInt();
		
		
		// SET UP INPUT FILE
		StringBuffer input = new StringBuffer();
		sc = new Scanner(fastaFile);
		
		sc.nextLine();
		while (sc.hasNextLine()) {
			input.append(sc.nextLine());
		}
		
		String dna = input.toString();
		if(stopPos == -1)
			stopPos = dna.length();

		//TODO loop through GFFs
		
		// FIND GENES AND COUNT NEUCLEOTIDES
		dna = dna.substring(startPos, stopPos);
		StringBuffer concat = new StringBuffer();
		for (GFF gff : exons) {
			System.out.println(gff.start + " " + gff.end + " " + gff.direction);
			if (gff.start < dna.length() && gff.end < dna.length() && gff.start < gff.end) {
				concat.append(dna.substring(gff.start, gff.end));
			}
		}

		String all = concat.toString();
		//System.out.println(exons.length);
		//System.out.println(all.length());

		int i = 0;
		while (i < all.length()-2) {
			freq[ntoi(all.charAt(i))][ntoi(all.charAt(i+1))][ntoi(all.charAt(i+2))]++;
			i+=3;
		}
		
		out.write(output.toString().toCharArray());
		out.flush();
		out.close();
		
		outFileName = path + "CodonFrequency"+fileSuffix;
		outFile = new FileWriter(outFileName);
		out = new BufferedWriter(outFile);
		output = new StringBuffer();

		// PRINT ALL CODONS AND COUNTS
		for (char one : Arrays.asList('T', 'C', 'A', 'G')) {
			for (char two : Arrays.asList('T', 'C', 'A', 'G')) {
				for (char three : Arrays.asList('T', 'C', 'A', 'G')) {
					output.append(""+one+three+two+","+getFreq(""+one+two+three)+"\n");
				}
				//System.out.println();
			}	
			//System.out.println();
		}
		
		out.write(output.toString().toCharArray());
		out.flush();
		out.close();
		
		outFileName = path + "CodonBias" + fileSuffix;
		outFile = new FileWriter(outFileName);
		out = new BufferedWriter(outFile);
		output = new StringBuffer();
	
		// PRINT CODON AND COUNT GROUPED BY AMINO ACID
		initAA();	
		for (String aa : aa2codon.keySet()) {
			output.append(aa + "," + chiSq(aa)+"\n");
			String[] codons = aa2codon.get(aa);
			int total = 0;
			for (String codon : codons) 
				total += getFreq(codon);
			
			for (String codon : codons){
				output.append(",," + codon +","+getFreq(codon) + "," + ((double) getFreq(codon)/total*100.0)+"\n");
			}
		}
		
		out.write(output.toString().toCharArray());
		out.flush();
		out.close();
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

	public static String reverse(String dna) {
		return null;
	}

	public static String complement(String dna) {
		return null;
	}

	public static String reverseComplement(String dna) {
		return reverse(complement(dna));
	}

	public static void findGenes(String dna) {
		/*
		 * old code. Will probably never use.
		
		// FIND GENES AND COUNT NEUCLEOTIDES
		dna = dna.substring(startPos, stopPos);
		int i = 0;
		boolean inGene = false;
		//System.out.println(dna);
		output.append("start,stop\n");
		while (i < dna.length()-2) {
			if (inGene) {
				//if we are looking at the contents of the gene
				if (STOP.contains(dna.substring(i, i+3))) {
					output.append(i+"\n");
					freq[ntoi(dna.charAt(i))][ntoi(dna.charAt(i+1))][ntoi(dna.charAt(i+2))]++;
					inGene = false;
					i+=3;
				} else {
					//this is a codon
					freq[ntoi(dna.charAt(i))][ntoi(dna.charAt(i+1))][ntoi(dna.charAt(i+2))]++;
					i+=3;
				}
			} else {
				//if we haven't found the start codon yet
				if (dna.substring(i, i + 3).equals(START)) {
					output.append(i + ",");
					inGene = true;
					i+=3;
				} else {
					i++;
				}
			}
		}
		*/
	}
}
