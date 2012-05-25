/**
 * DNALibrary.java
 * 
 * @author jtbiddle
 * @author jjculber
 * @author tshopshire
 * 
 * DNA analysis library for CPE448, Bioinformatics Algorithms
 * 
 * CSC Members: John Biddle, Justin Culbertson, and Tyler Shopshire
 * CHEM Members: Nora Goscinski and Matthew Randazzo
 */
package bioGUI.model;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.lang.StringBuilder;

import java.lang.Math;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Scanner;

import javax.swing.JOptionPane;
import javax.swing.event.ListSelectionEvent;

public class DNALibrary {
	public static File inputFile;
	//yes, this is a global variable.
	public static ArrayList<String> errors;
	
	public static void main(String[] args) throws IOException {
		System.out.println("MASTER CONCATENATOR");

		System.out.println("WHAT FOLDER WOULD YOU LIKE, BROTHER?");

		Scanner sc = new Scanner(System.in);

		String file = sc.next();

		concatChromosome(file);

		
	}

	/**
	 * Given a nucleotide as a char, returns its complement
	 * 
	 * @param nucleotide
	 *            'A' 'C' 'G' or 'T'
	 * @return 'A' 'C' 'G' or 'T'
	 */
	public static char complement(char nucleotide) {
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
			if (nucleotide != ' ')
				popupError("Non-nucleotide " + nucleotide + " found");
			return nucleotide;
		}
	}

	/**
	 * Returns the reverse complement of a strand of DNA
	 */
	public static String reverseComplement(String gene) {
		StringBuilder reverse = new StringBuilder();

		for (int i = gene.length(); i > 0; --i) {
			char result = gene.charAt(i - 1);
			reverse.append(complement(result));
		}
		return reverse.toString();
	}

	/**
	 * optimalFrequency Calculates the frequency of optimal codons
	 */
	public static double optimalFrequency(String dna, List<String> optimalCodons) {
		List<String> Xex = getExcludedCodons();
		List<String> codons = getAllCodons();

		List<String> Xnon = codons;
		Xnon.removeAll(optimalCodons);

		int optimal = 0;
		int nonoptimal = 0;

		for (int i = 3; i < dna.length(); i += 3) {
			String current = dna.substring(i - 3, i);
			if (optimalCodons.contains(current))
				optimal++;
			else if (Xnon.contains(current))
				nonoptimal++;
		}

		return (optimal / (nonoptimal + optimal));
	}

	/**
	 * 
	 * Relative Synonymous Codon Usage
	 * 
	 * @param dna
	 *            DNA sequence
	 * @param codon
	 *            The codon to compute the RCSU of
	 * @param acidCodons
	 *            All of the codons that code for this amino acid
	 * @return
	 */
	public static double RCSU(String dna, String codon, List<String> acidCodons) {
		int Xi = 0, sum = 0;
		int[] Xj = new int[acidCodons.size()];
		Arrays.fill(Xj, 0);

		for (int i = 3; i < dna.length(); i += 3) {
			String current = dna.substring(i - 3, i);
			if (current.equals(codon))
				Xi++;
			else if (acidCodons.contains(current))
				Xj[acidCodons.indexOf(current)]++;
		}

		for (int i = 0; i < Xj.length; i++) {
			sum += Xj[i];
		}

		return (Xi / ((1 / acidCodons.size()) * sum));
	}

	/**
	 * Codon Adaptation Index
	 * 
	 * @return
	 */
	public static double codonAdaptationIndex(String dna) {
		List<String> codons = getAllCodons();
		codons.removeAll(getExcludedCodons());

		Collection<String[]> aminoCodons = initAA().values();

		double[] rcsus = new double[codons.size()];

		for (int i = 0; i < codons.size(); i++) {
			for (String[] s : aminoCodons) {
				if (Arrays.asList(s).contains(codons.get(i))) {
					rcsus[i] = RCSU(dna, codons.get(i), Arrays.asList(s));
				}
			}
		}

		double rcsuMax = rcsus[arrayMaxIndex(rcsus)];

		double sum = 0;

		for (int i = 0; i < rcsus.length; i++) {
			sum += Math.log(rcsus[i] / rcsuMax);
		}

		return Math.exp((1 / rcsus.length) * sum);

	}

	public static List<String> getAllCodons() {
		String[] acids = { "TTT", "TTC", "TTA", "TTG", "TCT", "TCC", "TCA",
				"TCG", "TAT", "TAC", "TGT", "TGC", "TGG", "CTT", "CTC", "CTA",
				"CTG", "CCT", "CCC", "CCA", "CCG", "CAT", "CAC", "CAA", "CAG",
				"CGT", "CGC", "CGA", "CGG", "ATT", "ATC", "ATA", "ATG", "ACT",
				"ACC", "ACA", "ACG", "AAT", "AAC", "AAA", "AAG", "AGT", "AGC",
				"AGA", "AGG", "GTT", "GTC", "GTA", "GTG", "GCT", "GCC", "GCA",
				"GCG", "GAT", "GAC", "GAA", "GAG", "GGT", "GGC", "GGA", "GGG" };

		return Arrays.asList(acids);
	}

	public static List<String> getExcludedCodons() {
		List<String> Xex = new ArrayList<String>();
		// Stop codons
		Xex.add("TAG");
		Xex.add("TGA");
		Xex.add("TAA");
		// Methionine
		Xex.add("ATG");
		// Tryptophan
		Xex.add("TGG");

		return Xex;
	}

	/**
	 * calcGCContent Provides an interface to calculate the GC content of a
	 * FASTA file to our group members' specifications
	 * 
	 * Outputs a csv file of the gc content of the fasta file
	 * 
	 * @param fastaFilename
	 *            FASTA input filename
	 * @param start
	 *            Start base pair
	 * @param end
	 *            End base pair
	 * @param frameSize
	 *            The size of the sampling window
	 * @param frameShift
	 *            How many bp to move the sampling window by
	 */
	public static void calcGCContent(String fastaFilename, int start, int end,
			int frameSize, int frameShift) {
		try {
			String dna = readFASTA(fastaFilename, start, end);

			Strand gene = new Strand(start, end, '+', dna);
			Strand[] frames = gcPercentage(gene, frameSize, frameShift);

			File fastaFile = new File(fastaFilename);
			String path = fastaFile.getParent() + File.separator;

			String outFileName = path + fastaFile.getName() + "_GCContent.csv";
			FileWriter outFile = new FileWriter(outFileName);
			BufferedWriter out = new BufferedWriter(outFile);
			StringBuffer output = new StringBuffer();
			output.append("Start Pos, End Pos, CG Percent\n");

			for (Strand frame : frames) {
				output.append(frame.start + "," + frame.end + ","
						+ frame.gcPercent + "\n");
			}

			out.write(output.toString().toCharArray());
			out.flush();
			out.close();

			popupMessage("Ouput CSV file saved in " + path);
		} catch (Exception e) {
			e.printStackTrace();
			popupError("An Error occurred!\n" + e.getMessage());
		}
	}

	/**
	 * Calculates the GC Percent. Returns a double between 0.0 and 100.0
	 */
	public static double gcCount(String dna) {
		long total = 0;

		for (char c : dna.toCharArray()) {
			if (c == 'C' || c == 'G')
				total++;
		}
		return (double) total / dna.length() * 100;
	}

	/**
	 * Returns the gc percentage histogram of a DNA sequence as an array of
	 * Strands
	 * 
	 * @param gene
	 *            The DNA sequence to analyze
	 * @param frame
	 *            The sampling window size
	 * @param shift
	 *            Samping window shift distance
	 * @return An array of strands, each with a length of the frame size.
	 */
	public static Strand[] gcPercentage(Strand gene, int frame, int shift) {
		ArrayList<Strand> strands = new ArrayList<Strand>();

		double[] slices = new double[(int) (1 + ((gene.length / shift)))];

		if (shift > frame) {
			popupError("Shift " + shift + " is larger than frame " + frame);
			return null;
		}

		// Compute percentages of each slice
		int i, j;
		// Short slices are only one frame, don't need recalculating
		if (frame >= gene.length) {
			Strand[] single = new Strand[1];
			Strand shorty = new Strand();
			shorty.gcPercent = gcCount(gene.bases);
			shorty.start = 0;
			shorty.end = gene.length;
			single[0] = shorty;
			return single;
		}
		// If the shift divides the frame evenly, we can speed things up
		else if (frame % shift == 0) {
			for (i = 0; i < slices.length && (i - 1) * shift < gene.length; i++) {
				slices[i] = gcTotal(gene.bases.substring(i * shift,
						Math.min(gene.length, (i + 1) * shift)));
			}
			for (i = 0; (i - 1) * shift + frame < gene.length; i++) {
				Strand current = new Strand();
				current.start = i * shift;
				current.end = Math.min(gene.length, i * shift + frame);
				current.length = current.end - current.start;
				int numerator = 0;
				for (j = i; j < frame / shift + i; j++) {
					numerator += slices[j];
				}

				current.gcPercent = (numerator * 100.0 / current.length);
				strands.add(i, current);
			}
		} else {
			for (i = 0; i * shift + frame < gene.length; i++) {
				Strand current = new Strand();
				current.start = i * shift;
				current.end = i * shift + frame;

				current.gcPercent = gcCount(gene.bases.substring(current.start,
						current.end));
				strands.add(i, current);

				// last stubby frame
				if ((i + 1) * shift + frame > gene.length) {
					current.start = i * shift;
					current.end = gene.length;

					current.gcPercent = gcCount(gene.bases.substring(
							current.start, current.end));
					strands.add(i, current);

				}
			}
		}
		Strand[] ret = new Strand[strands.size()];
		return strands.toArray(ret);
	}

	/*
	 * Calculates the total number of Gs or Cs in a DNA String
	 */
	public static int gcTotal(String dna) {
		int total = 0;

		for (char c : dna.toCharArray()) {
			if (c == 'C' || c == 'G')
				total++;
		}
		return total;
	}

	/*
	 * Given a file name, calculate the codon bias and output to file
	 */
	public static void calcCodonBias(String filename) {
		try {
			File gffFile;
			FileWriter outFile;
			String outFileName;
			BufferedWriter out;
			StringBuffer output = new StringBuffer();
			Scanner sc;

			gffFile = new File(filename);

			Strand[] exons = DNALibrary.readStrandsFromGFF(filename);

			String path = gffFile.getParent() + File.separator;

			String fastaFileName = path + exons[0].name + ".txt";
			File fastaFile = new File(fastaFileName);

			// SET UP INPUT FILE
			StringBuffer input = new StringBuffer();
			sc = new Scanner(fastaFile);

			sc.nextLine();
			while (sc.hasNextLine()) {
				input.append(sc.nextLine());
			}

			String dna = input.toString();

			// FIND GENES AND COUNT NEUCLEOTIDES
			StringBuffer concat = new StringBuffer();
			for (Strand gff : exons) {
				System.out.println(gff.start + " " + gff.end + " "
						+ gff.direction);
				if (gff.start < dna.length() && gff.end < dna.length()) {
					if (gff.start < gff.end) {
						concat.append(dna.substring(gff.start, gff.end));
					} else {
						concat.append(DNALibrary.reverseComplement(dna
								.substring(gff.end, gff.start)));
					}
				}
			}

			String all = concat.toString();

			int[][][] freq = new int[4][4][4];

			int i = 0;
			while (i < all.length() - 2) {
				freq[ntoi(all.charAt(i))][ntoi(all.charAt(i + 1))][ntoi(all
						.charAt(i + 2))]++;
				i += 3;
			}

			outFileName = path + gffFile.getName() + "_CodonFrequency.csv";
			outFile = new FileWriter(outFileName);
			out = new BufferedWriter(outFile);
			output = new StringBuffer();
			output.append("Codon, frequency\n");

			// PRINT ALL CODONS AND COUNTS
			for (char one : Arrays.asList('T', 'C', 'A', 'G')) {
				for (char two : Arrays.asList('T', 'C', 'A', 'G')) {
					for (char three : Arrays.asList('T', 'C', 'A', 'G')) {
						output.append("" + one + three + two + ","
								+ getFreq("" + one + two + three, freq) + "\n");
					}
					// System.out.println();
				}
				// System.out.println();
			}

			out.write(output.toString().toCharArray());
			out.flush();
			out.close();
			outFileName = path + gffFile.getName() + "_CodonBias.csv";

			outFile = new FileWriter(outFileName);
			out = new BufferedWriter(outFile);
			output = new StringBuffer();
			output.append("AA, chi squared, codon, frequency, percent\n");

			// PRINT CODON AND COUNT GROUPED BY AMINO ACID
			Map<String, String[]> aa2codon = initAA();
			for (String aa : aa2codon.keySet()) {
				output.append(aa + "," + chiSq(aa, freq) + "\n");
				String[] codons = aa2codon.get(aa);
				int total = 0;
				for (String codon : codons)
					total += getFreq(codon, freq);

				for (String codon : codons) {
					output.append(",," + codon + "," + getFreq(codon, freq)
							+ ","
							+ ((double) getFreq(codon, freq) / total * 100.0)
							+ "\n");
				}
			}

			out.write(output.toString().toCharArray());
			out.flush();
			out.close();
			popupMessage("Saved CSV output files in " + path);
		} catch (Exception e) {
			popupError("An error occurred!\n" + e.getMessage());
		}
	}

	/**
	 * chiSq Calculates the chi square of a sequence of amino acids
	 * 
	 * @param aminoAcid
	 * @param counts
	 * @return The chi sq of the amino acid sequence
	 */
	public static double chiSq(String aminoAcid, int[][][] counts) {
		Map<String, String[]> aa2codon = initAA();
		String[] codons = aa2codon.get(aminoAcid);
		int degenerate = codons.length;
		int[] freqs = new int[degenerate];
		int total = 0;
		double e;
		double sum = 0.0;

		for (int i = 0; i < degenerate; i++) {
			freqs[i] = getFreq(codons[i], counts);
			total += freqs[i];
		}
		e = (double) total / degenerate;

		// Summation of ((# of occurences - E)^2)/E
		for (int i = 0; i < degenerate; i++) {
			sum += Math.pow(freqs[i] - e, 2.0) / e;
		}

		return sum;
	}

	public static int getFreq(String codon, int[][][] freqTable) {
		return freqTable[ntoi(codon.charAt(0))][ntoi(codon.charAt(2))][ntoi(codon
				.charAt(1))];
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

	public static Map<String, String[]> initAA() {
		Map<String, String[]> aa2codon = new HashMap<String, String[]>();

		aa2codon.put("I", new String[] { "ATT", "ATC", "ATA" });
		aa2codon.put("L", new String[] { "CTT", "CTC", "CTA", "CTG", "TTA",
				"TTG" });
		aa2codon.put("V", new String[] { "GTT", "GTC", "GTA", "GTG" });
		aa2codon.put("F", new String[] { "TTT", "TTC" });
		aa2codon.put("M", new String[] { "ATG" });
		aa2codon.put("C", new String[] { "TGT", "TGC" });
		aa2codon.put("A", new String[] { "GCT", "GCC", "GCA", "GCG" });
		aa2codon.put("G", new String[] { "GGT", "GGC", "GGA", "GGG" });
		aa2codon.put("P", new String[] { "CCT", "CCC", "CCA", "CCG" });
		aa2codon.put("T", new String[] { "ACT", "ACC", "ACA", "ACG" });
		aa2codon.put("S", new String[] { "TCT", "TCC", "TCA", "TCG", "AGT",
				"AGC" });
		aa2codon.put("Y", new String[] { "TAT", "TAC" });
		aa2codon.put("W", new String[] { "TGG" });
		aa2codon.put("Q", new String[] { "CAA", "CAG" });
		aa2codon.put("N", new String[] { "AAT", "AAC" });
		aa2codon.put("H", new String[] { "CAT", "CAC" });
		aa2codon.put("E", new String[] { "GAA", "GAG" });
		aa2codon.put("D", new String[] { "GAT", "GAC" });
		aa2codon.put("K", new String[] { "AAA", "AAG" });
		aa2codon.put("R", new String[] { "CGT", "CGC", "CGA", "CGG", "AGA",
				"AGG" });
		aa2codon.put("STOP", new String[] { "TAA", "TAG", "TGA" });
		return aa2codon;
	}

	/**
	 * Concatenate gffs and fastas in a folder Outputs a new FASTA file and GFF
	 * file (master files) Replaces concatFASTA and concatGFF
	 * 
	 * @throws IOException
	 */
	public static void concatChromosome(String folder) throws IOException {
		int i;
		File inFolder = new File(folder);
		String[] files = inFolder.list();
		ArrayList<FASTAFile> fastas = new ArrayList<FASTAFile>();
		ArrayList<Strand[]> gffs = new ArrayList<Strand[]>();
		List<Strand> allStrands = new ArrayList<Strand>();

		FASTAFile master;

		// Read all files and parse GFFS and FASTAS, begin constructing offset
		// array
		for (i = 0; i < files.length; i++) {
			// Change these to backslashes for Windows
			if(!files[i].contains("master"))
			{
				if (files[i].matches(".*\\.gff")) {
					System.out.println(files[i] + " read as GFF");
					Strand current[] = readStrandsFromGFF(folder + File.separator + files[i]);
					gffs.add(current);
				} else if (files[i].matches(".*\\.fna")) {
					System.out.println(files[i] + " read as FASTA");
					FASTAFile current = readFastaStrand(new File(folder + File.separator
							+ files[i]));
					fastas.add(current);
				} else if(files[i].charAt(0) != '.'){
					System.err.println("Filetype not recognized:" + files[i]);
					System.err.println("To process, rename as a .gff or .fna");
	
					popupError("Filetype not recognized:"
							+ files[i]
							+ "\n\nTo process, rename as a .gff or .fna\n(continuing anyways)");
				}
			}
		}
		Collections.sort(fastas);

		int[] offsets = new int[fastas.size()];
		master = fastas.get(0);
		offsets[0] = master.length;
		for (i = 1; i < fastas.size(); i++) {
			// concatenate master and FASTAFile[i]
			master = concatenate(master, fastas.get(i));
			offsets[i] = master.offset;
			// System.out.println("offset[" + i + "] = " + master.offset);
		}

		for (Strand[] current : gffs) {
			Strand[] buffer = current;
			int offset = 0;

			// Find the matching FASTA to get the offset
			// Matching is based on name, i.e. "fosmid21" == "fosmid21"
			for (i = 0; i < fastas.size(); i++) {
				// found a match
				if (fastas.get(i).fosmid.equals(buffer[0].name)) {
					System.out.println(buffer[0].name + " matched");

					// calculate gff gene offset

					if (i > 0)
						offset = offsets[fastas.get(i).number - 1];

					// System.out.println(fastas.get(i). + " offset " + offset);

					for (Strand s : buffer) {
						s.start += offset;
						s.end += offset;
					}
					allStrands.addAll(Arrays.asList(buffer));
				}
			}

		}
		
		HashSet<Strand> hs = new HashSet<Strand>();
		hs.addAll(allStrands);
		allStrands.clear();
		allStrands.addAll(hs);
		
		Collections.sort(allStrands);
		

		String bigError = "";
		
		for(String s : errors)
		{
			bigError += s + "\n"; 
		}
		
		popupMessage(bigError);
		
		outputGFF(allStrands, folder + "master.gff");
		outputFastaFile(master, folder + "master.fna");
	}

	/**
	 * Puts two overalpping strands together
	 * 
	 * <pre>
	 * strands must be in order
	 * @param s1
	 * @param s2
	 * @return
	 */
	public static Strand concatenate(Strand a, Strand b) {
		String s1 = a.bases;
		String s2 = b.bases;

		String input = s2 + '$' + s1;

		int[] indices = computePi(input);
		int overlap = arrayMaxIndex(indices);

		String whole = s1
				+ s2.substring(indices[overlap]
						+ (indices.length - 1 - overlap));

		Strand both = new Strand();

		both.bases = whole;

		both.end = b.end;
		both.start = indices[overlap] + (indices.length - 1 - overlap);
		both.id = a.id;
		both.length = whole.length();
		both.direction = a.direction;
		both.name = a.name;
		both.transcriptId = a.transcriptId;
		both.type = a.type;

		return both;
	}

	/**
	 * Puts two overalpping strands together
	 * 
	 * <pre>
	 * strands must be in order
	 * @param s1 The first FASTA strand
	 * @param s2 The second FASTA strand
	 * @return A new FASTA
	 */
	public static FASTAFile concatenate(FASTAFile a, FASTAFile b) {
		String s1 = a.data;
		String s2 = b.data;

		if(errors == null)
			errors = new ArrayList<String>(); 
		
		String whole = s1 + 'x' + s2;

		String input = s2 + '$' + s1;

		int[] indices = computePi(input);
		int maxIndex = arrayMaxIndex(indices);
		// System.out.println(Arrays.toString(computePi(input)));

		int overlap = input.length() - maxIndex + indices[maxIndex];

		// System.out.println(whole);
		// System.out.println(indices[maxIndex] + (indices.length - 1 -
		// maxIndex));
		// System.out.println(indices[maxIndex]);

		int incorrect = 0;

		// System.out.println(overlap);

		// System.out.println(Arrays.toString(computePi(input)));

		int offset = indices[maxIndex] + (indices.length - 1 - maxIndex);

		try {
			whole = s1
					+ s2.substring(indices[maxIndex]
							+ (indices.length - 1 - maxIndex));

			for (int i = 0; i < overlap; i++) {
				// System.out.println(s2.charAt(i) + " " + s1.charAt(s1.length()
				// - overlap + i));

				if (s2.charAt(i) != s1.charAt(s1.length() - overlap + i))
					incorrect++;
			}
			System.out.println("Base pair errors in " + a.fosmid + "/" + b.fosmid + ": " +  incorrect);
			errors.add("Base pair errors in " + a.fosmid + "/" + b.fosmid + ": " +  incorrect); 

		} catch (StringIndexOutOfBoundsException ex) {
			offset = 1;
			System.err.println("Strings " + a.fosmid + " and " + b.fosmid
					+ " couldn't be matched - just adding it to the end");
		}

		System.out.println("Offset of " + b.fosmid + " is " + offset);

		FASTAFile both = new FASTAFile();

		both.data = whole;

		both.offset = a.length - offset;
		both.length = whole.length();
		both.name = a.name;
		both.fosmid = b.fosmid;
		both.number = b.number;

		// System.out.println("Returned offset is " + both.offset);

		return both;
	}

	/**
	 * DEPRECATED Concatenates all of the FASTA files in a folder into one big
	 * file
	 * 
	 * @param folder
	 *            The folder name that the FASTA files are in
	 * @return A FASTAFile of the entire folder's chromosome
	 * @throws FileNotFoundException
	 */
	public static FASTAFile concatFASTA(String folder)
			throws FileNotFoundException {
		int i;
		File inFolder = new File(folder);
		String filename[] = inFolder.list();
		FASTAFile[] parts = new FASTAFile[filename.length];
		FASTAFile master;

		for (i = 0; i < filename.length; i++) {
			File current = new File(folder + "/" + filename[i]);
			parts[i] = readFastaStrand(current);
		}

		List<FASTAFile> sortMe = Arrays.asList(parts);

		Collections.sort(sortMe);

		parts = sortMe.toArray(parts);

		master = parts[0];
		for (i = 1; i < filename.length; i++) {
			// concatenate FASTAFile[i] and FASTAFile[i + 1]
			master = concatenate(master, parts[i]);
		}

		return master;
	}

	/**
	 * DEPRECATED Reads all of the GFFs in a folder
	 * 
	 * @param folder
	 *            The folder name containing GFFs
	 * @return A list
	 * @throws FileNotFoundException
	 */
	public static List<Strand> concatGFF(String folder)
			throws FileNotFoundException {
		int i;
		List<Strand> parts = new ArrayList<Strand>();
		File inFolder = null;
		try {
			inFolder = new File(folder);
		} catch (Exception ex) {
			System.err.println("Folder " + folder + " not found.");
			return null;
		}

		String filename[] = inFolder.list();

		for (i = 0; i < filename.length; i++) {
			if (!filename[i].matches(".*\\.gff"))
				break;

			List<Strand> newstrands = Arrays.asList(readStrandsFromGFF(folder
					+ File.separator + filename[i]));
			for (Strand s : newstrands) {
				if (s != null && !parts.contains(s))
					parts.add(s);
			}
		}
		Collections.sort(parts);
		return parts;
	}

	/**
	 * Returns the index of the largest int in an int array
	 * 
	 * @param array
	 * @return
	 */
	public static int arrayMaxIndex(int[] array) {
		if (array.length == 0)
			return 0;
		int max = 0;
		for (int i = 0; i < array.length; i++) {
			if (array[i] > array[max])
				max = i;
		}

		return max;
	}

	/**
	 * Returns the index of the largest int in an int array
	 * 
	 * @param array
	 * @return
	 */
	public static int arrayMaxIndex(double[] array) {
		if (array.length == 0)
			return 0;
		int max = 0;
		for (int i = 0; i < array.length; i++) {
			if (array[i] > array[max])
				max = i;
		}

		return max;
	}

	// Taken from Professor Staley's 349 lecture
	public static int[] computePi(String s) {
		char[] pattern = s.toCharArray();

		int[] pi = new int[pattern.length + 1];
		int nextChr, pfxLen;
		pfxLen = 0;

		for (nextChr = 1; nextChr < pattern.length; nextChr++) {
			while (pattern[pfxLen] != pattern[nextChr] && pfxLen > 0)
				pfxLen = pi[pfxLen];

			if (pattern[pfxLen] == pattern[nextChr])
				pfxLen++;

			pi[nextChr + 1] = pfxLen;
		}
		return pi;
	}

	/**
	 * Reads a FASTA file and outputs the nucleotide sequence in the specified
	 * range
	 * 
	 * @param filename
	 * @param start
	 * @param end
	 * @return
	 */
	public static String readFASTA(String filename, int start, int end) {
		int i;
		StringBuilder fileContents = new StringBuilder();

		File gene;

		try {
			gene = new File(filename);
			Scanner sc = new Scanner(gene);

			sc.nextLine(); // ignore first line

			for (i = 0; i < start - 50; i += 50) {
				sc.nextLine();
			}

			// Add first line, from the first codon wanted to the end of the
			// line
			if (end - start > 50) {
				fileContents.append(sc.nextLine().substring(start - i));
			} else if (end - i < 50) {
				fileContents
						.append(sc.nextLine().substring(start - i, end - i));
				return fileContents.toString();
			}
			// Read the rest of the lines until the line before the end
			for (; i < end - 50; i += 50)
				fileContents.append(sc.nextLine());
			if (sc.hasNextLine())
				fileContents.append(sc.nextLine().substring(0, end - i));

			sc.close();
		} catch (FileNotFoundException ex) {
			popupError("There's something wrong with the way you typed "
					+ filename + " or something bad happened");
		}

		return fileContents.toString();
	}

	/*
	 * Read GFF file and returns array of strands
	 */
	public static Strand[] readStrandsFromGFF(String filename) {
		ArrayList<Strand> strandsAL = new ArrayList<Strand>();

		Strand[] strands = new Strand[0];

		try {
			File gene = new File(filename);
			Scanner sc = new Scanner(gene);

			while (sc.hasNextLine()) {
				Strand current = new Strand();

				try {
					current.name = sc.next();
					sc.next(); // dot
					current.type = sc.next(); // mRNA or CDS
					current.start = sc.nextInt();
					current.end = sc.nextInt();
					sc.next(); // dot
					current.direction = sc.next().charAt(0);
					sc.next(); // dot
					sc.next(); // "gene_id" hopefully doesn't change
					current.id = sc.next();
					sc.next(); // transcript_id
					current.transcriptId = sc.nextLine();

					current.id = current.id.substring(1,
							current.id.length() - 2);
					current.transcriptId = current.transcriptId.substring(1,
							current.transcriptId.length() - 2);

					current.length = current.end - current.start;
					if (current.end < current.start) {
						int temp = current.start;
						current.start = current.end;
						current.end = temp;
					}
					strandsAL.add(current);
				} catch (Exception ex) {
					popupError("Problem reading GFF file" + ex);
				}

			}
			sc.close();
		} catch (Exception ex) {
			popupError(ex + " - Invalid Strand filename?");
		}

		strands = strandsAL.toArray(strands);

		return strands;
	}

	public static FASTAFile readFastaStrand(File fasta)
			throws FileNotFoundException {
		FASTAFile out = new FASTAFile();
		Scanner sc = new Scanner(fasta);

		out.name = sc.next().substring(1);

		String info = sc.next();
		out.fosmid = info.substring(info.indexOf('=') + 1, info.indexOf(':'));
		out.number = Integer.parseInt(out.fosmid.substring("contig".length()));
		out.length = Integer.parseInt(info.substring(info.indexOf('-') + 1)
				.trim());

		while (sc.hasNext()) {
			String temp = sc.next();
			if (out.data != null)
				out.data += temp;
			else
				out.data = temp;
		}

		return out;
	}

	/**
	 * Takes a list of GFFs and outputs one big GFF file
	 * 
	 * @param gffs
	 * @param outputname
	 *            the name of the file to write. Add the extension, this method
	 *            doesn't.
	 * @throws IOException
	 */
	public static void outputGFF(List<Strand> gffs, String outputname)
			throws IOException {
		PrintWriter out = new PrintWriter(new FileWriter(outputname));
		for (Strand s : gffs) {
			out.printf(
					"%s\t.\t%s\t%d\t%d\t.\t%c\t.\tgene_id \"%s\"; transcript_id \"%s\";\n",
					"master", s.type, s.start, s.end, s.direction, s.id,
					s.transcriptId);
		}
		out.close();
		System.out.println("GFF file " + outputname + " written");
		popupMessage("GFF file " + outputname + " written");
		
	}

	// Outputs a fasta file with the name of the input
	public static void outputFastaFile(FASTAFile input, String name)
			throws IOException {
		try {
			PrintWriter fna = new PrintWriter(new FileWriter(name));
			fna.printf(">%s range=%s:1-%d    \n", input.name, input.fosmid,
					input.length);

			for (int i = 0; i < input.length; i += 50) {
				fna.print("      ");
				fna.println(input.data.substring(i,
						Math.min((i + 50), input.data.length())));
			}
			System.out.println("FASTA file " + name + " written");
			popupMessage("FASTA file " + name + " written");
			fna.close();
		
			
		} catch (Exception ex) {
			ex.printStackTrace();
			System.err.println("FASTA file " + name + " NOT written");
			popupError("FASTA file " + name + " NOT written"); 
		}

		
	}

	public static void printGCtoCSV(int frame, int shift, String gffFile,
			String csvFile) throws IOException {

		PrintWriter csv;
		csv = new PrintWriter(new FileWriter(csvFile));

		csv.println("Start Codon,End Codon,GC Percentage");

		Strand[] files = readStrandsFromGFF(gffFile);

		for (Strand gene : files) {
			Strand[] frames;
			// System.out.println(gene.name + " " + gene.start + " " +
			// gene.end);

			String input = readFASTA(gene.name + ".txt", gene.start, gene.end);
			gene.bases = input;
			frames = gcPercentage(gene, frame, shift);

			for (Strand section : frames) {
				csv.printf("%d,%d,%.1f%%\n", gene.start + section.start,
						gene.start + section.end, section.gcPercent);

				// System.out.printf("%d,%d,%.1f%%\n", gene.start +
				// section.start, gene.start + section.end, section.gcPercent);
			}
			csv.println();
		}
		csv.close();

		// System.out.println("Done! File " +filename.substring(0,
		// filename.lastIndexOf('.')) + ".csv" + " written.");
		// System.out.println("Done. Goodbye.");
	}

	/*
	 * Pops up an error Message
	 */
	public static void popupError(String message) {
		JOptionPane.showMessageDialog(null, message, "Error",
				JOptionPane.ERROR_MESSAGE);
	}

	/*
	 * Convenience method for a warning popup dialog
	 */
	public static void popupMessage(String message) {
		JOptionPane.showMessageDialog(null, message, "",
				JOptionPane.PLAIN_MESSAGE);
	}

	// Representation of a FASTA file, basically a long string of bases
	public static class FASTAFile implements Comparable<FASTAFile> {
		String name;
		int number;
		String fosmid;
		int length;
		String data;
		int offset;

		public int compareTo(FASTAFile arg0) {
			return this.number - arg0.number;
		}

	}

	/**
	 * Class representation of a strand of DNA
	 * 
	 * @author johntbiddle
	 * 
	 */

	public static class Strand implements Comparable<Strand> {
		String bases;
		// Name of the file this strand came from
		String name;
		// number of nucleotides long
		int length;
		// start position on gene
		int start;
		// end position on gene
		int end;
		// mRNA or CDS
		String type;
		// + or -
		char direction;
		// GC percentage if calculated
		double gcPercent;
		// ID of this gene
		String id;
		// Transcript ID
		String transcriptId;

		public Strand() {
		}

		public Strand(int start, int end, char direction, String bases) {
			this.start = start;
			this.end = end;
			this.direction = direction;
			this.length = end - start;
			this.bases = bases;

			if (this.direction != '+' && this.direction != '-')
				popupError("Strand made with odd direction " + direction);
		}

		public Strand(String filename, int start, int end,
				boolean reverseComplement) {
			int i;
			StringBuilder fileContents = new StringBuilder();

			this.start = start;
			this.end = end;
			this.length = end - start;
			this.name = filename.substring(0, filename.indexOf("."));

			File gene;

			try {
				gene = new File(filename);
				Scanner sc = new Scanner(gene);

				sc.nextLine(); // ignore first line

				for (i = 0; i < start - 50; i += 50) {
					sc.nextLine();
				}

				// Add first line, from the first codon wanted to the end of the
				// line
				if (end - start > 50) {
					fileContents.append(sc.nextLine().substring(start - i));
				} else if (end - i < 50) {
					fileContents.append(sc.nextLine().substring(start - i,
							end - i));

					return;
				}
				// Read the rest of the lines until the line before the end
				else {
					for (; i < end - 50; i += 50)
						fileContents.append(sc.nextLine());

					if (sc.hasNextLine())
						fileContents
								.append(sc.nextLine().substring(0, end - i));
				}
				sc.close();
			} catch (FileNotFoundException ex) {
				popupError("There's something wrong with the way you typed "
						+ filename + " or something bad happened");
			}
			this.bases = fileContents.toString();
		}

		public boolean equals(Object other) {
			if (other == null || other.getClass() != this.getClass())
				return false;

			return (this.start == ((Strand) other).start
					&& this.end == ((Strand) other).end
					&& this.direction == ((Strand) other).direction && this.transcriptId
						.equals(((Strand) other).transcriptId));

		}

		public int compareTo(Strand o) {
			// int number =
			// Integer.parseInt(this.name.substring("contig".length()))
			// - Integer.parseInt(o.name.substring("contig".length()));
			// if(number ==0)
			return this.start - o.start;
			// return number;
		}
	}

}
m