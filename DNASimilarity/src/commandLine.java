import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.StringReader;
import java.util.Scanner;

import bioGUI.model.DNALibrary;

import neobio.alignment.*;

public class commandLine extends DNALibrary{

	public static ScoringMatrix score;
	
	public static void main(String[] args) throws FileNotFoundException, IOException, InvalidScoringMatrixException, InvalidSequenceException, IncompatibleScoringSchemeException
	{
		Scanner sc = new Scanner(System.in);
		System.out.println("What gap penalty would you like?");
		score = new ScoringMatrix(new FileReader("blosum62.txt"));
		
		score.deletionCost = sc.nextInt();
		
		System.out.println("FASTA candidate 1?");
		String a = sc.next();
		
		System.out.println("FASTA candidate 2?");
		String b = sc.next(); 
		
		System.out.println(getScore(a, b));
		//getBestAlignment(ref, a, b); 
		
	}
	
	public static int getScore(String a, String b) throws IOException, InvalidSequenceException, IncompatibleScoringSchemeException
	{
		final NeedlemanWunsch algorithm = new NeedlemanWunsch(); 
		
		algorithm.loadSequences(new StringReader(a), new StringReader(b)); 
		
		algorithm.setScoringScheme(score);
		
		PairwiseAlignment alignment = algorithm.getPairwiseAlignment();
		System.out.println(alignment.getGappedSequence1());
		System.out.println(alignment.getGappedSequence2());
		
		return alignment.getScore();		
	}
	
	/**
	 * returns positive score if a is better
	 * returns negative score if b is better
	 * @param ref
	 * @param a
	 * @param b
	 * @return
	 * @throws IOException
	 * @throws InvalidSequenceException
	 * @throws IncompatibleScoringSchemeException
	 */
	public static int getBestAlignment(FASTAFile ref, FASTAFile a, FASTAFile b) throws IOException, InvalidSequenceException, IncompatibleScoringSchemeException
	{
		int scoreA = getScore(ref.data, a.data);
		int scoreB = getScore(ref.data, b.data);
		
		if(scoreA > scoreB)
		{
			System.out.println(a.fosmid + " is more accurate with score " + scoreA  );
			return scoreA;
		}
		else System.out.println(b.fosmid + " is more accurate with score " + scoreB  );
		return (-1) * scoreB; 
	}
}
