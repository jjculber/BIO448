package bioGUI.model;

import java.util.Scanner;
import java.io.File;
import java.lang.StringBuilder;
import java.awt.BorderLayout;
import java.awt.Container;
import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;

public class CGContent {
   public static File inputFile;

	public static void main(String[] args) throws Exception {

			System.out.println(args[1] + args[2]);

			System.out.printf("%.1f%%\n", gcSubCount(args[0], Integer.valueOf(args[1]), Integer.valueOf(args[2])));
	}

	public static double gcSubCount(String filename, int start, int end) {
		long total= 0;
		StringBuilder fileContents = new StringBuilder();

		try
		{
			Scanner sc = new Scanner(new File(filename));
			
			sc.nextLine(); //ignore first line
			
			while (sc.hasNextLine()) {
				fileContents.append(sc.nextLine());
			}
		}
		catch(Exception ex)
		{
			System.out.println("Invalid filename");
		}

		String sub = fileContents.toString().substring(start, end);
		
		
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
