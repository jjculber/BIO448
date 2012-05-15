package bioGUI;

import bioGUI.dialogs.CGContentDialog;
import bioGUI.dialogs.CodonBiasDialog;
import bioGUI.dialogs.GffConcatenatorDialog;
import bioGUI.dialogs.PalindromeDialog;
import bioGUI.dialogs.RepeatDialog;

import java.awt.event.ActionListener;

import javax.swing.JFrame;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;

public class MainWindow extends JFrame {
   private static final String WINDOW_TITLE = "CSC 448 - Bioinformatics";
   private static final int FRAME_WIDTH = 500,
                            FRAME_HEIGHT = 600;

   private static MainWindow mMainFrame = null;
   private JMenuBar mMenuBar;

   /**
    * This is a singleton constructor for the MainWindow. This means that
    * only one MainWindow is instantiated at a time. Whenever a component
    * wants to get the main window, it'll always get this one.
    */
   public static MainWindow getMainFrame() {
      if (mMainFrame == null) {
         mMainFrame = new MainWindow();
      }

      return mMainFrame;
   }
   
   private MainWindow() {
      super(WINDOW_TITLE);

      setSize(FRAME_WIDTH, FRAME_HEIGHT);
      //centers the component on the screen
      setLocationRelativeTo(null);
   }

   /**
    * Initialize the Menu Bar for this JFrame and any other relevant
    * components.
    */
   public void init() {
      mMenuBar = new JMenuBar();

      mMenuBar.add(initFileMenu());
      mMenuBar.add(initAnalysisMenu());

      setJMenuBar(mMenuBar);
      validate();
   }

   /**
    * Initialize menu items that will be present in the 'File' menu, then
    * return the initialized 'File' menu.
    */
   private JMenu initFileMenu() {
      JMenuItem exitProgram = new JMenuItem("Exit");
      exitProgram.addActionListener(new ActionListener() {
         public void actionPerformed(java.awt.event.ActionEvent e) {
            System.exit(0);
         }
      });


      JMenu fileMenu = new JMenu("File");
      fileMenu.add(exitProgram);

      return fileMenu;
   }

   /**
    * Initialize menu items that will be present in the 'Analysis' menu, then
    * return the initialized 'Analysis' menu.
    */
   private JMenu initAnalysisMenu() {
      JMenuItem gcContent = new JMenuItem("Calculate GC Content");
      gcContent.addActionListener(new ActionListener() {
         public void actionPerformed(java.awt.event.ActionEvent e) {
            CGContentDialog gcContentDialog = new CGContentDialog(mMainFrame, "GC Content Input Parameters");

            gcContentDialog.init();
            gcContentDialog.setVisible(true);
         }
      });

      JMenuItem condonBias = new JMenuItem("Calculate Codon Bias");
      condonBias.addActionListener(new ActionListener() {
         public void actionPerformed(java.awt.event.ActionEvent e) {
            CodonBiasDialog codonBiasDialog = new CodonBiasDialog(mMainFrame, "Codon Bias Input Parameters");

            codonBiasDialog.init();
            codonBiasDialog.setVisible(true);
         }
      });

      JMenuItem gffConcatenator = new JMenuItem("Concatenate GFF Files");
      gffConcatenator.addActionListener(new ActionListener() {
         public void actionPerformed(java.awt.event.ActionEvent e) {
        	 GffConcatenatorDialog gffConcatenatorDialog = new GffConcatenatorDialog(mMainFrame, "GFF Concatenator Input Parameters");

            gffConcatenatorDialog.init();
            gffConcatenatorDialog.setVisible(true);
         }
      });

      JMenuItem repeats = new JMenuItem("Find Repeats");
      repeats.addActionListener(new ActionListener() {
         public void actionPerformed(java.awt.event.ActionEvent e) {
        	 RepeatDialog repeatDialog = new RepeatDialog(mMainFrame, "Find Repeats");

        	 repeatDialog.init();
        	 repeatDialog.setVisible(true);
         }
      });

      JMenuItem palindromes = new JMenuItem("Find Palindromes");
      palindromes.addActionListener(new ActionListener() {
         public void actionPerformed(java.awt.event.ActionEvent e) {
        	 PalindromeDialog palindromeDialog = new PalindromeDialog(mMainFrame, "Find Palindromes");

        	 palindromeDialog.init();
        	 palindromeDialog.setVisible(true);
         }
      });
      
      
      JMenu analysisMenu = new JMenu("Analysis");
      analysisMenu.add(gcContent);
      analysisMenu.add(condonBias);
      analysisMenu.add(gffConcatenator);
      analysisMenu.add(repeats);
      analysisMenu.add(palindromes);

      return analysisMenu;
   }

   /**
    * Convenience method for ensuring the MainWindow is visible and repainted.
    */
   public void showWindow() {
      setVisible(true);
      repaint();
   }
}
