package bioGUI.dialogs;

import java.awt.Component;
import java.awt.ComponentOrientation;
import java.awt.Container;
import java.awt.Frame;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.FlowLayout;

import javax.swing.JPanel;
import javax.swing.JOptionPane;
import javax.swing.JDialog;
import javax.swing.JFileChooser;
import javax.swing.JLabel;
import javax.swing.JTextField;
import javax.swing.JButton;
import javax.swing.BoxLayout;

import bioGUI.model.DNALibrary;


import java.io.File;

public class CodonBiasDialog extends JDialog {
   /*
    * CONSTANTS
    */
   private final int DIALOG_HEIGHT = 200, DIALOG_WIDTH = 500;

   /*
    * GUI Components
    */
   private Container mPane = null, mOwner = null;
   private JDialog mDialog = null;
   private JTextField mGFF;

   public CodonBiasDialog(Frame owner, String title) {
      super(owner, title);

      this.setSize(DIALOG_WIDTH, DIALOG_HEIGHT);
      this.setResizable(false);
      this.setLocationRelativeTo(null);

      mOwner = owner;
      mDialog = this;

      mPane = this.getContentPane();
      mPane.setLayout(new BoxLayout(mPane, BoxLayout.Y_AXIS));
      mPane.setSize(DIALOG_WIDTH, DIALOG_HEIGHT);

      setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);

      mGFF = new JTextField(20);
   }

   public static void main(String[] args) {
      CodonBiasDialog dialog = new CodonBiasDialog(null, "Codon Bias");

      dialog.init();
      dialog.setVisible(true);
   }

   /**
    * Method for initializing a default Dialog window ready for GC Content
    * input
    */
   public void init() {
	  JLabel gffFileLabel1 = new JLabel("(required FASTA files must be in the same directory)");
      JLabel gffFileLabel2 = new JLabel("Select GFF File:");
      JPanel gffFileField = prepareGffField(mGFF);

      mPane.add(gffFileLabel1);
      mPane.add(gffFileLabel2);
      mPane.add(gffFileField);

      mPane.add(initControls());

      mPane.validate();
   }

   /**
    * Convenience method for constructing a JPanel that contains the JTextField
    * and file browse button used for selecting a FASTA file.
    */
   private JPanel prepareGffField(JTextField fastaField) {
      JPanel fastaFileField = new JPanel();

      fastaFileField.setLayout(new FlowLayout(FlowLayout.LEADING));

      fastaFileField.add(fastaField);
      fastaFileField.add(prepareBrowseButton(fastaField));

      return fastaFileField;
   }

   /**
    * Convenience method for creating a file browse button. This is abstracted
    * so that it is not necessarily associated with the fasta file field.
    */
   private JButton prepareBrowseButton(final JTextField fastaField) {
      JButton fileBrowse = new JButton("Browse");

      fileBrowse.addActionListener(new ActionListener() {

         public void actionPerformed(ActionEvent e) {
            JFileChooser chooser = new JFileChooser();
            int returnVal = chooser.showOpenDialog(chooser);
            
            if (returnVal == JFileChooser.CANCEL_OPTION) {
               System.out.println("cancelled");
            }

            else if (returnVal == JFileChooser.APPROVE_OPTION) {
               File fastaFile = chooser.getSelectedFile();
               fastaField.setText(fastaFile.getAbsolutePath());
            }

            else {
               System.out.println("Encountered Unknown Error");
               System.exit(0);
            }
         }
      });

      return fileBrowse;
   }

   /**
    * Convenience method for creating JTextFields for the start or end
    * nucleotide range parameters. This is abstracted so that it may be called
    * for both the start and end text fields.
    */
   @SuppressWarnings("unused")
   private JPanel prepareRangeField(String labelPrefix, JTextField rangeInput) {
      JPanel rangeField = new JPanel();

      rangeField.setLayout(new BoxLayout(rangeField, BoxLayout.Y_AXIS));

      rangeField.add(new JLabel(labelPrefix + " Position:"));
      rangeField.add(rangeInput);

      return rangeField;
   }

   /**
    * Convenience method for creating a JPanel containing an okay and cancel
    * button at the bottom of the Dialog window. If cancel is pressed then the
    * dialog window is disposed and nothing happens. If okay is pressed then
    * the appropriate input parameters should be grabbed and handled
    * appropriately. Some slight error checking is also implemented.
    */
   public JPanel initControls() {
      JPanel dialogControls = new JPanel();

      dialogControls.setComponentOrientation(ComponentOrientation.LEFT_TO_RIGHT);

      dialogControls.add(createOkayButton());
      dialogControls.add(createCancelButton());

      dialogControls.setAlignmentX(Component.CENTER_ALIGNMENT);

      return dialogControls;
   }

   /**
    * Convenience method for creating the okay button at the bottom of the
    * Dialog window. Anything that should be added to the okay button's
    * functionality should probably be added here, or separated into an
    * ActionListener and added as an ActionListener for the button.
    */
   private JButton createOkayButton() {
      JButton okayButton = new JButton("Okay");

      okayButton.addActionListener(new ActionListener(){
         public void actionPerformed(ActionEvent e) {
            if (mGFF.getText().equals("")) {
               JOptionPane.showMessageDialog(mOwner,
                "No GFF file was selected",
                "Invalid File", JOptionPane.ERROR_MESSAGE);
               return;
            }
            else
            {
              // TODO Codon bias work
            	DNALibrary.calcCodonBias(mGFF.getText());
            }
            dispose();
         }
      });

      return okayButton;
   }

   /**
    * Convenience method for creating the cancel button at the bottom of the
    * Dialog window. This currently just disposes of the dialog window, though
    * some more complex behavior may be desired.
    */
   private JButton createCancelButton() {
      JButton cancelButton = new JButton("Cancel");

      cancelButton.addActionListener(new ActionListener() {
         public void actionPerformed(ActionEvent e) {
            dispose();
            return;
         }
      });

      return cancelButton;
   }
}
