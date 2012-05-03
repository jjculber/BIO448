package bioGUI.dialogs;

import java.awt.Component;
import java.awt.ComponentOrientation;
import java.awt.Container;
import java.awt.FlowLayout;
import java.awt.Frame;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;

import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JDialog;
import javax.swing.JFileChooser;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JTextField;

import bioGUI.model.DNALibrary;

public class CGContentDialog extends JDialog {
   /*
    * CONSTANTS
    */
   private final int DIALOG_HEIGHT = 450, DIALOG_WIDTH = 500;

   /*
    * GUI Components
    */
   private Container mPane = null, mOwner = null;
   private JDialog mDialog = null;
   private JTextField mFasta, mRangeBegin, mRangeEnd, mFrameSize, mFrameShift;

   public CGContentDialog(Frame owner, String title) {
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

      mFasta = new JTextField(20);
      mRangeBegin = new JTextField(20);
      mRangeEnd = new JTextField(20);
      
      mFrameSize = new JTextField(20);
      mFrameShift = new JTextField(20);
   }

   public static void main(String[] args) {
      CGContentDialog dialog = new CGContentDialog(null, "GC Content");

      dialog.init();
      dialog.setVisible(true);
   }

   /**
    * Method for initializing a default Dialog window ready for GC Content
    * input
    */
   public void init() {
      JLabel fastaFileLabel = new JLabel("Select FASTA File:");
      JPanel fastaFileField = prepareFastaField(mFasta);

      JPanel beginField = prepareRangeField("Start Position:", mRangeBegin);
      JPanel endField = prepareRangeField("End Position:", mRangeEnd);

      JPanel nucleotideRangeField = new JPanel();
      nucleotideRangeField.setLayout(new FlowLayout(FlowLayout.LEADING));
      nucleotideRangeField.add(beginField);
      nucleotideRangeField.add(endField);

      JPanel frameSize = prepareRangeField("Frame Size:", mFrameSize);
      JPanel frameShift = prepareRangeField("Frame Shift:", mFrameShift);
      
      JPanel frameField = new JPanel();
      frameField.setLayout(new FlowLayout(FlowLayout.LEADING));
      frameField.add(frameSize);
      frameField.add(frameShift);
      

      mPane.add(fastaFileLabel);
      mPane.add(fastaFileField);

      mPane.add(nucleotideRangeField);
      
      mPane.add(frameField);

      mPane.add(initControls());

      mPane.validate();
   }

   /**
    * Convenience method for constructing a JPanel that contains the JTextField
    * and file browse button used for selecting a FASTA file.
    */
   private JPanel prepareFastaField(JTextField fastaField) {
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
   private JPanel prepareRangeField(String labelPrefix, JTextField rangeInput) {
      JPanel rangeField = new JPanel();

      rangeField.setLayout(new BoxLayout(rangeField, BoxLayout.Y_AXIS));

      rangeField.add(new JLabel(labelPrefix));
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
            if (mFasta.getText().equals("")) {
               JOptionPane.showMessageDialog(mOwner,
                "No FASTA file was selected",
                "Invalid File", JOptionPane.ERROR_MESSAGE);
               return;
            }
            else
            {
            	int start, end, size, shift;
            	
            	try {
            		start = Integer.parseInt(mRangeBegin.getText());
            	} catch (Exception ex) {
            		start = 0;
            	}
            	
            	try {
            		end = Integer.parseInt(mRangeEnd.getText());
            	} catch (Exception ex) {
            		end = 0;
            	}

            	try {
            		size = Integer.parseInt(mFrameSize.getText());
            		shift = Integer.parseInt(mFrameShift.getText());
            		
                	DNALibrary.calcGCContent(mFasta.getText(), start, end, size, shift);

            	} catch (Exception ex) {
            		DNALibrary.popupError("Bad Frame size or shift.");
            	}
            	
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
