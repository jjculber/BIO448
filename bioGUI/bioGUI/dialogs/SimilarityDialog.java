package bioGUI.dialogs;

import java.awt.Component;
import java.awt.ComponentOrientation;
import java.awt.Container;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.Frame;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.io.FileReader;

import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JComponent;
import javax.swing.JDialog;
import javax.swing.JFileChooser;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.JTextField;

import neobio.alignment.ScoringMatrix;
import bioGUI.model.DNALibrary;
import bioGUI.model.Similarity;

public class SimilarityDialog extends JDialog {
	/*
	 * CONSTANTS
	 */
	private final int DIALOG_HEIGHT = 450, DIALOG_WIDTH = 550;

	/*
	 * GUI Components
	 */
	private Container mPane = null, mOwner = null;
	private JDialog mDialog = null;

	private JTextArea mDNA1, mDNA2, mResults;
	private JScrollPane mDNA1Scroll, mDNA2Scroll, mResultsScroll;

	public SimilarityDialog(Frame owner, String title) {
		super(owner, title);

		this.setSize(DIALOG_WIDTH, DIALOG_HEIGHT);
		// this.setResizable(false);
		this.setLocationRelativeTo(null);

		mOwner = owner;
		mDialog = this;

		mPane = this.getContentPane();
		mPane.setLayout(new BoxLayout(mPane, BoxLayout.Y_AXIS));
		mPane.setSize(DIALOG_WIDTH, DIALOG_HEIGHT);

		setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);

		mDNA1 = new JTextArea(100, 100);
		mDNA2 = new JTextArea(100, 100);
		mResults = new JTextArea(100, 100);
		mResults.setEnabled(false);
		mResults.setLineWrap(false);
		mDNA1.setLineWrap(true);
		mDNA2.setLineWrap(true);
		Font mono = new Font(Font.MONOSPACED, Font.PLAIN, 12);
		mResults.setFont(mono);

		mDNA1Scroll = new JScrollPane(mDNA1, JScrollPane.VERTICAL_SCROLLBAR_ALWAYS, JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
		mDNA2Scroll = new JScrollPane(mDNA2, JScrollPane.VERTICAL_SCROLLBAR_ALWAYS, JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
		mResultsScroll = new JScrollPane(mResults, JScrollPane.VERTICAL_SCROLLBAR_NEVER, JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);
	}

	public static void main(String[] args) {
		SimilarityDialog dialog = new SimilarityDialog(null, "DNA Similarity");

		dialog.init();
		dialog.setVisible(true);
	}

	/**
	 * Method for initializing a default Dialog window ready for GC Content
	 * input
	 */
	public void init() {
		JLabel feild1Label = new JLabel("Enter Amino Acid String:");
		JLabel feild2Label = new JLabel("Enter Amino Acid String:");
		JLabel resultsLabel = new JLabel("Results:");

		mPane.add(feild1Label);
		mPane.add(mDNA1Scroll);

		mPane.add(feild2Label);
		mPane.add(mDNA2Scroll);
		
		mPane.add(resultsLabel);
		mPane.add(mResultsScroll);

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
	private JPanel prepareRangeField(String labelPrefix, JComponent rangeInput) {
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

		dialogControls
				.setComponentOrientation(ComponentOrientation.LEFT_TO_RIGHT);

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

		okayButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				if (mDNA1.getText().equals("") || mDNA2.getText().equals("")) {
					JOptionPane.showMessageDialog(mOwner,
							"No DNA file was input", "Invalid Input",
							JOptionPane.ERROR_MESSAGE);
					return;
				} else {

					// TODO call similarity code
					try {
					ScoringMatrix score = new ScoringMatrix(new FileReader("blosum62.txt"));
					
					score.deletionCost = -3; // TODO change to user input later
					
					String a = mDNA1.getText();
					String b = mDNA2.getText(); 
					
					String results = Similarity.computeGlobalAlignment(a, b);
							
					mResults.setText(results);

					} catch (Exception ex) {
						DNALibrary.popupError(ex.getMessage());
					}
				}
//				dispose();
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
		JButton cancelButton = new JButton("Close");

		cancelButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				dispose();
				return;
			}
		});

		return cancelButton;
	}
	
}
