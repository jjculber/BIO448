package bioGUI;

import bioGUI.MainWindow;

public class BioGUI {

   public static void main(String[] args) {
      MainWindow gui = MainWindow.getMainFrame();
      gui.init();
      gui.showWindow();
   }
}
