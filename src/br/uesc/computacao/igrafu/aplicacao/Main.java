package br.uesc.computacao.igrafu.aplicacao;

import br.uesc.computacao.igrafu.apresentacao.Igrafu_UI;
import javax.swing.UIManager;
import javax.swing.UnsupportedLookAndFeelException;

/**
 *Classe que inicia o programa.
 * Configura Look and Feel e exibe a janela inicial.
 *
 * @author Anderson Carlos Sousa e Santos </br>
 * Orientadora: Martha Ximena Torres Delgado
 */
public class Main {

    /**
     * @param args the command line arguments
     */
    public static void main(String args[]) {

        java.awt.EventQueue.invokeLater(new Runnable() {
            public void run() {
                try {
                    UIManager.setLookAndFeel(
                            "com.sun.java.swing.plaf.gtk.GTKLookAndFeel");

		} catch (ClassNotFoundException e) {
			e.printStackTrace();

		} catch (InstantiationException e) {
			e.printStackTrace();
		} catch (IllegalAccessException e) {
			e.printStackTrace();
		} catch (UnsupportedLookAndFeelException e) {
			e.printStackTrace();
		}

                new Igrafu_UI().setVisible(true);

            }
        });
    }
}
