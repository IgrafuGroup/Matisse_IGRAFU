package br.uesc.computacao.igrafu.aplicacao;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import javax.swing.JButton;
import javax.swing.JLabel;
import javax.swing.JProgressBar;
import javax.swing.SwingUtilities;

/**
 * Classe que realiza execução do MrBayes, tratando dados de saida para serem
 * exibidos na GUI.
 *
 * @author Anderson Carlos Sousa e Santos </br>
 * Orientadora: Martha Ximena Torres Delgado
 */
public class ExecucaoMrBayes implements Runnable{

    private int ngen;
    private int count = 1;
    private JLabel estimateLabel;
    private JLabel outputLabel1;
    private JLabel outputLabel2;
    private JProgressBar progressBar;
    private String estimative;
    private String output = "";
    private String last_output = "";
    private JButton cancel_okButton;
    private Process p;


    /**
     * 
     * @param ngen numero de ciclos
     * @param estimate Rotulo para exibir a estimativa  de termino da execução
     * @param output1 Rotulo para exibir estatisticas de convergencia
     * @param output2 Rotulo para exibir estatisticas de convergencia </br>
     * (exibe a convergência que foi exibida anteriormente para comparação)
     * @param progress - barra de progresso para ser incrementada.
     * @param button butão que terá texto alterado para representar fim da execução;
     * 
     */
    public ExecucaoMrBayes(int ngen, JLabel estimate, JLabel output1,
                            JLabel output2, JProgressBar progress, JButton button){
        this.ngen = ngen;
        this.estimateLabel = estimate;
        this.outputLabel1 = output1;
        this.outputLabel2 = output2;
        this.progressBar = progress;
        this.cancel_okButton = button;
    }


    /**
     * Metodo herdado de Runnable, para execução de uma nova thread.
     * <p> Inicia execução do programa externo mrbayes-3.1.2 com o "batch.nex"</br>
     * como parametro de entrada. Lê a saida do programa para obter atraves </br>
     * de strings o tempo estimado para encerramento do mcmc e a quantidade de </br>
     * ciclos realizados.</p>
     *
     * <p> Utiliza-se de SwingUtilities.invokeLater(Runnable) para atualizar </br>
     * a barra de progresso (JProgressBar) e e os rotulos que exibem as estatisticas</p>
     *
     * @see java.lang.Runnable
     * @see javax.swing.SwingUtilities
     * @see br.uesc.computacao.igrafu.apresentacao.MrBayesUI#jButtonExecuteActionPerformed(java.awt.event.ActionEvent) 
     */
    public void run() {

        try
        {
            String line = "";

            String tRemain = "not";
            int index;
            
            p = Runtime.getRuntime().exec ("./Programs/mrbayes-3.1.2/mb batch.nex");

            BufferedReader input =
            new BufferedReader(new InputStreamReader(p.getInputStream()));

            while ((line = input.readLine()) != null) {
                System.out.println(line);

                //if( (index = line.indexOf("--", line.indexOf("]"))) != -1 ){
                if(line.matches("[\\s]*[0-9]+[\\s][-][-].+[*].+[-][-].+")){
                    //tRemain = line.substring(index+1);
                    tRemain = line.split("--")[2];
                    count++;
                    
                    estimative = "Completed iteration " + count + " of " + ngen +
                             " (Estimated " + tRemain + " remaining)";

                    SwingUtilities.invokeLater
                    (
                        new Runnable()
                        {
                            public void run()
                            {
                                estimateLabel.setText(estimative);
                                
                                progressBar.setValue( ((count*100) / ngen) );
                                
                            }
                        }
                    );
                }
                if( (index = line.indexOf("Average")) != -1 ){
                    last_output = output;
                    output = line;

                    SwingUtilities.invokeLater
                    (
                        new Runnable()
                        {
                            public void run()
                            {
                                outputLabel1.setText(output);
                                outputLabel2.setText(last_output);
                            }
                        }
                    );
                    
                }
            }
        input.close();
        cancel_okButton.setText("OK");
        } catch (Exception err) {
            err.printStackTrace();
        }
    }

    
    /**
     * Interrompe a execução a aplicação externa, matando o subprocesso.
     */
    public void interroper(){
        p.destroy();

    }
}
