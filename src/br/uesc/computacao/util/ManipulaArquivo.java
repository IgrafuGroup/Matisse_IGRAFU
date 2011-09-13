package br.uesc.computacao.util;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.RandomAccessFile;

import javax.swing.JOptionPane;

/**
 * Classe utilizada para manipular arquivos
 * @author Zilton José Maciel Cordeiro Junior - Orientadora: Martha Ximena Torres Delgado
 * @version 2.0
 */

public class ManipulaArquivo {
    private String arquivo;
    private String str;
    private static String conteudoDoArquivo = null;

    /**
     * 
     * @param arquivo
     * @param str
     */
    public ManipulaArquivo(String arquivo, String str) {
        this.arquivo = arquivo;
        this.str = str;
        gravaArquivo(this.arquivo, this.str);
    }

    /**
     *
     * @param arquivo
     */
    public ManipulaArquivo(String arquivo) {
    	leArquivo(arquivo);
    }

    /**
     *
     */
    public ManipulaArquivo(){};

    /**
     * grava str na proxima linha do 'arquivo'
     * @param arquivo
     * @param str
     */
    public static void gravaArquivo( String arquivo, String str ) {
        File file = new File( arquivo );
        try {
            FileWriter writer = new FileWriter( file, true );
            PrintWriter saida = new PrintWriter( writer, true );
            saida.println(str);
            saida.close();
            writer.close();
         }
         catch( IOException ioe ) {
             ioe.printStackTrace();
             JOptionPane.showMessageDialog( null, " Erro - Não foi possível gravar no arquivo ", "Erro", JOptionPane.ERROR_MESSAGE );
         }
         catch( Exception e ) {
             e.printStackTrace();
             JOptionPane.showMessageDialog( null, " Erro - Não foi possível gravar no arquivo ", "Erro", JOptionPane.ERROR_MESSAGE );
         }
    }

    /**
     * retorna todo o conteudo do 'arquivo'
     * @param arquivo
     * @return
     */
    public static int contaLinhas( String arquivo ) {
    	int cont = 0;
        try {
            BufferedReader in = new BufferedReader( new FileReader( arquivo ) );
            while( (in.readLine()) != null )
                cont++;
            in.close();
        }
        catch( IOException ioe ) {
            ioe.printStackTrace();
            JOptionPane.showMessageDialog( null, " Erro - Não foi possível ler o arquivo '" + arquivo +"' ", "Erro ao ler arquivo", JOptionPane.ERROR_MESSAGE );
        }
        catch( Exception e ) {
            e.printStackTrace();
            JOptionPane.showMessageDialog( null, " Erro - Não foi possível ler o arquivo '" + arquivo +"' ", "Erro ao ler arquivo", JOptionPane.ERROR_MESSAGE );
        }
        return cont;
    }

    /**
     * retorna todo o conteudo do 'arquivo'
     * @param arquivo
     * @return
     */
    public static String leArquivo( String arquivo ) {
        try {
            BufferedReader in = new BufferedReader( new FileReader( arquivo ) );
            String line;
            StringBuffer buffer = new StringBuffer();
            while( (line = in.readLine()) != null )
                buffer.append( line + "\n");
            in.close();
            conteudoDoArquivo = buffer.toString();
        }
        catch( IOException ioe ) {
            ioe.printStackTrace();
            JOptionPane.showMessageDialog( null, " Erro - Não foi possível ler o arquivo '" + arquivo +"' ", "Erro ao ler arquivo", JOptionPane.ERROR_MESSAGE );
        }
        catch( Exception e ) {
            e.printStackTrace();
            JOptionPane.showMessageDialog( null, " Erro - Não foi possível ler o arquivo '" + arquivo +"' ", "Erro ao ler arquivo", JOptionPane.ERROR_MESSAGE );
        }
        return conteudoDoArquivo;
    }

    /**
     * apaga todo o conteudo do 'arquivo'
     * @param arquivo
     */
    public static void limpaArquivo( String arquivo )
    {
        File file = new File( arquivo );
        if( file.exists() ) {
            FileOutputStream output;
            try {
                output = new FileOutputStream( file );
                output.flush();
                output.close();
            }
            catch( IOException ioe ) {
                ioe.printStackTrace();
                JOptionPane.showMessageDialog( null, " Erro - Não foi possível apagar o conteúdo do arquivo ", "Erro", JOptionPane.ERROR_MESSAGE );
            }
            catch( Exception e ) {
                e.printStackTrace();
                JOptionPane.showMessageDialog( null, " Erro - Não foi possível apagar o conteúdo do arquivo ", "Erro", JOptionPane.ERROR_MESSAGE );
            }
        }
    }

    /**
     * deleta o 'arquivo'
     * @param arquivo
     */
    public static void deletaArquivo( String arquivo ) {
        boolean del = true;
        File file = new File( arquivo );
        if( file.exists() ) {
            del = file.delete();
        }
        if( !del ) {
            //JOptionPane.showMessageDialog( null, " Erro - Não foi possível deletar o arquivo ", "Erro", JOptionPane.ERROR_MESSAGE );
        }
    }

    /**
     * retorna a linha do 'arquivo' indicada por 'linha'
     * @param arquivo
     * @param linha
     * @return
     */
    public static String leLinha( String arquivo, int linha ) {
        String aux = "";
        int cont = 0;
        try{
            RandomAccessFile random = new RandomAccessFile( arquivo, "r" );
            while( cont <= linha ){
                aux = random.readLine();

                cont++;
            }
        }
        catch( IOException ioe ) {
            ioe.printStackTrace();
            JOptionPane.showMessageDialog( null, " Erro - Não foi possível ler o arquivo ", "Erro", JOptionPane.ERROR_MESSAGE );
        }
        catch( Exception e ) {
            e.printStackTrace();
            JOptionPane.showMessageDialog( null, " Erro - Não foi possível ler o arquivo ", "Erro", JOptionPane.ERROR_MESSAGE );
        }
        return aux;
    }

    /**
     *
     * @param arquivo
     * @return
     */
    public static boolean existeArquivo(String arquivo) {
        File file = new File( arquivo );
        if( file.exists() )
            return true;
        return false;
    }

    /**
     *
     * @param arquivo
     * @param diretorio
     */
    public static void moveArquivo(String arquivo, String diretorio){
        Process mover;
        try {
            mover = Runtime.getRuntime().exec(("mv ")+(arquivo)+(" ")+(diretorio));
            mover.waitFor();
          } catch (IOException ex) {
            System.out.println("Erro na execução do comando:\n" + ex.getMessage());
        } catch (InterruptedException ex) {
                System.out.println("Erro na execução do comando:\n" + ex.getMessage());
        }
        //return "";
    }

}
