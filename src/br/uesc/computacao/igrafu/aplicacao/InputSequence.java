package br.uesc.computacao.igrafu.aplicacao;


import br.uesc.computacao.util.ManipulaArquivo;
import java.io.File;

/**
 * Classe que abstrai o arquivo de entrada.
 * Possui informações uteis sobre o arquivo de entrada, para serem usadas por 
 * outras rotinas que operam sobre o arquivo.
 *
 * @author Anderson Carlos Sousa e Santos </br>
 * Orientadora: Martha Ximena Torres Delgado
 */
public class InputSequence {

    private static String dataType;

    private static File file;

    private static String format;

    /**
     * Construtor que recebe inicializa todas as variaveis.
     *
     * @param dataType - tipo de dado utilzado na sequencia (Ex.: DNA)
     * @param filePath - caminho absoluto do arquivo de sequencia
     * @param format - formato utilizado no arquivo. (Ex.: Nexus)
     */
    public InputSequence(String dataType, String filePath, String format) {
        setDataType(dataType);
        setFile(filePath);
        setFormat(format);
    }

    /**
     * Cria arquivo (objeto da Classe File) através do caminho(path) do arquivo
     * e seta variavel file com esse novo objeto
     *
     * @param filePath - caminho absoluto do arquivo de sequencia
     * @return <ul><li> true - se o arquivo recebido existe
     * <li> false - se o arquivo recebido não existe.
     */
    public static boolean setFile(String filePath) {
        if(ManipulaArquivo.existeArquivo(filePath)){
            file = new File(filePath);
            return true;
        }
        else
            return false;

    }

    /**
     * Seta vaviavel file com o objeto newFile
     *
     * @param newFile - Arquivo de sequencia
     * @return <ul><li> true - se o arquivo recebido não é nulo
     * <li> false - se o arquivo recebido é nulo.
     */
    public boolean setFile(File newFile){
        if(newFile != null){
            file = newFile;
            return true;
        }
        else
            return false;
    }

    /**
     * Atribui o tipo de dados que está sendo descrito no arquivo.
     * Guarda o tipo de dados sempre em minusculo.
     *
     * @param newDataType - tipo de dado utilzado na sequencia (Ex.: DNA)
     */
    public static void setDataType(String newDataType) {
        dataType = newDataType.toLowerCase();
    }

    /**
     *
     * @return tipo de dado utilzado na sequencia (Ex.: DNA)
     */
    public static String getDataType() {
            return dataType;
    }

    /**
     *
     * @return arquivo contendo sequencias de dna ou proteina
     */
    public static File getFile() {
            return file;
    }

    /**
     * Metodo para definir o formato do arquivo de sequencias.
     * Guarda o formato sempre em minusculo.
     *
     * @param format configura o formato de arquivos de sequencia (Ex.: Nexus)
     */
    public static void setFormat(String format) {
            format = format.toLowerCase();
    }

    /**
     *
     * @return formato utilizado no arquivo.
     */
    public static String getFormat() {
            return format;
    }
}
