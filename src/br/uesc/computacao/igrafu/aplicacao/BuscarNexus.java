package br.uesc.computacao.igrafu.aplicacao;

import br.uesc.computacao.util.ManipulaArquivo;
/**
 * Classe que realiza Busca em arquivo Nexus. 
 *
 * @author Anderson Carlos Sousa e Santos </br>
 * Orientadora: Martha Ximena Torres Delgado
 */

 public class BuscarNexus {


    /**
     *
     * @param file arquivo com formato nexus
     * @return <ul> 
     * 
     * <li> "dna" -  caso o parametro datatype do nexus seja dna
     * <li> "protein" -  caso o parametro datatype do nexus seja protein
     * <li> "tipo_fail" - caso datatype exista no arquivo, mas não está
     * definido como dna ou protein;
     * <li> "no_datatype" - caso em que não foi encontrado o parametro datatype.</br>
     * Este retorno pode ser usado para definir que o formato não é nexus.
     *
     */
    public static String defineDataType(String file){
        String datatype = "notEmpty";

            datatype = ManipulaArquivo.leArquivo(file);
            datatype.toLowerCase();

            if(datatype.indexOf("datatype=") != -1 ){
                if ( datatype.indexOf("datatype=dna") != -1 ){
                   return "dna";
                }
                else if( datatype.indexOf("datatype=protein") != -1 ){
                   return "protein";
                }
                else
                    return "tipo_fail";
            }

            return "no_datatype";
        }

}
