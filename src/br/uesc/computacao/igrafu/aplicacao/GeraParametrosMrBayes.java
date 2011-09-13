package br.uesc.computacao.igrafu.aplicacao;

import java.util.ArrayList;

/**
 * Classe que gera os parametros de forma separada por objetivo:
 * Configurar modelo, mcmc, arvore, sumarização de arvores e/ou de parametros.
 * Divide também parametros dos tipos de dados Dna ou Proteina.
 * Essa divisão corresponde a divisão grafica presente em MrBayesUI
 *
 * Os diversos metodos formatam os dados recebidos para o formato aceitavel
 * pelo MrBayes_3.2.1 colocando os parametros com seus determinados comandos.
 * Para criar a forma do arquivo batch use o metodo batchContent().
 * Ele utiliza a formatação gerada pelos outros metodos e alinha os parametros 
 * por comandos gerais como o lset e o mcmc.
 *
 * Portando para usar o batchContent() você deverá ter chamado os outros metodos.
 *
 * @author Anderson Carlos Sousa e Santos </br>
 * Orientadora: Martha Ximena Torres Delgado
 *
 * @see br.uesc.computacao.igrafu.apresentacao.MrBayesUI
 */
public class GeraParametrosMrBayes {

    private String mcmcParam = "mcmc";

    private String sumpParam = "";

    private String sumtParam = "";

    private String lsetParam = "lset";

    private String prsetParam = "prset";

    private String userTreeParam = "";

    private final String setParam = "set autoclose=yes nowarnings=yes";



    /**
     * Formata os parametros comuns ao modelo de Dna e Proteina.
     *
     * @param covarion - rececebe yes/no
     * @param rates - representa parametro de mesmo nome
     * @param ngammacat - representa parametro de mesmo nome
     * @param ratecorrpr - representa parametro de mesmo nome
     * @param pinvarpr - representa parametro de mesmo nome
     */
    public void setParamModelo(boolean covarion,
                                String rates,
                                int ngammacat,
                                ArrayList<Integer> ratecorrpr,
                                ArrayList<Integer> pinvarpr) {

            rates = rates.toLowerCase();
            String covar = " covarion=";

            if(covarion)
                    covar += "yes";
            else
                    covar +="no";

            lsetParam += covar;

            lsetParam += " rates="+rates;

            if(!rates.equals("propinv") && !rates.equals("equal"))
                    lsetParam += " ngammacat="+ngammacat;

            if(rates.equals("adgamma")){
                    if(ratecorrpr.size() == 1)
                            prsetParam += " ratecorrpr=fixed(" + ratecorrpr.get(0) + ")";
                    else
                            prsetParam += " ratecorrpr=uniform(" + ratecorrpr.get(0) + "," + ratecorrpr.get(1) + ")" ;
            }
            if(rates.equals("propinv")){
                    if(pinvarpr.size() == 1)
                            prsetParam += " pinvarpr=fixed(" + pinvarpr.get(0) + ")";
                    else
                            prsetParam += " pivarrpr=uniform(" + pinvarpr.get(0) + "," + pinvarpr.get(1) + ")";
            }
    }



    /**
     * Formata os parametros especificos para DNA.
     *
     * @param nucmodel
     * @param substitutionModelIndex
     * @param codeIndex
     * @param omegavarIndex
     */
    public void setParamDna(String nucmodel,
                            int substitutionModelIndex,
                            int codeIndex, int omegavarIndex) {

        nucmodel = nucmodel.toLowerCase();
        
        lsetParam += " nucmodel="+nucmodel;

        switch(substitutionModelIndex){
            case 0://GTR
                lsetParam += " nst=6";
                break;
            case 1://SYM
                lsetParam += " nst=6";
                prsetParam += " statefreqpr=fixed(0.25, 0.25, 0.25, 0.25)";
                break;
            case 2://HKY
                lsetParam += " nst=2";
                break;
            case 3://K2P
                lsetParam += " nst=2";
                prsetParam += " statefreqpr=fixed(equal)";
                break;
            case 4://F81
                lsetParam += " nst=1";
                break;
            case 5://JC
                lsetParam += " nst=1";
                prsetParam += " statefreqpr=fixed(equal)";
                break;
            default:
                System.out.println("Erro no index de Substitution Model");
                break;
        }

        if(nucmodel.equals("codon")){
            switch(codeIndex){
                case 0://universal
                    lsetParam += " code=universal";
                    break;
                case 1://vertmt
                    lsetParam += " code=vertmt";
                    break;
                case 2://mycoplasma
                    lsetParam += " code=mycoplasma";
                    break;
                case 3://yeast
                    lsetParam += " code=yeast";
                    break;
                case 4://ciliates
                    lsetParam += " code=ciliates";
                    break;
                case 5://metmt
                    lsetParam += " code=metmt";
                    break;
                default:
                    System.out.println("ERRO no index de Genetic Code");
                    break;
            }

            switch(omegavarIndex){
                case 0://equal
                    lsetParam += " omegavar=equal";
                    break;
                case 1://ny98
                    lsetParam += " omegavar=ny98";
                    break;
                case 3://m3
                    lsetParam += " omegavar=m3";
                    break;
                default:
                    System.out.println("ERRO no index de Omega Variation");
                    break;
            }
        }
            

    }


    /**
     * Formata os parametros especificos para Proteina.
     *
     * @param matrixRate
     * @param fixedModelIndex
     * @param varModelIndex
     */
    public void setParamProteina(String matrixRate,
                                int fixedModelIndex,
                                int varModelIndex) {

        lsetParam = "";
        matrixRate = matrixRate.toLowerCase();
        if(matrixRate.equals("fixed")){
            switch(fixedModelIndex){
                case 0://poisson
                    prsetParam += " aamodelpr=fixed(poisson)";
                    break;
                case 1://jones
                   prsetParam += " aamodelpr=fixed(jones)";
                   break;
                case 2://dyhoff
                   prsetParam += " aamodelpr=fixed(dyhoff)";
                   break;
                case 3://mtrev
                   prsetParam += " aamodelpr=fixed(mtrev)";
                   break;
                case 4://mtmam
                   prsetParam += " aamodelpr=fixed(mtmam)";
                   break;
                case 5://wag
                   prsetParam += " aamodelpr=fixed(wag)";
                   break;
                case 6://rtrev
                   prsetParam += " aamodelpr=fixed(rtrev)";
                   break;
                case 7://cprev
                   prsetParam += " aamodelpr=fixed(cprev)";
                   break;
                case 8://vt
                   prsetParam += " aamodelpr=fixed(vt)";
                   break;
                case 9://blosum
                   prsetParam += " aamodelpr=fixed(blosum)";
                   break;
                default:
                   System.out.println("ERRO no index de Fixed Rate Models");
                   break;

            }
        }//End if

        else if(matrixRate.equals("variable"))
        {
            switch(varModelIndex){
                case 0://equalin
                   prsetParam += " aamodelpr=fixed(equalin)";
                   break;
                case 1://gtr
                    prsetParam += " aamodelpr=fixed(gtr)";
                   break;
                default:
                   System.out.println("ERRO no index de Variable Rate Models");
                   break;

            }
        }else
            prsetParam += " mixed";


    }

    /**
     * Formata os parametros mais comuns do MCMC.
     *
     * @param ngen
     * @param nrun
     * @param samplefreq
     * @param nchain
     * @param swapfreq
     * @param nswap
     * @param diagnfreq
     * @param minpartfreq
     * @param stoprule
     * @param stopval
     */
    public void setParamMcmc(int ngen,
                             int nrun,
                             int samplefreq,
                             int nchain,
                             int swapfreq,
                             int nswap,
                             int diagnfreq,
                             double minpartfreq,
                             boolean stoprule,
                             double stopval) {

        mcmcParam += " ngen="+ngen;
        mcmcParam += " nrun="+nrun;
        mcmcParam += " samplefreq="+samplefreq;
        mcmcParam += " nchain="+nchain;
        if(nchain > 1){
            mcmcParam += " swapfreq="+swapfreq;
            mcmcParam += " nswap="+nswap;
        }
        mcmcParam += " diagnfreq="+diagnfreq;
        mcmcParam += " minpartfreq="+minpartfreq;
        if(nrun > 1){
            if(stoprule){
                mcmcParam += " stoprule=yes";
                mcmcParam += " stopval="+stopval;
            }
            else
                mcmcParam += " stoprule=no";
        }

        //Parametros fixos
        mcmcParam += " mcmcdiagn=yes";
    }


    /**
     * Formata os parametros mais avançados do MCMC.
     *
     * @param seed
     * @param swapseed
     * @param savebrlens
     * @param ordertaxa
     * @param reweightDecrease
     * @param reweightIncrease
     * @param reweightIncrement
     * @param relburnin
     * @param burninfraq
     * @param burnin
     */
    public void setParamAdvancedMcmc(long seed,
                                    long swapseed,
                                    boolean savebrlens,
                                    boolean ordertaxa,
                                    double reweightDecrease,
                                    double reweightIncrease,
                                    double reweightIncrement,
                                    String relburnin,
                                    double burninfraq,
                                    int burnin) {

        relburnin = relburnin.toLowerCase();

        mcmcParam += " seed="+seed;
        mcmcParam += " swapseed="+swapseed;
        if(savebrlens)
            mcmcParam += " savebrlens=yes";
        else
            mcmcParam += " savebrlens=no";

        if(ordertaxa)
            mcmcParam += " ordertaxa=yes";
        else
            mcmcParam += " ordertaxa=no";

        mcmcParam += " reweight=("+reweightDecrease+","
                                  +reweightIncrease+","
                                  +reweightIncrement+")";
        if(relburnin.equals("proportional")){
            mcmcParam += " relburnin=yes";
            mcmcParam += " burninfrac="+burninfraq;
        }
        else{
            mcmcParam += " relburnin=no";
            mcmcParam += " burnin="+burnin;
        }
    
    }


    /**
     * Formata os parametros relativos a arvore.
     *
     * @param startingTree
     * @param userTree
     * @param nperts
     * @param brlensprClock
     * @param clockModelIndex
     * @param brlensNum
     * @param treeHeightpr
     * @param thetapr
     * @param ploidy
     * @param speciationpr
     * @param sampleprob
     * @param extinctionpr
     */
    public void setParamArvore(String startingTree,
                               String userTree,
                               int nperts,
                               String brlensprClock,
                               int clockModelIndex,
                               int brlensNum,
                               int treeHeightpr,
                               int thetapr,
                               String ploidy,
                               int speciationpr,
                               int sampleprob,
                               int extinctionpr) {

        startingTree = startingTree.toLowerCase();
        userTree = userTree.toLowerCase();
        brlensprClock = brlensprClock.toLowerCase();
        ploidy = ploidy.toLowerCase();

        if(startingTree.equals("user")){
            userTreeParam = "UserTree=";
            userTreeParam += userTree;
            mcmcParam += " nperts="+nperts;
        }
        mcmcParam += " startingtree=" + startingTree;

        if(brlensprClock.equals("non-clock")){
            prsetParam += " brlenspr=unconstrained:exponential("+brlensNum+")";
        }
        else
        switch(clockModelIndex){
        case 0://simple clock
            prsetParam += " brlenspr=clock:uniform";
            prsetParam += " treeheightpr=exponential("+treeHeightpr+")";
            break;
        case 1://coalescence
            prsetParam += " brlenspr=clock:coalescence";
            prsetParam += " thetapr=exponential("+thetapr+")";
            lsetParam += " ploidy="+ploidy;
        break;
        case 2://birth death
            prsetParam += " brlenspr=clock:birthdeath";
            prsetParam += " speciationpr=exponential("+speciationpr+")";
            prsetParam += " extinctionpr=exponential("+extinctionpr+")";
            prsetParam += " sampleprob="+sampleprob;
            break;
        default:
            System.out.println("ERRO no index de Clock");
        break;

        }



                                }

    /**
     * Formata os parametros relativos ao comando sump
     *
     * @param burnin
     * @param plot
     * @param marglike
     * @param table
     */
    public void setParamSump(int burnin,
                             boolean plot,
                             boolean marglike,
                             boolean table
                             ) {
        sumpParam = "sump";

        sumpParam += " burnin="+burnin;
        if(plot)
            sumpParam += " plot=yes";
        else
            sumpParam += " plot=no";

        if(marglike)
            sumpParam += " marglike=yes";
        else
            sumpParam += " marglike=no";

        if(table)
            sumpParam += " table=yes";
        else
            sumpParam += " table=no";

        //Parametros fixos
        sumpParam += " printtofile=yes";

    }

    /**
     * Formata os parametros relativos ao comando sumt
     *
     * @param burnin
     * @param displaygeq
     * @param contype
     * @param calctreeprobs
     */
    public void setParamSumt(int burnin,
                             double displaygeq,
                             String contype,
                             boolean calctreeprobs) {

        sumtParam = "sumt ";
        
        sumtParam += " burnin=" + burnin;

        sumtParam += " displaygeq="+ displaygeq;

        sumtParam += " contype=" + contype;

        if(calctreeprobs)
            sumtParam += " calctrprobs=yes";
        else
            sumtParam += " calctrprobs=no";

        //Parametros fixos
        sumtParam += " showtreeprobs=no";
    }


    /**
     * Organiza e formata a os comandos seguidos de subcomando os parametros 
     * assumindo que estes foram formatados pelos metodos anteriores. 
     * Caso não tenham sido o retorno será um arquivo com os comandos gerais 
     * sem parametros ou subcomandos.
     * 
     * Este metodo adiciona também o comando para executar o arquivo de entrada.
     *
     * @return Conteudo de um arquivo batch para execução do mrbayes
     */
    public String batchContent(){
        String batch = "begin mrbayes;\n\t";

            batch += this.getSetParam() + ";\n";
            batch += "\t" + "execute " + InputSequence.getFile().getPath() + ";\n";
            batch += "\t mcmcp printfreq=1;\n";

            if(!lsetParam.equals(""))
                batch += "\t" + this.getLsetParam() + ";\n";

            batch += "\t" + this.getPrsetParam() + ";\n";

            if(!getUserTreeParam().equals(""))
                batch += "\t" + this.getUserTreeParam() + ";\n";

            batch += "\t" + this.getMcmcParam() + ";\n";

            if(!getSumpParam().isEmpty())
                batch += "\t" + this.getSumpParam() + ";\n";

            if(!getSumtParam().isEmpty())
                batch += "\t" + this.getSumtParam() + ";\n";

       batch += "end;";


        return batch;
    }



    /**
     * @return the mcmcParam
     */
    public String getMcmcParam() {
        return mcmcParam;
    }

    /**
     * @return the sumpParam
     */
    public String getSumpParam() {
        return sumpParam;
    }

    /**
     * @return the sumtParam
     */
    public String getSumtParam() {
        return sumtParam;
    }

    /**
     * @return the lsetParam
     */
    public String getLsetParam() {
        return lsetParam;
    }

    /**
     * @return the prsetParam
     */
    public String getPrsetParam() {
        return prsetParam;
    }

    /**
     * @return the userTreeParam
     */
    public String getUserTreeParam() {
        return userTreeParam;
    }

    /**
     * @return the setParam
     */
    public String getSetParam() {
        return setParam;
    }


}