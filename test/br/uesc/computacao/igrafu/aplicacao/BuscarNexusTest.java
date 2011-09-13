package br.uesc.computacao.igrafu.aplicacao;

import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author Anderson Carlos Sousa e Santos </br>
 * Orientadora: Martha Ximena Torres Delgado
 */
public class BuscarNexusTest {


    /**
     * 
     * @throws Exception
     */
    @BeforeClass
    public static void setUpClass() throws Exception {
    }

    /**
     *
     * @throws Exception
     */
    @AfterClass
    public static void tearDownClass() throws Exception {
    }

    /**
     * Test of defineDataType method, of class BuscarNexus.
     */
    @Test
    public void testDefineDataType() {
        System.out.println("defineDataType");
        String file = "/home/acarlos/mrbayes-3.1.2/cynmix.nex";
        String expResult = "tipo_fail";
        String result = BuscarNexus.defineDataType(file);
        assertEquals(expResult, result);

        file = "/home/acarlos/mrbayes-3.1.2/gpl.txt";
        expResult = "no_datatype";
        result = BuscarNexus.defineDataType(file);
        assertEquals(expResult, result);

        file = "/home/acarlos/mrbayes-3.1.2/primates.nex";
        expResult = "dna";
        result = BuscarNexus.defineDataType(file);
        assertEquals(expResult, result);

        file = "/home/acarlos/mrbayes-3.1.2/avian_ovomucoids.nex";
        expResult = "protein";
        result = BuscarNexus.defineDataType(file);
        assertEquals(expResult, result);
        
    }

}