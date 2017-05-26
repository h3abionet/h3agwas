/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package h3agwas.config;


import java.io.FileInputStream;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import org.apache.poi.xssf.usermodel.XSSFSheet;
import org.apache.poi.xssf.usermodel.XSSFWorkbook;
import org.apache.poi.ss.usermodel.Cell;
import org.apache.poi.ss.usermodel.Row;

/**
 *
 * @author scott
 */
public class H3agwasConfig {

    
    static HashMap<String,String> harg = new HashMap<>();
    
    private static String showArg(String arg) {
        String out[] = {arg,harg.get(arg)};
        return
        String.format("    %-20s = \"%s\"\n",out);
    }
    
    private static String common_docker_options = 
                  "        process.$removeDuplicateSNPs.container = \"$plinkImage\"\n"
                + "        process.$identifyIndivDiscSexinfo.container = \"$plinkImage\"\n"
                + "        process.$calculateSampleMissing.container = \"$plinkImage\"\n"
                + "        process.$calculateSampleHetrozygosity.container = \"$plinkImage\"\n"
                + "        process.$compPCA.container = \"$plinkImage\"\n"
                + "        process.$pruneForIBD.container = \"$plinkImage\"\n"
                + "        process.$removeQCIndivs.container = \"$plinkImage\"\n"
                + "        process.$calculateMaf.container = \"$plinkImage\"\n"
                + "        process.$calculateSnpMissingness.container = \"$plinkImage\"\n"
                + "        process.$calculateSnpSkewStatus.container = \"$plinkImage\"\n"
                + "        process.$removeQCPhase1.container = \"$plinkImage\"\n"
                + "        process.$computePhase0.container = \"$plinkImage\"\n"
                + "        process.$drawPCA.container = \"$rEngineImage\"\n"
                + "        process.$generateIndivMissingnessPlot.container=\"$rEngineImage\"\n"
                + "        process.$generateMissHetPlot.container = \"$rEngineImage\"\n"
                + "        process.$generateMafPlot.container = \"$rEngineImage\"\n"
                + "        process.$generateSnpMissingnessPlot.container = \"$rEngineImage\"\n"
                + "        process.$generateDifferentialMissingnessPlot.container = \"$rEngineImage\"\n"
                + "        process.$generateHwePlot.container = \"$rEngineImage\"\n"
                + "        process.$produceReports.container =\"$latexImage\"\n"
                + "        process.$lgen2ped.container = \"plinkImage\"\n"
                + "        process.$create_bed = \"plinkImage\"\n"
                + "        process.$createPCA = \"plinkImage\"\n"
                + "        process.$computeTest = \"gemmaImage\"\n"
    + "         + "        
                + "\n"
                + "        docker.remove = true\n"
                + "        docker.runOptions = '--rm'\n"
                + "\t      docker.registry = 'quay.io'\n"
                + "        docker.enabled = true\n"
                + "        docker.temp = 'auto'\n"
                + "        docker.fixOwnership = true\n";
    
    private static String template() {
        String template = 
                  "plinkImage = '"+harg.get("plinkImage")+"'\n"
                + "rEngineImage = '"+harg.get("rEngineImage")+"'\n"
                + "latexImage = '"+harg.get("latexImage")+"'\n"
                + "swarmPort = '2376'\n"
                + "\n"
                + "manifest {\n"
                + "    homePage = 'http://github.com/h3abionet/h3agwas'\n"
                + "    description = 'GWAS Pipeline for H3Africa'\n"
                + "    mainScript = 'plink-qc.nf'\n"
                + "}\n"
                + "\n\n"
                + "aws {\n"
                + "    accessKey ='"+harg.get("accessKey")+"'\n"
                + "    secretKey ='"+harg.get("secretKey")+"'\n"
                + "    region    ='"+harg.get("region")+"'\n"
                + "}\n\n"
                + "    cloud {\n" +
                "\n" +
                "            imageId = \""+harg.get("AMI")+"\"      // specify your AMI id here\n" +
                "            instanceType = \""+harg.get("instanceType")+"\"\n" +
                "            subnetId = \""+harg.get("subnetid")+"\"\n" +
                "            sharedStorageId   = \""+harg.get("sharedStorageId")+"\"\n"+
                "   	     sharedStorageMount = \""+harg.get("sharedStorageMount")+"\"\n"+
                "            bootStorageSize = \""+harg.get("bootStorageSize")+"\"     // Size of disk for images spawned\n" +
                "//          instanceStorageMount = \"\"   // Set a common mount point for images\n" +
                "//          instanceStorageDevice = \"\"  // Set a common block device for images\n" +
                "            autoscale {\n" +
                "               enabled = true\n" +
                "               maxInstances = "+harg.get("maxInstances")+"\n" +
                "               terminateWhenIdle = true\n" +
                "             }\n" +
                "\n" +
                "    }\n\n\n"
                + "params {\n"
                + "\n"
                + "    // Directories\n"
                + "    work_dir                = \""+harg.get("work_dir")+"\"\n"
                + "    input_dir               = \""+harg.get("input_dir")+"\"\n"
                + "    output_dir              = \""+harg.get("output_dir")+"\"\n"
                + "    scripts                 = \""+harg.get("scripts")+"\"\n"
                + "\n"
                + "    // Data\n";
        String parms [] = {"input_pat","high_ld_regions_fname","sexinfo_available",
                           "cut_het_high","cut_het_low","cut_miss","cut_diff_miss",
                           "cut_maf","cut_mind","cut_geno","cut_hwe","pi_hat",
                           "plink_process_memory","other_process_memory",
                           "max_plink_cores","accessKey","secretKey","region","AMI","instanceType","bootStorageSize","maxInstances"};
        for (String arg : parms)
           template = template + showArg(arg);
        template = template
                + "\n"
                + "}\n"
                + "profiles {\n"
                + "\n"
                + "    // For execution on a local machine, no containerization. -- Default\n"
                + "    standard {\n"
                + "        process.executor = 'local'\n"
                + "    }\n"
                + "\n"
                + "    // For execution on a PBS scheduler, no containerization.\n"
                + "    pbs {\n"
                + "        process.executor = 'pbs'\n"
                + "        process.queue = '"+harg.get("queue")+"'\n"
                + "    }\n"
                + "\n"
                + "    // For execution on a PBS scheduler with containerization.\n"
                + "    pbsDocker {\n"
                + "\n"
                + "        process.executor = 'pbs'\n"
                + common_docker_options
                + "\n"
                + "    }\n"
                + "\n"
                + "    // Execute pipeline with Docker locally\n"
                + "    docker {\n"
                + common_docker_options
                + "        docker.process.executor = 'local'\n"
                + "    }\n"
                + "\n"
                + "    dockerpbs {\n"
                + "        process.executor = 'pbs'\n"
                + common_docker_options
                + "        docker.process.executor = 'local'\n"
                + "        docker.fixOwnership = true\n"
                + "    }\n"
                + "\n"
                + "\n"
                + "    // Execute pipeline with Docker Swarm setup\n"
                + "    dockerSwarm {\n"
                + "\n"
                + common_docker_options
                + "        docker.process.executor = 'local'\n"
                + "        docker.engineOptions = \"-H :$swarmPort\"\n"
                + "    }\n"
                + "\n"
                + "\n"
                + "\n"
                + "}\n"
                + "";       
       return template;
        
    }
    
    
    private static void getOptions(Path conffile) throws IOException {
        XSSFSheet sheet;
        Row row;
        Iterator<Row> rows;
        String var_name,def_val, use_val;
        FileInputStream f = new FileInputStream(conffile.toString());
        rows = new XSSFWorkbook(f).getSheetAt(0).rowIterator();
        while (rows.hasNext()) {
            row = rows.next();
            if (row.getCell(0).getStringCellValue().charAt(0)=='#') continue;
            var_name = row.getCell(1).getStringCellValue();
            def_val  = row.getCell(2) != null ? row.getCell(2).getStringCellValue() :  "";
            use_val  = (row.getCell(4) != null) ? row.getCell(4).getStringCellValue() : "";
            if (use_val.length() == 0 ) use_val = def_val;
            harg.put(var_name, use_val);
        }     
    }
    /**
     * @param args the command line arguments
     * @throws java.io.IOException
     */
    public static void main(String[] args) throws IOException {
        getOptions(Paths.get(args[0]));
        System.out.println(template());
    }
    
}
