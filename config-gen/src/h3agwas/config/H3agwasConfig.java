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
import java.util.LinkedHashMap;
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

    
    static LinkedHashMap<String,String> harg = new LinkedHashMap<>();
    
    private static String showArg(String arg) {
        String out[] = {arg,harg.get(arg)};
        return
        String.format("    %-20s = \"%s\"\n",out);
    }
    

    private static String template() {
        String template = "";
        for (String key : harg.keySet()) {
            template=  template+key+"=\""+harg.get(key)+"\"\n";
        }     
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
            if (row.getCell(0).getStringCellValue().charAt(0)=='#') {
                 System.out.println();
                 continue;
            }
            var_name = row.getCell(1).getStringCellValue();
            def_val  = row.getCell(2) != null ? row.getCell(2).getStringCellValue() :  "";
            use_val  = (row.getCell(4) != null) ? row.getCell(4).getStringCellValue() : "";
            if (use_val.length() == 0 ) use_val = def_val;
            if (!(use_val.matches("^[+-]?\\d+\\.?\\d*$"))) 
                use_val = String.format("\"%s\"", use_val);
            //harg.put(var_name, use_val);
            System.out.println(String.format("%-20s = %s",var_name, use_val));
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
