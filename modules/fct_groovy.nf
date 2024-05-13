
def strmem(val){
 return val as nextflow.util.MemoryUnit
}


// Checks if the file exists
def  fileexist(filec) { 
   File fn =  new File(filec)
   if (fn.exists())
       return fn;
    else
       error("\n\n-----------------\nFile $fn does not exist\n\n---\n")
}

def checkresume(){
 if (!workflow.resume) {
    def dir = new File(params.output_dir)
    if (dir.exists() && dir.directory && (!(dir.list() as List).empty)) {
       println "\n\n============================================"
       println "Unless you are doing a -resume, the output directory should be empty"
       println "We do not want to overwrite something valuable in "+params.output_dir
       println "Either clean your output directory or check if you meant to do a -resume"
       System.exit(-1)
    }
 }
}

def is_nullfile(filename){
	nullfile = [false,"False","false", "FALSE",0,"","0","null",null]
	return nullfile.contains(filename)
}

//This method first checks that the data file has the stated column
// If so, it creates a channel for it
// NB: if the file is in S3 we cannot do the test since Groovy does not
// allow us to access the file directly
def fileheader_create_ch = { parm, parm_name, col_name ->
  if (parm.toString().contains("s3://")) {
    println "The file <$parm> is in S3 so we cannot do a pre-check";
    return Channel.fromPath(parm, checkIfExists:true);
  }
  if (parm.toString().contains("az://")) {
    println "The file <$parm> is in Azure so we cannot do a pre-check";
    return Channel.fromPath(parm, checkIfExists:true);
  }
  if (is_nullfile(parm)) {
    filename = "emptyZ0${parm_name}.txt";
    new File(filename).createNewFile()
    new_ch = Channel.fromPath(filename,checkIfExists:true);

  } else {
    if (! file(parm).exists()) {
     error("\n\nThe file <$parm> given for <params.${parm_name}> does not exist")
    } else {
      def line
      new File(parm).withReader { line = it.readLine() }
      fields = line.split()
      if (! fields.contains(col_name))
          error("\n\nThe file <$parm> given for <params.${parm_name}> does not have a column <${col_name}>\n")
    }
    new_ch = Channel.fromPath(parm,checkIfExists:true);
  }
  return new_ch;
}


def checkparams(param, namesparam, type, min=null, max=null, possibleval=null, notpossibleval=null) {
  messageerror=""
  if(param==null){
    messageerror+="error :--"+namesparam+" is null "
  } else {
    if(!(param.getClass() in type)){
   messageerror+="error :--"+namesparam+" must be a "+ type
     if(params.getClass()==Boolean)messageerror+=", but no parameters given"
     else messageerror+=" but type is "+param.getClass()+" value "+ param
   }else{
   if(min && param<min)messageerror+="\nerror : --"+namesparam+" < min value :"+param +" < "+min
   if(max && param>max)messageerror+="\nerror : --"+namesparam +"> maxvalue :" + param+" > "+max
   if(possibleval && !(param in possibleval))messageerror+="\nerro : --"+namesparam +" must be one the value :"+possibleval.join(',')
   }
   }
    errormess(messageerror,2)
}


//This method first checks that the data file has the stated column
// If so, it creates a channel for it
// NB: if the file is in S3 we cannot do the test since Groovy does not
// allow us to access the file directly
def fileheader_create_ch ( parm, parm_name, col_name ){
  if (parm.toString().contains("s3://")) {
    println "The file <$parm> is in S3 so we cannot do a pre-check";
    return Channel.fromPath(parm, checkIfExists:true);
  }
  if (parm.toString().contains("az://")) {
    println "The file <$parm> is in Azure so we cannot do a pre-check";
    return Channel.fromPath(parm, checkIfExists:true);
  }
  if ((parm==0) || (parm=="0") || (parm==false) || (parm=="false")) {
    filename = "emptyZ0${parm_name}.txt";
    new File(filename).createNewFile()
    new_ch = Channel.fromPath(filename,checkIfExists:true);

  } else {
    if (! file(parm).exists()) {
     error("\n\nThe file <$parm> given for <params.${parm_name}> does not exist")
    } else {
      def line
      new File(parm).withReader { line = it.readLine() }
      fields = line.split()
      if (! fields.contains(col_name))
          error("\n\nThe file <$parm> given for <params.${parm_name}> does not have a column <${col_name}>\n")
    }
    new_ch = Channel.fromPath(parm,checkIfExists:true);
  }
  return new_ch;
}



def checkmultiparam(params, listparams, type, min=null, max=null, possibleval=null, notpossibleval=null){
 messageerror=""
 for(param in listparams){
   if(params.containsKey(param)){
     checkparams(params[param], param, type, min=min, max=max, possibleval=possibleval, notpossibleval=notpossibleval)
   }else{
     messageerror+="param :"+param+" not initialize\n"
   }
 }
 errormess(messageerror, 2)
}

def checkColumnHeader(fname, columns) {
  if (workflow.profile == "awsbatch") return;
  if (fname.toString().contains("s3://")) return;
  if (fname.contains("az://") ) return;
  if (is_nullfile(fname)) return;
  print(fname)
  new File(fname).withReader { line = it.readLine().tokenize() }
  problem = false;
  columns.each { col ->
    if (! line.contains(col) ) {
      println "The file <$fname> does not contain the column <$col>";
      problem=true;
    }
    if (problem)
      System.exit(2)
  }
}




def getConfig = {
  all_files = workflow.configFiles.unique()
  text = ""
  all_files.each { fname ->
      base = fname.baseName
      curr = "\n\n*-subsection{*-protect*-url{$base}}@.@@.@*-footnotesize@.@*-begin{verbatim}"
      file(fname).eachLine { String line ->
        if (line.contains("secretKey")) { line = "secretKey='*******'" }
        if (line.contains("accessKey")) { line = "accessKey='*******'" }
        curr = curr + "@.@"+line
      }
      curr = curr +"@.@*-end{verbatim}\n"
      text = text+curr
  }
  return text
}


