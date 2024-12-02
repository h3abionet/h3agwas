def configfile_analysis(file, sep){
   File theInfoFile = new File( file )
   if(file.contains("s3://")){
     println "File " + file + "is in S3 so its existence can't be checked."
   }
   if(file.contains("az://")){
     println "File " + file + "is in Azure so its existence can't be checked."
   }
   if( !theInfoFile.exists() && !file.contains("s3://") && !file.contains("az://") ) {
      println "File "+ file+" does not exist"
      exit(0)
   }
   def lines = theInfoFile.readLines()
   def SplH=lines[0].split(sep)
   if(SplH.size()==1){
      println SplH
      println "problem of separator"
      exit(1)
   }
   resInfo=[]
   resFile=[]
   resIndex=[]
   lines.remove(0)
   PosCmp=-1
   CmtL=0
   listnum=[]
   NumRef=-1
   for(line in lines){
      def splLine=line.split(sep)
      if(splLine.size()>1){
       def cmtelem=0
       def SubRes=[]
       while (cmtelem < SplH.size()){
            if(SplH[cmtelem]=='File'){
               resFile.add(splLine[cmtelem])
            }
            if(SplH[cmtelem]=='Index'){
               resIndex.add(splLine[cmtelem])
            }
            else if(SplH[cmtelem]=='IsRefFile' && splLine[cmtelem]=='1'){
               NumRef=CmtL
               cmtelem2=0
            }
            if(splLine[cmtelem]!='NA' && splLine[cmtelem]!='' && SplH[cmtelem]!='File' && SplH[cmtelem]!='IsRefFile'){
                 SubRes.add(SplH[cmtelem]+':'+splLine[cmtelem])
            }
            cmtelem+=1
       }
       resInfo.add(SubRes.join(','))
       listnum.add(CmtL)
       CmtL+=1
     }
   }
 if(CmtL<2){
   println "no enough file found for meta analyse, check your column header and sep in your csv file\n exit"
   exit(1)
 }
 if(resIndex.size()==0)  resIndex=listnum
 return([resFile,resInfo, resIndex, listnum])
}
 

