class MSC_TOOLS():

  def exit_number_from_outFile(outFile=None)
     fid_out=open(outFile,'r')
     for ln in fid_out:
       if (string.find(ln,'tress iteration') is not -1):
         print ln
         fid_out.close()
       elif (string.find(ln,'Exit number') is not -1):
         substr=ln[string.find(ln,'Exit number'):len(ln)]
         exitnumber=substr[12:16]
         fid_out.close()
         return exitnumber
     return -1    
