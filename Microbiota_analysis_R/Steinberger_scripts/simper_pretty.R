"
simper.pretty: automates simper exectution for comparisons of interest

Andrew Steinberger, Kim Dill-Mcfarland, Madison Cox
asteinberger@wisc.edu
Suen Lab
University of Wisconsin-Madison

      Copyright (C) 2016 Andrew Steinberger
  
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
  
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
  
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

Simper (OTU correlations within inver-metric correlations)

Makes list with the names of all the columns from metrics file that you want to
analyze with SIMPER and combines all results into csv file

x=            OTU.table   
metrics=      metadata.table  
interesting=  list the colum headers of interest from metrics file ex. c('sloth_sp','type','sp.type')
perc_cutoff=  % cutoff desired for SIMPER output, as decimal (i.e. write 50% as 0.5)
low_cutoff=   'y' if you want to REMOVE OTUs that contribute to less than 1% of significance
low_val=      value of low cutoff (0.01)
output_name=  name to append to the simper and clean_simper output files"
###########################################################################################

simper.pretty = function(x, metrics, interesting, perc_cutoff, low_cutoff, low_val, output_name){
  library(vegan)
  #handling otu tables for taxa levels
  save=list(0)
  if(grepl("Otu", colnames(x)[1])!=TRUE){
    #converts output from A.Neumann Taxonomy script
    save=list(58)
    x=as.data.frame(t(x))
    orig_names=colnames(x)
    new_names=list()
    l=1
    for(n in colnames(x)){
      ifelse((l<10), (colnames(x)[l]=c(paste0('Otu000',c(l)))), (colnames(x)[l]=c(paste0('Otu00',c(l))))) 
      new_names=append(new_names, colnames(x)[l])
      l=l+1
    }
    orig_names=as.data.frame(orig_names, row.names = as.character(new_names))
  }
  #running simper
  for(variables in interesting){
    test_1=with(metrics, simper(x, metrics[[variables]]))
    #parsing through simper output, saving desired info to table
    for(name in names(test_1)){
      testmx=matrix(ncol=length(interesting))
      testmx=cbind(test_1[[name]]$ord,test_1[[name]]$cusum)
      sorted=testmx[order(testmx[,1]),]
      sorted=cbind(sorted,test_1[[name]]$species)
      sorted=sorted[order(sorted[,2]),]
      t=matrix(sorted[sorted[,2]<=perc_cutoff,],ncol=3)
      i=nrow(t)
      #converting percents to percent of whole
      while(i>1){
        t[i,2]=as.character(as.numeric(t[i,2])-as.numeric(t[i-1,2]))
        i=i-1
      }
      t[,1]=name
      write.table(t,file=paste(output_name,'_simper.csv',sep=""), append=TRUE, sep=",", col.names = FALSE)
    }}
  y=read.table(paste(output_name,'_simper.csv',sep=""), header=FALSE,sep=",",fill = TRUE,row.names = NULL)
  file.remove(paste(output_name,'_simper.csv',sep = ""))
  y=y[-c(1)]
  colnames(y) = c("Comparison", "SIMPER", "OTU")
  #removing results lower than low cutoff
  if(low_cutoff=='y'){
    y=y[!(as.numeric(as.character(y$SIMPER))<low_val),]
  }
  #prevents changing of colnames if OTU table
  if(58 %in% save){
    y$OTU=orig_names[match(y$OTU, rownames(orig_names)),1]
  }
  write.csv(y,file=paste(output_name,'_clean_simper.csv', sep=''))
}

#Using the function
#simper.pretty(bacteria, metrics, c('interesting'), perc_cutoff=0.5, low_cutoff = 'y', low_val=0.01, 'name')

