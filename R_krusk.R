'
kruskal.pretty: automates kruskal test execution on simper_pretty.R output

Andrew Steinberger
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

Changelog
8/2016     Start of script development
10/14/2016 Developed method to split Variables in Comparison column of X_clean_simper.csv
           by the "_" seperating them, storing each as variable
11/30/2016 Utilized the comparison variables to pull data out of metrics and bacteria 
           files by comparison and save as data table
12/5/2016  Utilized dplyr for otu table generation, implemented double loop system to 
           generate otu and metrics tables for variables of interest
12/6/2016  Removed multi-topic checking as SIMPER req. compar.made w/in same topic,
           kruskal.test implemented, output p.values as list, write to csv,
           First working version!! Added fdr pvalue correctn, converted to function
12/14/2016 Output includes "Taxonomy" column with the taxonomic classification of each
           otu, added breaks in nested loops to shorten processing time
12/15/2016 Fixed bug causing crash when only one topic in interesting
1/20/2017  Output includes mean relative abundances of each significant OTU for each
           of the conditions compared
1/24/2017  Handles simper.pretty output at taxonomic lvls, adds rel abund stdev to
           output file
2/28/2017  Fixed crash when rownames and csv$X unequal, optional taxonomy argument

Mothur output:
  otu=      otu table
  metrics=  metadata table
  taxonomy= .taxonomy output file from classify.otu command in mothur (optional)
simper.pretty output:
  csv=         _clean_simper.csv (*Must be imported as data.frame)
    (i.e. csv= data.frame(read.csv("PATH to .csv")))
  interesting= columns of var of interest in metadata (same as simper.pretty input)
  output_name= desired name of output (i.e. outputname_krusk_simper.csv)'
###########################################################################################

kruskal.pretty = function(otu, metrics, csv, interesting, output_name, taxonomy){
  library(vegan)
  library(dplyr)
  if(grepl("Otu", colnames(otu)[1])!=TRUE){
    #converts output from A.Neuman Taxonomy script
    otu=as.data.frame(t(otu))
  }
  #changing csv$X to rownames to allow proper splitting of comparisons
  csv$X=as.integer(rownames(csv))
  L=list()
  R=list()
  mean_L=c()
  sd_L=c()
  mean_R=c()
  sd_R=c()
  L_mean=c()
  R_mean=c()
  L_sd=c()
  R_sd=c()
  krusk=c()
  tax=c()
  L_abund=c()
  R_abund=c()
  L_abund_sd=c()
  R_abund_sd=c()
  abund=as.matrix(otu)
  abund=abund/rowSums(abund)
  for(b in levels(csv$Comparison)){
    otu_list=dplyr::filter(csv, Comparison==b) #saves otu list for current comparison
    for(i in csv$X){
      if(as.character(csv$Comparison[i])==b){  ##splitting comparisons so can call individually for table generation
        splt=as.data.frame(matrix(unlist(strsplit(as.character(csv$Comparison[i]),'_')), nrow=1, byrow=T))
        cola=as.character(splt[1,1])
        colb=as.character(splt[1,2])
        break
      }
    }
    #saving topic containing var of interest (cola/colb) (less memory intensive)
    for(topic in interesting){
      #preventing crash if there is only one topic in interesting
      if(is.null(levels(metrics[[topic]]))==TRUE){
        topic1=topic
        break
      }
      for(sbtpic in levels(metrics[[topic]])){
        if(sbtpic==cola){
          topic1=topic
          break
        }
      } 
    }
    #iterate thru rows in tpics of intrst til matches cola and colb, generates otu and metrics tbl  ##!Processing can be reduced!##
    for(rowe1 in metrics[[topic1]]){
      for(rowe2 in metrics[[topic1]]){ 
        if(rowe1==cola & rowe2==colb){ 
          listbact=otu[c(metrics[[topic1]]==cola|metrics[[topic1]]==colb),]
          listmet=metrics[c(metrics[[topic1]]==cola|metrics[[topic1]]==colb),]
          break
        }
      }
    }
    #collecting differential abundances
    sample_L=row.names(subset(metrics, metrics[[topic1]] == c(cola)))
    sample_R=row.names(subset(metrics, metrics[[topic1]] == c(colb)))
    #collecting abund values, perform/save mean and stdev calculations
    for(otus in otu_list$OTU){
      for(sample in sample_L){
        L=append(L,abund[sample,otus])
        mean_L[[otus]]=mean(as.numeric(L))
        sd_L[[otus]]=sd(as.numeric(L))
      }
      for(sample in sample_R){
        R=append(R,abund[sample,otus])
        mean_R[[otus]]=mean(as.numeric(R))
        sd_R[[otus]]=sd(as.numeric(R))
      }
      L=list()
      R=list()
    }
    #runs kruskal.test for each otu in simper csv, stores as list, also stores abundances
    for(otus in otu_list$OTU){
      result=kruskal.test(listbact[[otus]]~listmet[[topic1]])
      krusk=append(krusk, result$p.value)
      #stores taxonomic classification for each otu as list
      if(missing(taxonomy)){
        tax=append(tax, c("NA"))
      } else {
        tax=append(tax, as.character(taxonomy[otus, "Taxonomy"]))
      }
      L_mean=append(L_mean, as.character(mean_L[[otus]]))
      R_mean=append(R_mean, as.character(mean_R[[otus]]))
      L_sd=append(L_sd, as.character(sd_L[[otus]]))
      R_sd=append(R_sd, as.character(sd_R[[otus]]))
      
    }
  }
  #adjusted p-values for multiple comparisons
  fdr=p.adjust(krusk, method='fdr')
  #order csv to match 'krusk'/'fdr' list, add p.val, add taxonomy, re-ord to match orig csv, write to csv
  o_csv=dplyr::arrange(csv, Comparison)
  o_csv[,5]=krusk
  o_csv[,6]=fdr
  o_csv[,7]=tax
  o_csv[,8]=L_mean
  o_csv[,9]=L_sd
  o_csv[,10]=R_mean
  o_csv[,11]=R_sd
  o_csv=dplyr::arrange(o_csv, X)
  colnames(o_csv)[which(names(o_csv) == "V5")] <- "krusk_p.val" #changes column header
  colnames(o_csv)[which(names(o_csv) == "V6")] <- "fdr_krusk_p.val"
  colnames(o_csv)[which(names(o_csv) == "V7")] <- "Taxonomy"
  colnames(o_csv)[which(names(o_csv) == "V8")] <- "Left mean abund"
  colnames(o_csv)[which(names(o_csv) == "V9")] <- "Left stdev"
  colnames(o_csv)[which(names(o_csv) == "V10")] <- "Right mean abund"
  colnames(o_csv)[which(names(o_csv) == "V11")] <- "Right stdev"
  o_csv[,1]=NULL
  write.csv(o_csv, file=paste(output_name,"_krusk_simper.csv", sep=""))
}

#example of use
#kruskal.pretty(bacteria, metrics, csv, c('sloth_sp','type','sp.type'), 'sloth', taxonomy)
