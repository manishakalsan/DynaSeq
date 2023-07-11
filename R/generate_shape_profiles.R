# DynaSeq
# .
# .
# It predicts molecular dynamics-derived ensembles of 
# 13 DNA conformational parameters namely Shear, Stretch, 
# Stagger, Buckle, Prop-Tw, Opening, Shift, Slide, Rise, 
# Tilt, Roll, Twist, and MGW.
# .
# .
# takes a list of fasta files and generates shape profiles
# .
# .
# .
# @USAGE: generate_shape(list_files)
# @parameter list_files A list containing names of .fasta files
# @return writes a tab separated file containing DNA conformational 
# parameters (static and ensemble)
# requires rds files in the working directory
# 
# the function will report DNA shape profiles for DNA sequences without N. 
# Those containing N will be removed
# 
# @example
# test_fasta<-system("ls test*fasta",intern=T)
# generate_shape(test_fasta)
#

generate_shape<-function(list_files){
   norm_bin<-readRDS("dictionary_ensemble_5bin_5mer.rds")
   avg_shape<-readRDS("dictionary_static_5mer.rds")
   kmer_size<-5
   all_shape<-lapply(1:length(list_files),function(x){
      shape_res<-c()
      data<-read.table(list_files[x],header=F)
      r_seqid<-as.data.frame(data[seq(1,dim(data)[1],2),])
      #remove duplicate sequences if any
      uniq_seqid<-unique(r_seqid[,1]) #unique genomic indices
      r_seq<-lapply(1:length(uniq_seqid), function(id){
         genomic_seq<-data[(grep(uniq_seqid[id],data[,1]))+1,]
         if(length(genomic_seq) == 1){
            genomic_seq<-genomic_seq
         } else {
            genomic_seq<-genomic_seq[1]
         }
         genomic_seq
      })
      r_seq<-as.data.frame(unlist(r_seq))
      #remove sequences with N 
      if(length(grep("N", r_seq[,1]))>=1){
         n_seq<-grep("N", r_seq[,1]) # sequences with N
         r_seq<-as.data.frame(r_seq[-(n_seq),1]) 
         uniq_seqid<-uniq_seqid[-n_seq]
      }
      #generating pentamers of all sequences
      all_seq_sub_shape<-lapply(1:dim(r_seq)[1],function(y){
         aseq<-as.character(r_seq[y,1])
         aseq_sub<-lapply(1:(nchar(aseq)-(4)),function(l){
            one_sub<-toupper(substr(aseq,l,as.numeric(l+(4))))
            one_sub
         })
         aseq_sub<-unlist(aseq_sub)
         #generating DNA shape profiles for all sub sequences
         one_seq_shape<-lapply(1:length(aseq_sub),function(seq_id){
            ens_data<-norm_bin[grep(toupper(aseq_sub[seq_id]),rownames(norm_bin)),]
            avg_data<-avg_shape[grep(toupper(aseq_sub[seq_id]),rownames(avg_shape)),]
            seq_shape<-c(avg_data,ens_data)
            seq_shape
         })
         one_seq_shape<-do.call(rbind,one_seq_shape)
         #naming of rows
         shape_pos<-c(3:(nchar(aseq)-2))
         all_rows<-lapply(1:length(shape_pos), function(pos){
            pos_index<-paste(as.character(r_seqid[y,1]),"_pos",shape_pos[pos],"_",aseq_sub[pos],sep="")
            pos_index
         })
         rownames(one_seq_shape)<-all_rows
         one_seq_shape
      })
      shape_res<-do.call(rbind,all_seq_sub_shape)
      write.table(shape_res,paste("5bin_dna_shape_",strsplit(list_files[x],".fasta")[[1]][1],sep=""),quote=F)
   })
}

