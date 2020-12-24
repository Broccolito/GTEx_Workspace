library(dplyr)

split_gtex = function(gtex_data_dir, gtex_tissue){
  
  cat(paste0("Reading GTEx file from ", gtex_tissue, " ... \n"))
  gtex_data_dir = gzfile(gtex_data_dir)
  gtex_data = read.delim(gtex_data_dir)
  
  cat("Splitting GTEx file by Chromosomes ... \n")
  chr = unlist(lapply(strsplit(gtex_data$variant_id, split = "_"),
                      function(x){return(x[1])}))
  
  cat("Annotating GTEx file ... \n")
  gtex_data = gtex_data %>%
    mutate(chr = chr)
  
  cat("Generating Split Objects ... \n")
  gtex_data_chr1 = filter(gtex_data,chr == "chr1")
  gtex_data_chr2 = filter(gtex_data,chr == "chr2")
  gtex_data_chr3 = filter(gtex_data,chr == "chr3")
  gtex_data_chr4 = filter(gtex_data,chr == "chr4")
  gtex_data_chr5 = filter(gtex_data,chr == "chr5")
  gtex_data_chr6 = filter(gtex_data,chr == "chr6")
  gtex_data_chr7 = filter(gtex_data,chr == "chr7")
  gtex_data_chr8 = filter(gtex_data,chr == "chr8")
  gtex_data_chr9 = filter(gtex_data,chr == "chr9")
  gtex_data_chr10 = filter(gtex_data,chr == "chr10")
  gtex_data_chr11 = filter(gtex_data,chr == "chr11")
  gtex_data_chr12 = filter(gtex_data,chr == "chr12")
  gtex_data_chr13 = filter(gtex_data,chr == "chr13")
  gtex_data_chr14 = filter(gtex_data,chr == "chr14")
  gtex_data_chr15 = filter(gtex_data,chr == "chr15")
  gtex_data_chr16 = filter(gtex_data,chr == "chr16")
  gtex_data_chr17 = filter(gtex_data,chr == "chr17")
  gtex_data_chr18 = filter(gtex_data,chr == "chr18")
  gtex_data_chr19 = filter(gtex_data,chr == "chr19")
  gtex_data_chr20 = filter(gtex_data,chr == "chr20")
  gtex_data_chr21 = filter(gtex_data,chr == "chr21")
  gtex_data_chr22 = filter(gtex_data,chr == "chr22")
  
  cat("Saving Split Files ... \n")
  save(gtex_data_chr1,file = paste0(gtex_tissue,"_chr1.RData"))
  save(gtex_data_chr2,file = paste0(gtex_tissue,"_chr2.RData"))
  save(gtex_data_chr3,file = paste0(gtex_tissue,"_chr3.RData"))
  save(gtex_data_chr4,file = paste0(gtex_tissue,"_chr4.RData"))
  save(gtex_data_chr5,file = paste0(gtex_tissue,"_chr5.RData"))
  save(gtex_data_chr6,file = paste0(gtex_tissue,"_chr6.RData"))
  save(gtex_data_chr7,file = paste0(gtex_tissue,"_chr7.RData"))
  save(gtex_data_chr8,file = paste0(gtex_tissue,"_chr8.RData"))
  save(gtex_data_chr9,file = paste0(gtex_tissue,"_chr9.RData"))
  save(gtex_data_chr10,file = paste0(gtex_tissue,"_chr10.RData"))
  save(gtex_data_chr11,file = paste0(gtex_tissue,"_chr11.RData"))
  save(gtex_data_chr12,file = paste0(gtex_tissue,"_chr12.RData"))
  save(gtex_data_chr13,file = paste0(gtex_tissue,"_chr13.RData"))
  save(gtex_data_chr14,file = paste0(gtex_tissue,"_chr14.RData"))
  save(gtex_data_chr15,file = paste0(gtex_tissue,"_chr15.RData"))
  save(gtex_data_chr16,file = paste0(gtex_tissue,"_chr16.RData"))
  save(gtex_data_chr17,file = paste0(gtex_tissue,"_chr17.RData"))
  save(gtex_data_chr18,file = paste0(gtex_tissue,"_chr18.RData"))
  save(gtex_data_chr19,file = paste0(gtex_tissue,"_chr19.RData"))
  save(gtex_data_chr20,file = paste0(gtex_tissue,"_chr20.RData"))
  save(gtex_data_chr21,file = paste0(gtex_tissue,"_chr21.RData"))
  save(gtex_data_chr22,file = paste0(gtex_tissue,"_chr22.RData"))
  
  cat("Done! \n")

}

# split_gtex("GTEx_Analysis_v8_QTLs_GTEx_Analysis_v8_eQTL_all_associations_Artery_Aorta.allpairs.txt",
#            gtex_tissue = "Artery_Aorta")

split_gtex("GTEx_Analysis_v8_QTLs_GTEx_Analysis_v8_eQTL_all_associations_Artery_Coronary.allpairs.txt",
           gtex_tissue = "Artery_Coronary")

split_gtex("GTEx_Analysis_v8_QTLs_GTEx_Analysis_v8_eQTL_all_associations_Artery_Tibial.allpairs.txt",
           gtex_tissue = "Artery_Tibial")

split_gtex("GTEx_Analysis_v8_QTLs_GTEx_Analysis_v8_eQTL_all_associations_Heart_Atrial_Appendage.allpairs.txt",
           gtex_tissue = "Heart_Atrial_Appendage")

split_gtex("GTEx_Analysis_v8_QTLs_GTEx_Analysis_v8_eQTL_all_associations_Heart_Left_Ventricle.allpairs.txt",
           gtex_tissue = "Heart_Left_Ventricle")

split_gtex("GTEx_Analysis_v8_QTLs_GTEx_Analysis_v8_eQTL_all_associations_Lung.allpairs.txt",
           gtex_tissue = "Lung")

split_gtex("GTEx_Analysis_v8_QTLs_GTEx_Analysis_v8_eQTL_all_associations_Muscle_Skeletal.allpairs.txt",
           gtex_tissue = "_Muscle_Skeletal")

split_gtex("GTEx_Analysis_v8_QTLs_GTEx_Analysis_v8_eQTL_all_associations_Spleen.allpairs.txt",
           gtex_tissue = "Spleen")

split_gtex("GTEx_Analysis_v8_QTLs_GTEx_Analysis_v8_eQTL_all_associations_Whole_Blood.allpairs.txt",
           gtex_tissue = "Whole_Blood")




