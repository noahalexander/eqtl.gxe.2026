library(stringr)




#note that here the crosses to parents is slightly different than in prior versions
crosses.to.parents=list(
  '375'=c("M22", "BYa"),           #1
  '3051'  =c("BYa", "RMx"),           #2
  '376'=c("RMx", "YPS163a"),       #3
  'B'  =c("YPS163a", "YJM145x"),   #4
  '377'=c("YJM145x", "CLIB413a"),  #5
  '393'=c("CLIB413a", "YJM978x"),  #6
  '381'=c("YJM978x", "YJM454a"),   #7
  '3008'=c("YJM454a", "YPS1009x"),  #8
  '2999'=c("YPS1009x", "I14a"),     #9
  '3000'=c("I14a", "Y10x"),         #10
  '3001'=c("Y10x", "PW5a"),         #11
  '3049'=c("PW5a", "273614xa"),     #12
  '3003'=c("273614xa", "YJM981x"),  #13
  '3004'=c("YJM981x", "CBS2888a"),  #14
  '3043'=c("CBS2888a", "CLIB219x"), #15
  '3028'=c("CLIB219x", "M22")       #16
)

vcf_annotation_16 = readRDS("/Users/noahalexander/vcf_annotation_16.20240110.RDS")

#combined_dir = "/Volumes/Extreme SSD/data_zip/data/out/combined"
#cross="3004"
#perturbation = "nacl"
#timepoint="t30"

make_index_gr <- function(combined_dir,cross,perturbation,timepoint) {
  sdf = readRDS(paste0(combined_dir,"/",cross,"/", perturbation,"/",timepoint,"/","segData.RDS"))
  Gsub = sdf$Gsub
  #start with column names of gsub and strip out the sequence position
  gcols = colnames(Gsub)
  index = data.frame(chr=word(gcols,1,sep = "\\_"), start= word(gcols,2,sep = "\\_"))
  #make a granges using the positions of each marker (not ranges but positions). 
  gr.index= with(index, GRanges(chr, IRanges(start=start)))
  return(gr.index)
}


make_state_granges <- function(combined_dir,cross,perturbation,timepoint) {
  
  df = readRDS(paste0(combined_dir,"/",cross,"/", perturbation,"/",timepoint,"/","state_pheno_LOD_peaks.RDS"))
  
  #subset for rp assignments 
  dfrp = df[[1]]
  #subset for iesr assignments 
  dfie = df[[2]]
  #subset for ribi assignments 
  dfrb = df[[3]]
  
  ##make data frame for each state and make a granges 
  drp = data.frame(chr=word(dfrp$fscan.markers,1,sep = "\\_"), 
    start= sapply(strsplit(dfrp$CI.l, "_"), function(x) x[2]), 
    end= sapply(strsplit(dfrp$CI.r, "_"), function(x) x[2] ),
    p=dfrp$p,
    q=dfrp$q,
    LOD=dfrp$LOD,
    fscan.markers= dfrp$fscan.markers
  )
  drp$start = as.numeric(drp$start)
  #rangesB$start= as.numeric(rangesB$start)
  drp$end = as.numeric(drp$end)
  
  drb = data.frame(chr=word(dfrb$fscan.markers,1,sep = "\\_"), 
    start= sapply(strsplit(dfrb$CI.l, "_"), function(x) x[2]), 
    end= sapply(strsplit(dfrb$CI.r, "_"), function(x) x[2] ),
    p=dfrb$p,
    q=dfrb$q,
    LOD=dfrb$LOD,
    fscan.markers= dfrb$fscan.markers
  )
  drb$start = as.numeric(drb$start)
  #rangesB$start= as.numeric(rangesB$start)
  drb$end = as.numeric(drb$end)
  
  die = data.frame(chr=word(dfie$fscan.markers,1,sep = "\\_"), 
    start= sapply(strsplit(dfie$CI.l, "_"), function(x) x[2]), 
    end= sapply(strsplit(dfie$CI.r, "_"), function(x) x[2] ),
    p=dfie$p,
    q=dfie$q,
    LOD=dfie$LOD,
    fscan.markers= dfie$fscan.markers
  )
  die$start = as.numeric(die$start)
  die$end = as.numeric(die$end)
  
  
  #make granges object
  gr.rp= with(drp, GRanges(chr, IRanges(start=start, end=end), LOD=LOD, p=p, q=q))
  gr.rp= gr.rp[order(gr.rp$LOD, decreasing = T),]
  
  gr.rb= with(drb, GRanges(chr, IRanges(start=start, end=end), LOD=LOD, p=p, q=q))
  gr.rb= gr.rb[order(gr.rb$LOD, decreasing = T),]
  
  gr.ie= with(die, GRanges(chr, IRanges(start=start, end=end), LOD=LOD, p=p, q=q))
  gr.ie= gr.ie[order(gr.ie$LOD, decreasing = T),]
  
  granges_list = list(gr.rb, gr.rp, gr.ie)
  return(granges_list)
  
}

#expand the range by one marker in both directions 
expand_interval <- function(gr.state, gr.index ) {
  #############---------------#################-------------------###################################3004 cross NaCl timecourse 
  #read in segdata and make granges object using gsub. this can be used for all time points in this experiment   
  #look for markers from index granges that are to the left or right of the input interval (state qtl CI)
  #the output is an index that can be used to search the index granges for the adjacent markers to widen the interval 
  #these are indices, not bp coordinates but they can be converted using index granges 
  nearest.markers = nearest(gr.state, gr.index)
  right.markers = precede(gr.state, gr.index)
  # vector of indices that precede the interval in spite of the function name 
  left.markers = follow(gr.state, gr.index)
  #add positions of markers on either side of returned interval 
  gr.state$left.markers = left.markers
  gr.state$nearest = nearest.markers
  gr.state$right.markers = right.markers
  #pull out the bp coordinates using the indices 
  index = data.frame(gr.index)
  
  isub.left = index[left.markers,]
  left.coorinate = as.numeric(isub.left$start)
  isub.nearest = index[nearest.markers,]
  nearest.coordinate = as.numeric(isub.nearest$start)
  isub.right = index[right.markers,]
  right.coorinate = as.numeric(isub.right$start)
  
  #add bp coordinates on either side of Josh's returned interval 
  gr.state$left.coorinate = left.coorinate
  gr.state$nearest.coordinate = nearest.coordinate
  gr.state$right.coorinate = right.coorinate
  
  #if not containing an na, the difference should reflect the size of the new CI
  gr.state$new_CI_length = right.coorinate - left.coorinate
  
  #deal with NAs and remake ranges to reflect the new interval
  df = as.data.frame(gr.state)
  for(i in 1:nrow(df)) {
    row <- df[i,]
    #print(row[2])
    if (is.na(row[12]) == T) {
      #set left.marker to the original start of the interval
      df[i,12] = df[i,2]
      #print(row[9])
    }
    if (is.na(row[14]) == T) {
      #set right.marker to the original end of the interval
      df[i,14] = df[i,3]
      #print(row[9])
    }
  }
  
  #remake granges adding the new/checked broader CIs
  #gr.state= with(drp.3004.nacl.t0, GRanges(chr, IRanges(start=start, end=end), LOD=LOD, p=p, q=q))
  #using column names for new data frame, make granges with wider range
  gr.state = with(df, GRanges(seqnames = seqnames, IRanges(start=left.coorinate, end=right.coorinate), LOD=LOD, p=p, q=q, left.markers = left.markers, nearest.markers=nearest, right.markers=right.markers, original.start=start, original.end=end))
  df = as.data.frame(gr.state)
  df$interval_sgd_browser = paste0(df$seqnames,":",df$start,"-",df$end)
  new_colnames = c("seqnames" ,"start" , "end" ,"width" ,  "LOD" , "interval_sgd_browser" ,"p",  "q" ,"left.markers","nearest.markers","right.markers","original.start","original.end","strand" )
  df = df[,new_colnames]
  dl= list(df, gr.state)
  
  return(dl)
  
}



#all intervals are in same spreadsheet in the output produced by the following function but 
#a mod could be making separate sheets for each interval...current way seems fine and it is easy to control f
#at the current few mb file sizes 

excel_region_anno_output <- function(cross, input_granges){
  
  parents =crosses.to.parents[[cross]]
  print(parents)
  #causal_hit= causal_hits%>% filter(shared.parent  %in% parents)  %>% mutate(seqnames=chr) %>% as_granges()
  #causal_hit = causal_hit_3004 %>% mutate(seqnames=chr) %>% as_granges()
  #qtl = qtl_bloom %>% filter(cross == !!cross)
  vcf_names = str_remove_all(parents,pattern="a|x")
  vcf_annotation = vcf_annotation_16 %>% filter((grepl(vcf_names[1],ref_strains) & grepl(vcf_names[2],alt_strains)) |  (grepl(vcf_names[2],ref_strains) & grepl(vcf_names[1],alt_strains)))
  ## Missense of STOP
  size = unlist(lapply(str_split(vcf_annotation$mutation,"\\|"), function(x){x[3]}))
  vcf_annotation$impact  = size 
  
  sub_anno = input_granges %>% join_overlap_intersect(vcf_annotation) %>% as_data_frame()
  openxlsx::write.xlsx(sub_anno,file=paste0(deparse(substitute(cross)), deparse(substitute(input_granges)), "intersect.xlsx"))
}



vcf_anno_gen <- function(cross){
  
  
  parents =crosses.to.parents[[cross]]
  #causal_hit= causal_hits%>% filter(shared.parent  %in% parents)  %>% mutate(seqnames=chr) %>% as_granges()
  #causal_hit = causal_hit_3004 %>% mutate(seqnames=chr) %>% as_granges()
  #qtl = qtl_bloom %>% filter(cross == !!cross)
  vcf_names = str_remove_all(parents,pattern="a|x")
  vcf_annotation = vcf_annotation_16 %>% filter((grepl(vcf_names[1],ref_strains) & grepl(vcf_names[2],alt_strains)) |  (grepl(vcf_names[2],ref_strains) & grepl(vcf_names[1],alt_strains)))
  ## Missense of STOP
  size = unlist(lapply(str_split(vcf_annotation$mutation,"\\|"), function(x){x[3]}))
  vcf_annotation$impact  = size 
  #sub_anno = input_granges %>% join_overlap_intersect(vcf_annotation) %>% as_data_frame()
 return(vcf_annotation)
   
}



#####make bed from granges 

make_bed_from_granges <- function(input.granges) {
  
  df = as.data.frame(input.granges)
  print(head(df))
  dp= data.frame(chr=df$seqnames, start=df$start, end= df$end, freq = df$freq)
  #daf= data.frame(chr=df$chrom, start=df$pos, end= df$pos, af=df$joseph_af)
  head(dp)
  #dp=subset(dp, is.na(dp$prov)==F)
  #dp=subset(dp, abs(dp$prov) >=3)
  write.table(dp, file=paste0(deparse(substitute(input.granges)),".bed"), quote=F, sep="\t", row.names=F, col.names=T)
  
}



########################################---------------------------------------------------------



#for testing 
#gr.state = gr.list.3004.t0[[2]]

#expand the range by one marker in both directions 
expand_interval3 <- function(gr.state, gr.index ) {
  #############---------------#################-------------------###################################3004 cross NaCl timecourse 
  #read in segdata and make granges object using gsub. this can be used for all time points in this experiment   
  #look for markers from index granges that are to the left or right of the input interval (state qtl CI)
  #the output is an index that can be used to search the index granges for the adjacent markers to widen the interval 
  #these are indices, not bp coordinates but they can be converted using index granges 
  nearest.markers = nearest(gr.state, gr.index)
  right.markers = precede(gr.state, gr.index)
  # vector of indices that precede the interval in spite of the function name 
  left.markers = follow(gr.state, gr.index)
  #add positions of markers on either side of returned interval 
  gr.state$left.markers = left.markers
  gr.state$nearest = nearest.markers
  gr.state$right.markers = right.markers
  #pull out the bp coordinates using the indices 
  index = data.frame(gr.index)
  
  isub.left = index[left.markers,]
  left.coorinate = as.numeric(isub.left$start)
  isub.nearest = index[nearest.markers,]
  nearest.coordinate = as.numeric(isub.nearest$start)
  isub.right = index[right.markers,]
  right.coorinate = as.numeric(isub.right$start)
  
  #add bp coordinates on either side of Josh's returned interval 
  gr.state$left.coorinate = left.coorinate
  gr.state$nearest.coordinate = nearest.coordinate
  gr.state$right.coorinate = right.coorinate
  
  #if not containing an na, the difference should reflect the size of the new CI
  gr.state$new_CI_length = right.coorinate - left.coorinate
  
 
###note: haven't done anything to the range yet, just calculating positions that the intervals could be expanded to  

#now expand in a simple way 
##-------------------check if interval is:
  #1)a single base, if so, expand to nearest markers on either side, otherwise check if
  #2)under 15k, if it is keep it, otherwise
  #3)create new, smaller interval from midpoint of original interval extended by 7.5k on each side 
  
  #deal with NAs and remake ranges to reflect the new interval
  df = as.data.frame(gr.state)
  df = subset(df, df$LOD > 10)
  
  #fix NAs from chr ends 
  for(i in 1:nrow(df)) {
    row <- df[i,]
  
    if (is.na(row[12]) == T) {
      #set left.marker to the original start of the interval
      df[i,12] = df[i,2]
      #print(row[9])
    }
    
    if (is.na(row[14]) == T) {
      #set right.marker to the original end of the interval
      df[i,14] = df[i,3]
      #print(row[9])
    }
    
    if (row$width == 1) {
      df$new.start[i]=df$left.coorinate[i]
      df$new.end[i]=df$right.coorinate[i]
      #print(row[2])
    } else {
      df$new.start[i]=df$start[i]
      df$new.end[i]=df$end[i]
    }
    df$new.width[i] = df$new.end[i]-df$new.start[i]
      
  }
  
  #dfo = df
  #reduce intervals above 15k to 15k as an ad hoc way to cut them down a bit 
  for(i in 1:nrow(df)) {
    row <- df[i,]
    
    if (row$new.width > 15000) {
      #find midpoint of new start and end coordinates and establish a new start coordinate that is 7.5k upstream of that point 
      midpoint = (df$new.start[i]+((df$new.end[i]-df$new.start[i])/2))
      df$new.start[i]=midpoint -7500
      #find midpoint of new start and end coordinates and establish a new end coordinate that is 7.5k downstream of that point 
      df$new.end[i]=midpoint +7500
      #print(row[2])
    } else {
      df$new.start[i]=df$new.start[i]
      df$new.end[i]=df$new.end[i]
    }
    df$new.width[i] = df$new.end[i]-df$new.start[i]
    
  }
  
  
  #clean up NAs again 
  for(i in 1:nrow(df)) {
    row <- df[i,]
    
    if (is.na(row[12]) == T) {
      #set left.marker to the original start of the interval
      df[i,12] = df[i,2]
      #print(row[9])
    }
    
    if (is.na(row[14]) == T) {
      #set right.marker to the original end of the interval
      df[i,14] = df[i,3]
      #print(row[9])
    }
    
  }
  #remake granges adding the new/checked broader CIs
  #gr.state= with(drp.3004.nacl.t0, GRanges(chr, IRanges(start=start, end=end), LOD=LOD, p=p, q=q))
  #using column names for new data frame, make granges with wider range
  gr.state = with(df, GRanges(seqnames = seqnames, IRanges(start=new.start, end=new.end), LOD=LOD, p=p, q=q, interval_sgd_browser = paste0(df$seqnames,":",df$new.start,"-",df$new.end), left.markers = left.markers, nearest.markers=nearest, right.markers=right.markers, original.start=start, original.end=end))
  df = as.data.frame(gr.state)
  df$interval_sgd_browser = paste0(df$seqnames,":",df$new.start,"-",df$new.end)
  new_colnames = c("seqnames" ,"start" , "end" ,"width" ,  "LOD" , "interval_sgd_browser" ,"p",  "q" ,"left.markers","nearest.markers","right.markers","original.start","original.end","strand" )
  df = df[,new_colnames]
  dl= list(df, gr.state)
  
  return(dl)
  
}
