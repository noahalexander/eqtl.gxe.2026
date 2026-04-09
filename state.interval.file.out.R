##--------3004 NaCl t0
combined_dir = "/Users/noahalexander/repeat_fine_mapping/combined/"
cross="3004"
perturbation = "nacl"
timepoint="t0"

#make index granges of variant positions and extract/generate granges of state qtl confidence intervals 
gr.3004.t0.index = make_index_gr(combined_dir, cross, perturbation, timepoint)
gr.list.3004.t0 = make_state_granges(combined_dir, cross, perturbation, timepoint)

#make state granges objects for each of the iesr states with expanded confidence intervals 
#(here by 1 marker in each direction unless the existing border of the CI is chr end)
gr.rp.3004.nacl.t0 = expand_interval3((gr.list.3004.t0[[1]]), gr.3004.t0.index)
gr.rb.3004.nacl.t0 = expand_interval3(gr.list.3004.t0[[2]], gr.3004.t0.index)
gr.ie.3004.nacl.t0 = expand_interval3(gr.list.3004.t0[[3]], gr.3004.t0.index)


##--potentially also subset for regions that intersect with hotspots in relevant experiments 

#make excel output of subset of vcf for parents and CIs. 
excel_region_anno_output(cross = "3004", gr.rp.3004.nacl.t0[[2]])
excel_region_anno_output(cross = "3004", gr.rb.3004.nacl.t0[[2]])
excel_region_anno_output(cross = "3004", gr.ie.3004.nacl.t0[[2]])


#make bed tracks
dd = data.frame(chr = gr.ie.3004.nacl.t0[[1]]$seqnames, start = gr.ie.3004.nacl.t0[[1]]$start, end=gr.ie.3004.nacl.t0[[1]]$end)
dd$start = as.numeric(dd$start)
dd$end = as.numeric(dd$end)
write.table(dd, file="3004.nacl.t0.ie.aucell.CIs.bed", quote=F, sep="\t", row.names=F, col.names=T)

dd = data.frame(chr = gr.rp.3004.nacl.t0[[1]]$seqnames, start = gr.rp.3004.nacl.t0[[1]]$start, end=gr.rp.3004.nacl.t0[[1]]$end)
dd$start = as.numeric(dd$start)
dd$end = as.numeric(dd$end)
write.table(dd, file="3004.nacl.t0.rp.aucell.CIs.bed", quote=F, sep="\t", row.names=F, col.names=T)

dd = data.frame(chr = gr.rb.3004.nacl.t0[[1]]$seqnames, start = gr.rb.3004.nacl.t0[[1]]$start, end=gr.rb.3004.nacl.t0[[1]]$end)
dd$start = as.numeric(dd$start)
dd$end = as.numeric(dd$end)
write.table(dd, file="3004.nacl.t0.rb.aucell.CIs.bed", quote=F, sep="\t", row.names=F, col.names=T)





##--------3004 NaCl t30
combined_dir = "/Users/noahalexander/repeat_fine_mapping/combined/"
cross="3004"
perturbation = "nacl"
timepoint="t30"

#make index granges of variant positions and extract/generate granges of state qtl confidence intervals 
gr.3004.t30.index = make_index_gr(combined_dir, cross, perturbation, timepoint)
gr.list.3004.t30 = make_state_granges(combined_dir, cross, perturbation, timepoint)

#make state granges objects for each of the iesr states with expanded confidence intervals 
#(here by 1 marker in each direction unless the existing border of the CI is chr end)
gr.rp.3004.nacl.t30 = expand_interval3((gr.list.3004.t30[[1]]), gr.3004.t30.index)
gr.rb.3004.nacl.t30 = expand_interval3(gr.list.3004.t30[[2]], gr.3004.t30.index)
gr.ie.3004.nacl.t30 = expand_interval3(gr.list.3004.t30[[3]], gr.3004.t30.index)

#make excel output of subset of vcf for parents and CIs. 
excel_region_anno_output(cross = "3004", gr.rp.3004.nacl.t30[[2]])
excel_region_anno_output(cross = "3004", gr.rb.3004.nacl.t30[[2]])
excel_region_anno_output(cross = "3004", gr.ie.3004.nacl.t30[[2]])



#make bed tracks
dd = data.frame(chr = gr.ie.3004.nacl.t30[[1]]$seqnames, start = gr.ie.3004.nacl.t30[[1]]$start, end=gr.ie.3004.nacl.t30[[1]]$end)
dd$start = as.numeric(dd$start)
dd$end = as.numeric(dd$end)
write.table(dd, file="3004.nacl.t30.ie.aucell.CIs.bed", quote=F, sep="\t", row.names=F, col.names=T)

dd = data.frame(chr = gr.rp.3004.nacl.t30[[1]]$seqnames, start = gr.rp.3004.nacl.t30[[1]]$start, end=gr.rp.3004.nacl.t30[[1]]$end)
dd$start = as.numeric(dd$start)
dd$end = as.numeric(dd$end)
write.table(dd, file="3004.nacl.t30.rp.aucell.CIs.bed", quote=F, sep="\t", row.names=F, col.names=T)

dd = data.frame(chr = gr.rb.3004.nacl.t30[[1]]$seqnames, start = gr.rb.3004.nacl.t30[[1]]$start, end=gr.rb.3004.nacl.t30[[1]]$end)
dd$start = as.numeric(dd$start)
dd$end = as.numeric(dd$end)
write.table(dd, file="3004.nacl.t30.rb.aucell.CIs.bed", quote=F, sep="\t", row.names=F, col.names=T)


##--------3051 NaCl t0

combined_dir = "/Users/noahalexander/repeat_fine_mapping/combined/"
cross="A"
perturbation = "nacl"
timepoint="t0"

#make index granges of variant positions and extract/generate granges of state qtl confidence intervals 
gr.3051.t0.index = make_index_gr(combined_dir, cross, perturbation, timepoint)
gr.list.3051.t0 = make_state_granges(combined_dir, cross, perturbation, timepoint)

#make state granges objects for each of the iesr states with expanded confidence intervals 
#(here by 1 marker in each direction unless the existing border of the CI is chr end)
gr.rp.3051.nacl.t0 = expand_interval3((gr.list.3051.t0[[1]]), gr.3051.t0.index)
gr.rb.3051.nacl.t0 = expand_interval3(gr.list.3051.t0[[2]], gr.3051.t0.index)
gr.ie.3051.nacl.t0 = expand_interval3(gr.list.3051.t0[[3]], gr.3051.t0.index)

#make excel output of subset of vcf for parents and CIs. 
excel_region_anno_output(cross = "3051", gr.rp.3051.nacl.t0[[2]])
excel_region_anno_output(cross = "3051", gr.rb.3051.nacl.t0[[2]])
excel_region_anno_output(cross = "3051", gr.ie.3051.nacl.t0[[2]])



#make bed tracks
dd = data.frame(chr = gr.ie.3051.nacl.t0[[1]]$seqnames, start = gr.ie.3051.nacl.t0[[1]]$start, end=gr.ie.3051.nacl.t0[[1]]$end)
dd$start = as.numeric(dd$start)
dd$end = as.numeric(dd$end)
write.table(dd, file="3051.nacl.t0.ie.aucell.CIs.bed", quote=F, sep="\t", row.names=F, col.names=T)

dd = data.frame(chr = gr.rp.3051.nacl.t0[[1]]$seqnames, start = gr.rp.3051.nacl.t0[[1]]$start, end=gr.rp.3051.nacl.t0[[1]]$end)
dd$start = as.numeric(dd$start)
dd$end = as.numeric(dd$end)
write.table(dd, file="3051.nacl.t0.rp.aucell.CIs.bed", quote=F, sep="\t", row.names=F, col.names=T)

dd = data.frame(chr = gr.rb.3051.nacl.t0[[1]]$seqnames, start = gr.rb.3051.nacl.t0[[1]]$start, end=gr.rb.3051.nacl.t0[[1]]$end)
dd$start = as.numeric(dd$start)
dd$end = as.numeric(dd$end)
write.table(dd, file="3051.nacl.t0.rb.aucell.CIs.bed", quote=F, sep="\t", row.names=F, col.names=T)

##--------3051 NaCl t30
combined_dir = "/Users/noahalexander/repeat_fine_mapping/combined/"
cross="A"
perturbation = "nacl"
timepoint="t30"

#make index granges of variant positions and extract/generate granges of state qtl confidence intervals 
gr.3051.t30.index = make_index_gr(combined_dir, cross, perturbation, timepoint)
gr.list.3051.t30 = make_state_granges(combined_dir, cross, perturbation, timepoint)

#make state granges objects for each of the iesr states with expanded confidence intervals 
#(here by 1 marker in each direction unless the existing border of the CI is chr end)
gr.rp.3051.nacl.t30 = expand_interval3((gr.list.3051.t30[[1]]), gr.3051.t30.index)
gr.rb.3051.nacl.t30 = expand_interval3(gr.list.3051.t30[[2]], gr.3051.t30.index)
gr.ie.3051.nacl.t30 = expand_interval3(gr.list.3051.t30[[3]], gr.3051.t30.index)

#make excel output of subset of vcf for parents and CIs. 
excel_region_anno_output(cross = "3051", gr.rp.3051.nacl.t30[[2]])
excel_region_anno_output(cross = "3051", gr.rb.3051.nacl.t30[[2]])
excel_region_anno_output(cross = "3051", gr.ie.3051.nacl.t30[[2]])



#make bed tracks
dd = data.frame(chr = gr.ie.3051.nacl.t30[[1]]$seqnames, start = gr.ie.3051.nacl.t30[[1]]$start, end=gr.ie.3051.nacl.t30[[1]]$end)
dd$start = as.numeric(dd$start)
dd$end = as.numeric(dd$end)
write.table(dd, file="3051.nacl.t30.ie.aucell.CIs.bed", quote=F, sep="\t", row.names=F, col.names=T)

dd = data.frame(chr = gr.rp.3051.nacl.t30[[1]]$seqnames, start = gr.rp.3051.nacl.t30[[1]]$start, end=gr.rp.3051.nacl.t30[[1]]$end)
dd$start = as.numeric(dd$start)
dd$end = as.numeric(dd$end)
write.table(dd, file="3051.nacl.t30.rp.aucell.CIs.bed", quote=F, sep="\t", row.names=F, col.names=T)

dd = data.frame(chr = gr.rb.3051.nacl.t30[[1]]$seqnames, start = gr.rb.3051.nacl.t30[[1]]$start, end=gr.rb.3051.nacl.t30[[1]]$end)
dd$start = as.numeric(dd$start)
dd$end = as.numeric(dd$end)
write.table(dd, file="3051.nacl.t30.rb.aucell.CIs.bed", quote=F, sep="\t", row.names=F, col.names=T)


##--------3051 SP t0


combined_dir = "/Users/noahalexander/repeat_fine_mapping/combined/"
cross="A"
perturbation = "sp"
timepoint="t0"

#make index granges of variant positions and extract/generate granges of state qtl confidence intervals 
gr.3051.t0.index = make_index_gr(combined_dir, cross, perturbation, timepoint)
gr.list.3051.t0 = make_state_granges(combined_dir, cross, perturbation, timepoint)

#make state granges objects for each of the iesr states with expanded confidence intervals 
#(here by 1 marker in each direction unless the existing border of the CI is chr end)
gr.rp.3051.sp.t0 = expand_interval3((gr.list.3051.t0[[1]]), gr.3051.t0.index)
gr.rb.3051.sp.t0 = expand_interval3(gr.list.3051.t0[[2]], gr.3051.t0.index)
gr.ie.3051.sp.t0 = expand_interval3(gr.list.3051.t0[[3]], gr.3051.t0.index)

#make excel output of subset of vcf for parents and CIs. 
excel_region_anno_output(cross = "3051", gr.rp.3051.sp.t0[[2]])
excel_region_anno_output(cross = "3051", gr.rb.3051.sp.t0[[2]])
excel_region_anno_output(cross = "3051", gr.ie.3051.sp.t0[[2]])



#make bed tracks
dd = data.frame(chr = gr.ie.3051.sp.t0[[1]]$seqnames, start = gr.ie.3051.sp.t0[[1]]$start, end=gr.ie.3051.sp.t0[[1]]$end)
dd$start = as.numeric(dd$start)
dd$end = as.numeric(dd$end)
write.table(dd, file="3051.sp.t0.ie.aucell.CIs.bed", quote=F, sep="\t", row.names=F, col.names=T)

dd = data.frame(chr = gr.rp.3051.sp.t0[[1]]$seqnames, start = gr.rp.3051.sp.t0[[1]]$start, end=gr.rp.3051.sp.t0[[1]]$end)
dd$start = as.numeric(dd$start)
dd$end = as.numeric(dd$end)
write.table(dd, file="3051.sp.t0.rp.aucell.CIs.bed", quote=F, sep="\t", row.names=F, col.names=T)

dd = data.frame(chr = gr.rb.3051.sp.t0[[1]]$seqnames, start = gr.rb.3051.sp.t0[[1]]$start, end=gr.rb.3051.sp.t0[[1]]$end)
dd$start = as.numeric(dd$start)
dd$end = as.numeric(dd$end)
write.table(dd, file="3051.sp.t0.rb.aucell.CIs.bed", quote=F, sep="\t", row.names=F, col.names=T)



##--------3051 SP t10

combined_dir = "/Users/noahalexander/repeat_fine_mapping/combined/"
cross="A"
perturbation = "sp"
timepoint="t10"

#make index granges of variant positions and extract/generate granges of state qtl confidence intervals 
gr.3051.t10.index = make_index_gr(combined_dir, cross, perturbation, timepoint)
gr.list.3051.t10 = make_state_granges(combined_dir, cross, perturbation, timepoint)

#make state granges objects for each of the iesr states with expanded confidence intervals 
#(here by 1 marker in each direction unless the existing border of the CI is chr end)
gr.rp.3051.sp.t10 = expand_interval3((gr.list.3051.t10[[1]]), gr.3051.t10.index)
gr.rb.3051.sp.t10 = expand_interval3(gr.list.3051.t10[[2]], gr.3051.t10.index)
gr.ie.3051.sp.t10 = expand_interval3(gr.list.3051.t10[[3]], gr.3051.t10.index)

#make excel output of subset of vcf for parents and CIs. 
excel_region_anno_output(cross = "3051", gr.rp.3051.sp.t10[[2]])
excel_region_anno_output(cross = "3051", gr.rb.3051.sp.t10[[2]])
excel_region_anno_output(cross = "3051", gr.ie.3051.sp.t10[[2]])




#make bed tracks
dd = data.frame(chr = gr.ie.3051.sp.t10[[1]]$seqnames, start = gr.ie.3051.sp.t10[[1]]$start, end=gr.ie.3051.sp.t10[[1]]$end)
dd$start = as.numeric(dd$start)
dd$end = as.numeric(dd$end)
write.table(dd, file="3051.sp.t10.ie.aucell.CIs.bed", quote=F, sep="\t", row.names=F, col.names=T)

dd = data.frame(chr = gr.rp.3051.sp.t10[[1]]$seqnames, start = gr.rp.3051.sp.t10[[1]]$start, end=gr.rp.3051.sp.t10[[1]]$end)
dd$start = as.numeric(dd$start)
dd$end = as.numeric(dd$end)
write.table(dd, file="3051.sp.t10.rp.aucell.CIs.bed", quote=F, sep="\t", row.names=F, col.names=T)

dd = data.frame(chr = gr.rb.3051.sp.t10[[1]]$seqnames, start = gr.rb.3051.sp.t10[[1]]$start, end=gr.rb.3051.sp.t10[[1]]$end)
dd$start = as.numeric(dd$start)
dd$end = as.numeric(dd$end)
write.table(dd, file="3051.sp.t10.rb.aucell.CIs.bed", quote=F, sep="\t", row.names=F, col.names=T)


