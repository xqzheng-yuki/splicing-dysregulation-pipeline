import_bigwig <- function(file,bw_selection) {
  bw_data <- BRGenomics::makeGRangesBRG(import.bw(file,selection=rtracklayer::BigWigSelection(bw_selection)))
  keepStandardChromosomes(bw_data,pruning.mode = "coarse")
}

import_bigwig_ori <- function(file) {
  bw_data <- import.bw(file)
  keepStandardChromosomes(bw_data,pruning.mode = "coarse")
}

import_bigwig("/mnt/gtklab01/xiaoqing/star/results/group/unmapped_CTX_154_m1.bw",
              bw_selection=gr)
import_bigwig_ori("/mnt/gtklab01/xiaoqing/star/results/group/unmapped_CTX_154_m1.bw")


foo <- control_bw_data[[1]]

sapply(map(foo,function(x) x$score),sum)

show_track<-list()
show_track <- append(show_track,list(idt))
show_track <- append(show_track,list(grTrack))
show_track <- append(show_track,control_datatrack)
show_track <- append(show_track,ko_datatrack)
show_track <- append(show_track,list(gaxis))



