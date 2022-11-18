

len_list = function(file_name, out_name, ex_method, lib_method){
    a = read.table(file = file_name, header = TRUE, 
                   sep = '\t', stringsAsFactors = FALSE)
    b = a
    b[,2] = b[,2]/sum(b[,2])
    b = cbind(b, rep(ex_method, nrow(b)),
              rep(lib_method, nrow(b)),
              rep(paste(out_name, ex_method, sep = '_'), nrow(b)),
              rep(paste(out_name, lib_method, sep = '_'), nrow(b)))
    colnames(b) = c('len', 'freq', 'extract','lib','sample_ex', 'sample_lib')
    return(b)
}

p = rbind(len_list('len_BG21_skin_PCI_BEMT_collapse_20', 'BG21_skin', 'PCI', 'BEMT'),
          len_list('len_BG21_skin_PCI_BEST_collapse_20', 'BG21_skin', 'PCI', 'BEST'),
          len_list('len_GMNH09_bone_SIS_BEMT_collapse_20', 'GMNH09_bone', 'SIS', 'BEMT'),
          len_list('len_GMNH13_powder_SIS_BEMT_collapse_20', 'GMNH13_powder', 'SIS', 'BEMT'),
          len_list('len_RUSA02_bone_SIS_BEST_collapse_20', 'RUSA02_bone', 'SIS', 'BEST'),
          len_list('len_Y15_hair_PCI_BEMT_collapse_20', 'Y15_hair', 'PCI', 'BEMT'),
          len_list('len_GMNH13_powder_SIS_BEST_collapse_20', 'GMNH13_powder', 'SIS', 'BEST'),
          len_list('len_RUSA08_bone_DAB_BEMT_collapse_20', 'RUSA08_bone', 'DAB', 'BEMT'),
          len_list('len_Y15_skin_PCI_BEST_collapse_20', 'Y15_skin', 'PCI', 'BEST'),
          len_list('len_NGH1_skin_PCI_BEMT_collapse_20', 'NGH1_skin', 'PCI', 'BEMT'),
          len_list('len_RUSA08_bone_DAB_BEST_collapse_20', 'RUSA08_bone', 'DAB', 'BEST'),
          len_list('len_Y3_skin_PCI_BEMT_collapse_20', 'Y3_skin', 'PCI', 'BEMT'),
          len_list('len_GMNH09_china_DAB_BEMT_collapse_20', 'GMNH09_china', 'DAB', 'BEMT'),
          len_list('len_NGH1_skin_PCI_BEST_collapse_20', 'NGH1_skin', 'PCI', 'BEST'),
          len_list('len_RUSA08_bone_SIS_BEMT_collapse_20', 'RUSA08_bone', 'SIS', 'BEMT'),
          len_list('len_Y3_skin_PCI_BEST_collapse_20', 'Y3_skin', 'PCI', 'BEST'),
          len_list('len_GMNH09_china_SIS_BEMT_collapse_20', 'GMNH09_china', 'SIS', 'BEMT'),
          len_list('len_RUSA02_bone_DAB_BEMT_collapse_20', 'RUSA02_bone', 'DAB', 'BEMT'),
          len_list('len_RUSA08_bone_SIS_BEST_collapse_20', 'RUSA08_bone', 'SIS', 'BEST'),
          len_list('len_GMNH13_powder_DAB_BEMT_collapse_20', 'GMNH13_powder', 'DAB', 'BEMT'),
          len_list('len_RUSA02_bone_DAB_BEST_collapse_20', 'RUSA02_bone', 'DAB', 'BEST'),
          len_list('len_RUSA08_soil_DAB_BEST_collapse_20', 'RUSA08_soil', 'DAB', 'BEST'),
          len_list('len_GMNH13_powder_DAB_BEST_collapse_20', 'GMNH13_powder', 'DAB', 'BEST'),
          len_list('len_RUSA02_bone_SIS_BEMT_collapse_20', 'RUSA02_bone', 'SIS', 'BEMT'),
          len_list('len_RUSA08_soil_SIS_BEMT_collapse_20', 'RUSA08_soil', 'SIS', 'BEMT')
)

library(ggplot2)
library(Cairo)

CairoPDF(file = 'plot_comp_extract_method', height = 9, width = 12)
ggplot(as.data.frame(p)) + geom_line(aes(x=len, y=freq, color=extract), size = 0.5) + theme_bw() +
    facet_wrap(~sample_lib) + scale_colour_manual(values=c('deepskyblue', 'Brown', 'Green3')) +
    xlab('Sequence Length') + ylab('Length Frequency')

dev.off()

CairoPDF(file = 'plot_comp_lib_method', height = 9, width = 12)
ggplot(as.data.frame(p)) + geom_line(aes(x=len, y=freq, color=lib)) + theme_bw() +
    facet_wrap(~sample_ex) + scale_colour_manual(values=c('deepskyblue', 'Brown', 'Green3')) +
    xlab('Sequence Length') + ylab('Length Frequency')

dev.off()


### GC ###
g = rbind(len_list('GC_BG21_skin_PCI_BEMT', 'BG21_skin', 'PCI', 'BEMT'),
          len_list('GC_BG21_skin_PCI_BEST', 'BG21_skin', 'PCI', 'BEST'),
          len_list('GC_GMNH09_bone_SIS_BEMT', 'GMNH09_bone', 'SIS', 'BEMT'),
          len_list('GC_GMNH13_powder_SIS_BEMT', 'GMNH13_powder', 'SIS', 'BEMT'),
          len_list('GC_RUSA02_bone_SIS_BEST', 'RUSA02_bone', 'SIS', 'BEST'),
          len_list('GC_Y15_hair_PCI_BEMT', 'Y15_hair', 'PCI', 'BEMT'),
          len_list('GC_GMNH13_powder_SIS_BEST', 'GMNH13_powder', 'SIS', 'BEST'),
          len_list('GC_RUSA08_bone_DAB_BEMT', 'RUSA08_bone', 'DAB', 'BEMT'),
          len_list('GC_Y15_skin_PCI_BEST', 'Y15_skin', 'PCI', 'BEST'),
          len_list('GC_NGH1_skin_PCI_BEMT', 'NGH1_skin', 'PCI', 'BEMT'),
          len_list('GC_RUSA08_bone_DAB_BEST', 'RUSA08_bone', 'DAB', 'BEST'),
          len_list('GC_Y3_skin_PCI_BEMT', 'Y3_skin', 'PCI', 'BEMT'),
          len_list('GC_GMNH09_china_DAB_BEMT', 'GMNH09_china', 'DAB', 'BEMT'),
          len_list('GC_NGH1_skin_PCI_BEST', 'NGH1_skin', 'PCI', 'BEST'),
          len_list('GC_RUSA08_bone_SIS_BEMT', 'RUSA08_bone', 'SIS', 'BEMT'),
          len_list('GC_Y3_skin_PCI_BEST', 'Y3_skin', 'PCI', 'BEST'),
          len_list('GC_GMNH09_china_SIS_BEMT', 'GMNH09_china', 'SIS', 'BEMT'),
          len_list('GC_RUSA02_bone_DAB_BEMT', 'RUSA02_bone', 'DAB', 'BEMT'),
          len_list('GC_RUSA08_bone_SIS_BEST', 'RUSA08_bone', 'SIS', 'BEST'),
          len_list('GC_GMNH13_powder_DAB_BEMT', 'GMNH13_powder', 'DAB', 'BEMT'),
          len_list('GC_RUSA02_bone_DAB_BEST', 'RUSA02_bone', 'DAB', 'BEST'),
          len_list('GC_RUSA08_soil_DAB_BEST', 'RUSA08_soil', 'DAB', 'BEST'),
          len_list('GC_GMNH13_powder_DAB_BEST', 'GMNH13_powder', 'DAB', 'BEST'),
          len_list('GC_RUSA02_bone_SIS_BEMT', 'RUSA02_bone', 'SIS', 'BEMT'),
          len_list('GC_RUSA08_soil_SIS_BEMT', 'RUSA08_soil', 'SIS', 'BEMT')
)

library(ggplot2)
library(Cairo)

CairoPDF(file = 'plot_GC_comp_extract_method', height = 9, width = 12)
ggplot(as.data.frame(g)) + geom_line(aes(x=len, y=freq, color=extract), size = 0.5) + theme_bw() +
    facet_wrap(~sample_lib) + scale_colour_manual(values=c('deepskyblue', 'Brown', 'Green3')) +
    xlab('GC Content') + ylab('GC Content Frequency')

dev.off()

CairoPDF(file = 'plot_GC_comp_lib_method', height = 9, width = 12)
ggplot(as.data.frame(g)) + geom_line(aes(x=len, y=freq, color=lib)) + theme_bw() +
    facet_wrap(~sample_ex) + scale_colour_manual(values=c('deepskyblue', 'Brown', 'Green3')) +
    xlab('GC Content') + ylab('GC Content Frequency')

dev.off()

### get len ###

len_get = function(file_name, out_name, ex_method, lib_method){
    a = read.table(file = file_name, header = FALSE, 
                   sep = ' ', stringsAsFactors = FALSE)
    b = a
    b[,2] = b[,2]/sum(b[,2])
    b_len = sum(b[,1] * b[,2])
    b = c(b_len, ex_method, lib_method, out_name,
          paste(out_name, ex_method, sep = '_'),
          paste(out_name, lib_method, sep = '_')
          )
    return(b)
}

l = rbind(len_get('len_BG21_skin_PCI_BEMT_collapse_20', 'BG21_skin', 'PCI', 'BEMT'),
          len_get('len_BG21_skin_PCI_BEST_collapse_20', 'BG21_skin', 'PCI', 'BEST'),
          len_get('len_GMNH09_bone_SIS_BEMT_collapse_20', 'GMNH09_bone', 'SIS', 'BEMT'),
          len_get('len_GMNH13_powder_SIS_BEMT_collapse_20', 'GMNH13_powder', 'SIS', 'BEMT'),
          len_get('len_RUSA02_bone_SIS_BEST_collapse_20', 'RUSA02_bone', 'SIS', 'BEST'),
          len_get('len_Y15_hair_PCI_BEMT_collapse_20', 'Y15_hair', 'PCI', 'BEMT'),
          len_get('len_GMNH13_powder_SIS_BEST_collapse_20', 'GMNH13_powder', 'SIS', 'BEST'),
          len_get('len_RUSA08_bone_DAB_BEMT_collapse_20', 'RUSA08_bone', 'DAB', 'BEMT'),
          len_get('len_Y15_skin_PCI_BEST_collapse_20', 'Y15_skin', 'PCI', 'BEST'),
          len_get('len_NGH1_skin_PCI_BEMT_collapse_20', 'NGH1_skin', 'PCI', 'BEMT'),
          len_get('len_RUSA08_bone_DAB_BEST_collapse_20', 'RUSA08_bone', 'DAB', 'BEST'),
          len_get('len_Y3_skin_PCI_BEMT_collapse_20', 'Y3_skin', 'PCI', 'BEMT'),
          len_get('len_GMNH09_china_DAB_BEMT_collapse_20', 'GMNH09_china', 'DAB', 'BEMT'),
          len_get('len_NGH1_skin_PCI_BEST_collapse_20', 'NGH1_skin', 'PCI', 'BEST'),
          len_get('len_RUSA08_bone_SIS_BEMT_collapse_20', 'RUSA08_bone', 'SIS', 'BEMT'),
          len_get('len_Y3_skin_PCI_BEST_collapse_20', 'Y3_skin', 'PCI', 'BEST'),
          len_get('len_GMNH09_china_SIS_BEMT_collapse_20', 'GMNH09_china', 'SIS', 'BEMT'),
          len_get('len_RUSA02_bone_DAB_BEMT_collapse_20', 'RUSA02_bone', 'DAB', 'BEMT'),
          len_get('len_RUSA08_bone_SIS_BEST_collapse_20', 'RUSA08_bone', 'SIS', 'BEST'),
          len_get('len_GMNH13_powder_DAB_BEMT_collapse_20', 'GMNH13_powder', 'DAB', 'BEMT'),
          len_get('len_RUSA02_bone_DAB_BEST_collapse_20', 'RUSA02_bone', 'DAB', 'BEST'),
          len_get('len_RUSA08_soil_DAB_BEST_collapse_20', 'RUSA08_soil', 'DAB', 'BEST'),
          len_get('len_GMNH13_powder_DAB_BEST_collapse_20', 'GMNH13_powder', 'DAB', 'BEST'),
          len_get('len_RUSA02_bone_SIS_BEMT_collapse_20', 'RUSA02_bone', 'SIS', 'BEMT'),
          len_get('len_RUSA08_soil_SIS_BEMT_collapse_20', 'RUSA08_soil', 'SIS', 'BEMT')
)

colnames(l) = c('len', 'extract','lib','sample', 'sample_ex', 'sample_lib')


write.table(l, file='len_trim', quote = F, sep = '\t')