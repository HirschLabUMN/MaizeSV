

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} 

require(data.table)
require(dplyr)
require(tidyr)
require(parallel)
require(stringr)

outfile = args[4]
counts = read.table(args[1], head=T)
splits = read.table(args[2])
fsplits = read.table(args[3])

oldNames = c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10") 
newNames = c("exon","chrom","pos.exon","end.exon","gene","pos.gene","end.gene","parent","prog","famnum")

setnames(splits, old = oldNames, new = newNames)
setnames(fsplits, old = c(oldNames,"V11"), new = c(newNames, "source"))
setnames(counts, old = c("Genes"), new = c("exon"))

splits['source']="real"

head(counts)
head(splits)
head(fsplits)



vR = function(df){
  Theta = df$disp[1]
  aa = glmer(value ~ (1|tissue) + (1|sample) + (1 | gene) + (rep | sample), data=df, family=MASS::negative.binomial(theta=Theta))
  bb = as.data.frame(VarCorr(aa))
  ratio = bb[bb$grp=="gene",4] / sum(bb[,4])
  return(ratio)
}

varRatio2 = function(i){
  Prog = as.character(unique(pp$prog)[i])
  dat = pp %>% filter(as.character(prog)==Prog) %>% as.data.frame()
  dat$tissue = as.factor(dat$tissue)
  dat$sample = as.factor(dat$sample)
  dat$rep = as.factor(dat$rep)
  ratio = tryCatch({vR(dat)}, error = function(e) {return(-9)})
  #ratio = vR(dat)
  return(c(ratio, Prog))
}

munge = function(counts, realSplits, fakeSplits){
  splits = rbind(realSplits, fakeSplits)
  d = merge(splits, counts, by = c("exon"), all.y=TRUE)
  d.m = melt(d, id.vars = c("exon","chrom" ,"pos.exon","end.exon","gene","pos.gene","end.gene","parent","prog","famnum", "source"))
  d.n = d.m %>% mutate(ref = str_split_fixed(as.character(variable), "_", 3)[,3], sample = str_split_fixed(as.character(variable), "[.]", 3)[,1], tissue = str_split_fixed(as.character(variable), "[.]", 3)[,2], rep = str_split_fixed(as.character(variable), "[.]", 3)[,3]) %>% mutate(rep = str_split_fixed(as.character(rep),"_",3)[,1]) %>% as.data.frame()
  return(d.n)
}

p = munge(counts, splits, fsplits)

nn = colnames(counts)[2:61]
mm = c(rep("B",20), rep("P",20), rep("W",20))
tt = c(rep(c("A","A","Em","Em","En","En","IE","IE","I","I","L10","L10","L","L","R","R","SC","SC","T","T"),3))
rr = c(rep(c("R1","R2"),30))
coldata = data.frame(row.names = nn, geno = mm, tissue = tt, rep = rr)

nsg = dcast(p, gene + prog ~ variable, value.var = "value", fun.aggregate = sum)

dds = DESeqDataSetFromMatrix(countData = nsg[,3:62], colData = coldata, design = ~ geno + tissue + rep)

dds = DESeq(dds)

disp = data.frame(gene = nsg$gene, disp = dispersions(dds))

pp = merge(p, disp, by=c("gene"), all.x = TRUE)

pp = pp %>% filter(as.character(parent) != as.character(gene)) %>% as.data.frame()

#ff = mclapply(1:length(unique(pp$prog)), FUN = varRatio2, mc.cores=3)

ff = mclapply(1:10, FUN = varRatio2, mc.cores=3)

ff = as.data.frame(do.call(rbind,ff))

write.table(ff, outfile, col.names = F, quote=F)