

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
require(DESeq2)
require(lme4)
require(inline)
require(tibble)

outfile = args[4]
outfile2 = args[5]
counts = read.table(args[1], head=T)
splits = read.table(args[2])
fsplits = read.table(args[3])

oldNames = c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10") 
newNames = c("exon","chrom","pos.exon","end.exon","gene","pos.gene","end.gene","parent","prog","famnum")

setnames(splits, old = oldNames, new = newNames)
setnames(fsplits, old = c(oldNames,"V11"), new = c(newNames, "source"))
setnames(counts, old = c("Genes"), new = c("exon"))

splits['source']="real"

error.df <- data.frame("Df" = c(-9, -9), "AIC" = c(-9, -9), "BIC" = c(-9, -9), "logLik" = c(-9, -9), "deviance" = c(-9, -9), "Chisq" = c(-9, -9), "pval" = c(-9, -9), "varRatio" = c(-9, -9), "prog" = c(-9, -9), "parent" = c(-9, -9), "status" = c("error", "error"))

getStats <- function(df, Theta){
  full = glmer(value ~ (1|tissue) + (1|sample) + (1 | gene), data=df, family=MASS::negative.binomial(theta=Theta))
  reduced = glmer(value ~ (1|tissue) + (1|sample), data=df, family=MASS::negative.binomial(theta=Theta))
  cc = as.data.frame(VarCorr(full))
  ratio = cc[cc$grp=="gene",4] / sum(cc[,4])
  dd = as.data.frame(anova(full, reduced))[,-7]
  colnames(dd)[7] = "pval"
  dd['varRatio'] = ratio
  dd['prog'] = df$prog[1]
  dd['parent'] = df$parent[1]
  dd['status'] = "Good"
  dd['source'] = df$source[1]
  return(dd)
}

runNB <- function(i){
  Prog = as.character(unique(pp$prog)[i])
  df = pp %>% filter(as.character(prog)==Prog & as.character(parent) == as.character(gene)) %>% as.data.frame()
  Theta = df$disp[1]
  dat = pp %>% filter(as.character(prog)==Prog & as.character(parent) != as.character(gene)) %>% as.data.frame()
  results = tryCatch({
    getStats(dat, Theta)
  }, warning = function(w) {
    print(w)
    rr = getStats(dat, Theta)
    rr['status'] = "warning"
    return(rr)
  }, error = function(e) {
    errD = cbind(error.df)
    errD$parent = dat$parent[1]
    errD$prog = dat$prog[1]
    errD$source = dat$source[1]
    return(errD)}) 
  return(results)
}

runOne = function(i){
  Prog = as.character(unique(pp$prog)[i])
  df = pp %>% filter(as.character(prog)==Prog & as.character(parent) == as.character(gene)) %>% as.data.frame()
  Theta = df$disp[1]
  dat = pp %>% filter(as.character(prog)==Prog & as.character(parent) != as.character(gene)) %>% as.data.frame()
  aa = glmer(value ~ (1|tissue) + (1|sample) + (1 | gene), data=dat, family=MASS::negative.binomial(theta=Theta))
  return(aa)
}

munge = function(counts, realSplits, fakeSplits){
  splits = rbind(realSplits, fakeSplits)
  d = merge(splits, counts, by = c("exon"), all.y=TRUE)
  d.m = melt(d, id.vars = c("exon","chrom" ,"pos.exon","end.exon","gene","pos.gene","end.gene","parent","prog","famnum", "source"))
  d.n = d.m %>% mutate(ref = str_split_fixed(as.character(variable), "_", 3)[,3], sample = str_split_fixed(as.character(variable), "[.]", 3)[,1], tissue = str_split_fixed(as.character(variable), "[.]", 3)[,2], rep = str_split_fixed(as.character(variable), "[.]", 3)[,3]) %>% mutate(rep = str_split_fixed(as.character(rep),"_",3)[,1]) %>% as.data.frame()
  d.n$tissue = as.factor(d.n$tissue)
  d.n$sample = as.factor(d.n$sample)
  d.n$rep = as.factor(d.n$rep)
  d.n = d.n %>% filter(!exon %in% c("__alignment_not_unique", "__ambiguous", "__no_feature", "__not_aligned", "__too_low_aQual", "__alignment_not_unique")) %>% filter(nlevels(droplevels(gene)) > 1) %>% as.data.frame()
  return(d.n)
}

munge2 = function(fpkm){
  d.m = melt(fpkm, id.vars = c("gene","parent","prog", "famnum", "source","pos.gene","end.gene"))
  d.n = d.m %>% mutate(ref = str_split_fixed(as.character(variable), "_", 3)[,3], sample = str_split_fixed(as.character(variable), "[.]", 3)[,1], tissue = str_split_fixed(as.character(variable), "[.]", 3)[,2], rep = str_split_fixed(as.character(variable), "[.]", 3)[,3]) %>% mutate(rep = str_split_fixed(as.character(rep),"_",3)[,1]) %>% as.data.frame()
  d.n$tissue = as.factor(d.n$tissue)
  d.n$sample = as.factor(d.n$sample)
  d.n$rep = as.factor(d.n$rep)
  d.n$gene = as.factor(d.n$gene)
  d.n = d.n %>% filter(nlevels(droplevels(gene)) > 1 & source != "unchanged") %>% as.data.frame()
  return(d.n)
}

print("Munging...")
p = munge(counts, splits, fsplits)
print("Munging...done")

p['basepairs'] = abs(p$end.exon - p$pos.exon)
gene_lengths = p %>% filter(sample=="B" & tissue=="A" & rep=="R1") %>% group_by(gene) %>% summarize(basepairs = sum(basepairs)) %>% as.data.frame()

nn = colnames(counts)[2:61]
mm = c(rep("B",20), rep("P",20), rep("W",20))
tt = c(rep(c("A","A","Em","Em","En","En","IE","IE","I","I","L10","L10","L","L","R","R","SC","SC","T","T"),3))
rr = c(rep(c("R1","R2"),30))
coldata = data.frame(row.names = nn, geno = mm, tissue = tt, rep = rr)

print("Casting count data...")
nsg = dcast(p, gene + prog ~ variable, value.var = "value", fun.aggregate = sum)
nsg = nsg %>% filter(!is.na(gene) & !duplicated(gene)) %>% as.data.frame()
print("Casting count data...done")

nsg = merge(nsg, gene_lengths, by="gene")

dds = DESeqDataSetFromMatrix(countData = nsg[,3:62], colData = coldata, design = ~ geno + tissue + rep)
rownames(dds) = nsg$gene

mcols(dds) <- cbind(mcols(dds), nsg[63])
dds = DESeq(dds)

print("Extracting DESeq2 data...")
disp = data.frame(gene = nsg$gene, disp = dispersions(dds))
qq = p %>% group_by(gene, parent, prog) %>% summarize(pos.gene = min(pos.exon), end.gene = max(end.exon), famnum=min(famnum), source = min(source)) %>% as.data.frame()
ww = as.data.frame(fpkm(dds, robust=T))
ww = rownames_to_column(ww, var = "gene")
ww = merge(ww, qq, by="gene", all.x=T)
print("Extracting DESeq2 data...done")

print("Re-munging...")
ee = munge2(ww)
print("Re-munging...done")

print("Merging dispersion estimates...")
pp = merge(ee, disp, by=c("gene"), all.x = TRUE)
print("Merging dispersion estimates...done")

print("Filtering for valid progeny sets...")
pp = pp %>% filter(disp != 60) %>% as.data.frame() ## FILTER FOR GENES WITH DISPERSION ==60; MEANS VERY LOW COUNTS
pp = pp %>% group_by(parent, prog) %>% filter(nlevels(droplevels(gene)) > 2) %>% ungroup() %>% as.data.frame()
print("Filtering for valid progeny sets...done")

print("Writing DEseq2 results...")
write.table(pp, outfile2, row.names = F, quote=F)
print("Writing DEseq2 results...done")

rm(ee)
rm(nsg)
rm(ww)
rm(dds)
rm(qq)
rm(p)

includes <- '#include <sys/wait.h>'
code <- 'int wstat; while (waitpid(-1, &wstat, WNOHANG) > 0) {};'
wait <- cfunction(body=code, includes=includes, convention='.C')

results = data.frame("Df" = as.numeric(), "AIC" = as.numeric(), "BIC" = as.numeric(), "logLik" = as.numeric(), "deviance" = as.numeric(), "Chisq" = as.numeric(), "Chi Df" = as.numeric(), "Pr(>Chisq)" = as.numeric(), "varRatio" = as.numeric(), "prog" = as.character(), "parent" = as.character(), "status" = as.character(), "source" = as.character())

for (i in 1:(length(unique(pp$prog))/10)){
  # for (i in 1:10){
  j = i * 10
  print(paste0("Running jobs: ", j - 9, "-", j))
  ff = mclapply((j - 9):j, FUN = runNB, mc.cores=2)
  wait()
  ff = as.data.frame(do.call(rbind,ff))
  results = rbind(results, ff)
}

write.table(results, outfile, row.names = F, quote=F)