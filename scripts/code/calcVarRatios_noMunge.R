args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} 

outfile = args[2]
pp = read.table(args[1], head=T)

require(lme4)
require(inline)
require(parallel)
require(dplyr)
require(tibble)

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
  print(Prog)
  df = pp %>% filter(as.character(prog)==Prog & as.character(parent) == as.character(gene)) %>% as.data.frame()
  Theta = df$disp[1]
  dat = pp %>% filter(as.character(prog)==Prog & as.character(parent) != as.character(gene)) %>% as.data.frame()
  results = tryCatch({
    getStats(dat)
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

includes <- '#include <sys/wait.h>'
code <- 'int wstat; while (waitpid(-1, &wstat, WNOHANG) > 0) {};'
wait <- cfunction(body=code, includes=includes, convention='.C')

pp = pp %>% filter(nlevels(droplevels(gene)) > 1) %>% as.data.frame()

results = data.frame("Df" = as.numeric(), "AIC" = as.numeric(), "BIC" = as.numeric(), "logLik" = as.numeric(), "deviance" = as.numeric(), "Chisq" = as.numeric(), "Chi Df" = as.numeric(), "Pr(>Chisq)" = as.numeric(), "varRatio" = as.numeric(), "prog" = as.character(), "parent" = as.character(), "status" = as.character(), "source" = as.character())

runOne = function(i){
  Prog = as.character(unique(pp$prog)[i])
  df = pp %>% filter(as.character(prog)==Prog & as.character(parent) == as.character(gene)) %>% as.data.frame()
  print(head(df))
  Theta = df$disp[1]
  dat = pp %>% filter(as.character(prog)==Prog & as.character(parent) != as.character(gene)) %>% as.data.frame()
  print(head(dat))
  aa = glmer(value ~ (1|tissue) + (1|sample) + (1 | gene), data=dat, family=MASS::negative.binomial(theta=Theta))
  return(list(aa, dat))
}

for (i in 1:(length(unique(pp$prog))/10)){
# for (i in 1:10){
  j = i * 10
  progress = round(i/(length(unique(pp$prog))/10), digits = 4)
  print(paste0("Running jobs: ", j - 9, "-", j, "; (", progress, "% complete)"))
  runNB(j)
  ff = mclapply((j - 9):j, FUN = runNB, mc.cores=2)
  wait()
  ff = as.data.frame(do.call(rbind,ff))
  results = rbind(results, ff)
  # print(tail(results))
}

write.table(results, outfile, row.names = F, quote=F)