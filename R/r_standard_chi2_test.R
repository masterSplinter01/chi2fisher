chi2fisher <- function(inputfilepath, outputfilepath,th){
  library(gdsfmt)
  library(foreach)
  library(parallel)
  library(doParallel)
  #inputfilepath <- "/home/dan/Workspace/R/1000G_gds/1000_genomes_FINGBR.gds"
  #outputfilepath <- "/home/dan/Workspace/R/R_standard_chi2/data"
  f <- openfn.gds(filename = inputfilepath, readonly = T)
 
  phenotype <-as.factor(as.numeric (as.factor(as.character(read.gdsn(index.gdsn(f, "phenotype")) ) )))
  phen_levels_number <- nlevels(phenotype)
  phenotype_length <- length(phenotype)
  phenotype  <- as.numeric(phenotype)
  
  genotype <- read.gdsn(index.gdsn(f,"genotype"), c(1,1), c(-1,-1))
  genotypes_num_cols <- ncol(genotype)
  closefn.gds(f)
  start.time <- Sys.time()

  genotype_d <- matrix (genotype, ncol = genotypes_num_cols)
  genotype_d[genotype_d == 0 | genotype_d == 1] <- 0
  genotype_d[genotype_d == 2] <- 1
  
  genotype_r <- matrix(genotype, ncol = genotypes_num_cols)
  genotype_r[genotype_r == 1 | genotype_r == 2] <- 1
  
  genotype_allele1 <- matrix( as.numeric(genotype == 2) + 3 * as.numeric(genotype == 3), nrow = phenotype_length)
  genotype_allele2 <-matrix( as.numeric(genotype == 2 | genotype == 1) + 3 * as.numeric(genotype == 3), nrow = phenotype_length)
  genotype_allele  <-  rbind(genotype_allele2, genotype_allele1)
  phenotype_allele <- c(phenotype, phenotype)
  
  cores <- th #detectCores() 
  cl <- makePSOCKcluster(cores) 
  registerDoParallel(cl)
  
  res <- data.frame(matrix(ncol = genotypes_num_cols, nrow = phenotype_length))
  
  res <- foreach (i =icount( genotypes_num_cols), .combine = rbind) %dopar% { 

    cur_genotype_cd <- genotype[,i]
  
    
    cd_table <- table(phenotype, cur_genotype_cd[cur_genotype_cd != 3])
    
    cd_stat <- chisq.test(cd_table, correct = F)
    stats_cd <-cd_stat$statistic 
    degrees_freedom_cd <- cd_stat$parameter
    p_values_cd <- cd_stat$p.value
    corr_cd <- cor(phenotype, cur_genotype_cd[cur_genotype_cd != 3])
    
    
    cur_genotype_d <- genotype_d[,i]
    
    d_table <- table(phenotype, factor(cur_genotype_d[cur_genotype_d != 3], levels = 0:1))
    expected_table_d <- as.array(margin.table(d_table,1)) %*% t(as.array(margin.table(d_table,2))) / margin.table(d_table)
    min_d <- min(expected_table_d, na.rm = T)
    corr_d <- cor(phenotype, cur_genotype_d[cur_genotype_d != 3])
    
    if (min_d > 5){
      d_stat <- chisq.test(d_table, correct = F)
      stats_d <-format(  d_stat$statistic, scientific = F)
      degrees_freedom_d <- d_stat$parameter
      p_values_d <- d_stat$p.value
    }
    else{
      d_stat <- fisher.test(d_table)
      stats_d <- d_table[1,1]
      p_values_d <- d_stat$p.value
      degrees_freedom_d <- -1
    }  
    
    
    cur_genotype_r <- genotype_r[, i]
    
    r_table <- table(phenotype, factor(cur_genotype_r[cur_genotype_r != 3], levels = 0:1))
    expected_table_r <- as.array(margin.table(r_table,1)) %*% t(as.array(margin.table(r_table,2))) / margin.table(r_table)
    min_r <- min(expected_table_r, na.rm = T)
    corr_r <- cor(phenotype, cur_genotype_r[cur_genotype_r != 3])
    
    if (min_r > 5){
      r_stat <- chisq.test(r_table, correct = F)
      stats_r <-r_stat$statistic
      degrees_freedom_r <- r_stat$parameter
      p_values_r <- r_stat$p.value
    }
    else {
      r_stat <- fisher.test(r_table)
      stats_r <- r_table[1,1] 
      p_values_r <- r_stat$p.value
      degrees_freedom_r <- -1 
    }
  
  
  cur_genotype_allele <- genotype_allele[,i]
  
  allele_table <- table(phenotype_allele, factor(cur_genotype_allele[cur_genotype_allele != 3], levels = 0:1))
  
  corr_allele <- cor(phenotype_allele, cur_genotype_allele[cur_genotype_allele != 3])  
  
  expected_table_allele <- as.array(margin.table(allele_table,1)) %*% t(as.array(margin.table(allele_table,2))) / margin.table(allele_table)
  
  min_allele <- min(expected_table_allele, na.rm = T)
  
  if (min_allele > 5){
    allele_stat <- chisq.test(allele_table, correct = F)
    stats_allele <- allele_stat$statistic
    degrees_freedom_allele <- allele_stat$parameter
    p_values_allele <- allele_stat$p.value
  }
  else {
    allele_stat <- fisher.test(allele_table)
    stats_allele <- allele_table[1,1]
    p_values_allele <- allele_stat$p.value
    degrees_freedom_allele <- -1 
  }
  
  return(c(cd = stats_cd,
                    df_cd = degrees_freedom_cd,
                    p_values_cd = p_values_cd,
                    corr_cd = corr_cd,
                    d = stats_d,
                    df_d = degrees_freedom_d,
                    p_values_d = p_values_d,
                    corr_d = corr_d,
                    r = stats_r,
                    df_r = degrees_freedom_r,
                    p_values_r = p_values_r,
                    corr_r = corr_r,
                    allele = stats_allele,
                    df_allele = degrees_freedom_allele,
                    p_values_allele = p_values_allele,
                    corr_allele = corr_allele ))
  }
  stopImplicitCluster()

  end.time <- Sys.time()
  time.taken <- end.time - start.time
  
  filename <- "chi2fisherOutput.csv"
  dir.create(outputfilepath, showWarnings = F)
  outcsv <- file.path(outputfilepath, filename)
  
  write.csv(res, file = outcsv)
  print (time.taken)
}
  
  
  
