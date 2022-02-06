#!/usr/bin/env Rscript

main <- function() {
  ## open libraries
  packages = c('ggplot2','stringr')
  # check.packages function: install and load multiple R packages.
  # Check to see if packages are installed. Install them if they are not, then load them into the R session.
  check.packages <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg)) 
      install.packages(new.pkg, dependencies = TRUE, repo="http://cran.rstudio.com/")
    sapply(pkg, require, character.only = TRUE)
  }
  msg.trap <- capture.output( suppressMessages( check.packages(packages)) ) # silent mode
  
  ## parse arguments
  args <- commandArgs(trailingOnly = TRUE)
  ## parse arguments
  file <- args[1]  ## needed for importing data
  # output option (argument #2)
  if (is.na(args[2])==FALSE){
    output_file=args[2]
  } else {
    output_file=paste(str_split(file, "[.]")[[1]][1],'.jpeg',sep="")
  }
  # graph full canonical toxin length (argument #3)
  if (is.na(args[3])==FALSE){
    toxlen = as.numeric(args[3])
  }
  # graph height option (argument #4)
  if (is.na(args[4])==FALSE){
    grheight = args[4]
  } else {
    grheight = 3
  }
  # graph width option (argument #5)
  if (is.na(args[5])==FALSE){
    grwidth = args[5]
  } else {
    grwidth = 5
  }
  
  
  #****************** MAIN SCRIPT *******************

  ## Read data from pre-formatted dataset.
  depth.df <- read.table(file, sep="\t", header = TRUE)
  names(depth.df) = c('sample','pos','depth')
  depth.df$depth[depth.df$depth == 0] <- NA
  depth.df = na.omit(depth.df)
  filename = sapply(str_split(file, "[/]"), tail, 1)
  filename_v = as.vector(unlist(str_split(filename, "_")))
  referencename = filename_v[length(filename_v)-2]
  samplename = paste(filename_v[1:(length(filename_v)-3)], collapse=" ")
  
  ############## AREA PLOT #################
  if (exists("toxlen")){
    ggplot(depth.df, aes(x=pos, y=depth)) +
      scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                    labels = scales::trans_format("log10", scales::math_format(10^.x))) +
      geom_area(fill = "lightblue") +
      theme_classic() +
      labs(x="Position (nt)", y="Coverage", 
           title=paste(samplename,"Reads Mapped to",referencename)) +
      scale_x_continuous(breaks = seq(0, toxlen, by=400), limits=c(0,as.numeric(toxlen)))
    ggsave(filename=output_file, height=as.numeric(grheight),
           width=as.numeric(grwidth))
  } else {
    ggplot(depth.df, aes(x=pos, y=depth)) +
      scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                    labels = scales::trans_format("log10", scales::math_format(10^.x))) +
      geom_area(fill = "lightblue") +
      theme_classic() +
      labs(x="Position (nt)", y="Coverage", 
           title=paste(samplename,"Reads Mapped to",referencename))
    ggsave(filename=output_file, height=as.numeric(grheight),
           width=as.numeric(grwidth))
  }

}
main()
