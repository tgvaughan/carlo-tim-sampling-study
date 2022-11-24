library(ggplot2)
library(gridExtra)

setwd("./figures/")
lf <- list.files("./")
type1 <- read.table(file = lf[1], header = T)
type2 <- read.table(file = lf[2], header = T)

names(type1)

lapply(type1[, 1:5], unique)
lapply(type2[, 1:5], unique)

table(type1[,1:4])
table(type2[,1:4])


make.plot <- function(df, param, truth=NULL) {
  library(ggplot2)
  p<-ggplot(df[df$name == param,], aes(as.factor(m_truth), median, col=as.factor(r))) + 
    geom_point() +
    geom_errorbar(aes(ymin=lower, ymax=upper), width=0.1) +
    geom_hline(yintercept = truth) +
    facet_grid(s1_truth~s2_truth) +
    ggtitle(param)
  print(p)
  return(p)
}


# make.plot(type1, "R0Values", c(1.1, 1.5))
# Row of the facet grid are s in deme1, col are s in deme2
lplots1 <- mapply(make.plot, param=c("origin", "R0Values", "rep", "sampPropValues"), 
                            truth=list(10,    c(1.1, 1.5),  3, c(0.1, 0.2)), 
                                       MoreArgs = list(df=type1),
                 SIMPLIFY = F)       
ggsave("type1_plots.png", plot = grid.arrange(grobs=lplots1), 
       width = 24, height = 15, units = "cm")

# Type2
origin_type2 <- make.plot(type2, "origin", 10)
ggsave("type2_origin_plot.png", plot = origin_type2,
       width = 24, height = 15, units = "cm")

# this is the easier solution I have found to get the function to work
type2$R0.1_truth <- 1.1
type2$R0.2_truth <- 1.5

make.plot2 <- function(df, param) {
  library(ggplot2)
  p<-ggplot(df[df$name == param,], aes(as.factor(m_truth), median, col=as.factor(r))) + 
    geom_point() +
    geom_errorbar(aes(ymin=lower, ymax=upper), width=0.1) +
    geom_point(aes(y=switch (param,
                            m.1 = m_truth,
                            m.2 = m_truth,
                            R0Values.1 = R0.1_truth,
                            R0Values.2 = R0.2_truth,
                            sampPropValues.1 = s1_truth,
                            sampPropValues.2 = s2_truth
    )), shape=4) +
    facet_grid(s1_truth~s2_truth) +
    ggtitle(param)
  print(p)
  return(p)
}

# Row of the facet grid are s in deme1, col are s in deme2
lplots2 <- lapply(c("m.1", "m.2", 
                    "R0Values.1", "R0Values.2", 
                    "sampPropValues.1", "sampPropValues.2"), make.plot2,
                  df=type2) 

ggsave("type2_plots.pdf", plot = marrangeGrob(grobs=lplots2, ncol = 2, nrow = 1), 
       width = 24, height = 15, units = "cm")
