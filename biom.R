#!/usr/bin/env Rscript
X = read.table("raw_zOTU2.txt", header = TRUE, sep = "\t")
Y <- read.table("EtsintaxR.txt", header = TRUE, sep = "\t", na.string = "")
table <- merge(X,Y, by = "OTUID")
write.table(table, "OTU_completed_EZtaxon.txt", sep="\t", row.names = F, quote =F)
