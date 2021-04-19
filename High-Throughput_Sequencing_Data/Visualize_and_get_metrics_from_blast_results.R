
# title: "Methods - Visualize and get integration site metrics from BLAST results"
# author: "Tanner Wiegand"
# date: "4/16/2020"
# Note 1: This script is meant to be run in RStudio (Written and tested in v1.2.5001)
# Note 2: This script processes one set of reads (i.e. one barcode) at a time
# Note 3: p1, p2, p3A, and p3B refer to primers listed in supplemental data of accompanying manuscript


# Initiate required R Packages
library("Biostrings")
library("ggplot2")
library("tidyverse")

# Clear saved variables
rm(list=ls())

# Set working directory (path to folder containing: High-Throughput_Sequencing_Data/)
setwd("High-Throughput_Sequencing_Data/")

# Set barcode you would like to analyze:
BC <- 1
  # WT          = 1
  # -5D         = 2
  # +5D         = 3
  # +10D        = 4
  # +10D, +10U  = 5
  # ∆M          = 6
  # WT∆IHF      = 9

##Import HTS Reads
  # Reads were processed using the methods described in HTS_Paired-End_Reads_Processing_and_BLAST_Parameters.Rmd
seq <- readDNAStringSet(paste0("Fasta_Files/BC", BC ,".masked.fna"))

##Function to quickly check the sequence of a read
zeq <- function(a) {as.character(seq[a])}


###Sort to figure out which primer each read starts from
mm <- 3 # Mismatches that will be tolerated in primer search

#Enter Primer Sequences [Note: primers 3A and 3B will be ignored, since they have so few hits]
p1 <- "ACCACCCGGCTTTCTTAG"
p2 <- "CCAATTGCCCGAAGCTTC"
#p3A <- "TGGTCTAGGAAAAGGACTTGGAC"
#p3B <- "CTCCTTGTCCAAGTCCTTTTCC"

#Get and check subseq of first 40 nucleotides for each read
short.seq <- subseq(seq, start =1, end =40)

##Count number of occurences of primer in each read, convert to true/false and select only those reads
p1.hits <- vcountPattern(p1, short.seq, max.mismatch = mm)
p1.hits <- p1.hits == 1
p1.seqs <- seq[p1.hits]

p2.hits <- vcountPattern(p2, short.seq, max.mismatch = mm)
p2.hits <- p2.hits == 1
p2.seqs <- seq[p2.hits]

### Only 10-100 hits from p3A and p3B, so we will ignore them
#p3A.hits <- vcountPattern(p3A, short.seq, max.mismatch = mm)
#p3A.hits <- p3A.hits == 1
#p3A.seqs <- seq[p3A.hits]

#p3B.hits <- vcountPattern(p3B, short.seq, max.mismatch = mm)
#p3B.hits <- p3B.hits == 1
#p3B.seqs <- seq[p3B.hits]

rm(p1.hits, p2.hits) 


###Import BLAST results
  # BLAST was done according to methods outlined in HTS_Paired-End_Reads_Processing_and_BLAST_Parameters.Rmd
bc.x <- read.delim(paste0("Blast-results/CRISPRs/BC", BC ,".masked.txt"), header = FALSE, sep = "\t")
bc.x.int <- read.delim(paste0("Blast-results/Int_Substrate/BC", BC, ".masked.int.txt"), header = FALSE, sep = "\t")

# Renames columns on blast results
colnames(bc.x) <- c("id.quer", "len.quer", "start.quer", "end.quer", "start.sub", "end.sub", "strand.sub", 
                    "e.value", "len.algn", "pident", "nident", "mmatch", "gaps")
colnames(bc.x.int) <- colnames(bc.x)

## Consolidate blast info for each read
bc.c.pos <- bc.x %>% filter(strand.sub == "plus")
bc.c.neg <- bc.x %>% filter(strand.sub == "minus")

# Summarize read data from p1
bc.c.pos <- bc.c.pos %>% group_by(id.quer) %>% summarize(len.quer = max(len.quer), 
                                                         start.quer = min(start.quer), 
                                                         end.quer = max(end.quer), 
                                                         start.sub = min(start.sub),
                                                         end.sub = max(end.sub), 
                                                         strand.sub = "plus", 
                                                         len.algn = max(len.algn))

# Summarize read data from p2
bc.c.neg <- bc.c.neg %>% group_by(id.quer) %>% summarize(len.quer = max(len.quer), 
                                                         start.quer = min(start.quer), 
                                                         end.quer = max(end.quer), 
                                                         start.sub = max(start.sub),
                                                         end.sub = min(end.sub), 
                                                         strand.sub = "minus", 
                                                         len.algn = max(len.algn))

# Concatenate read data summaries from positivie and negative strands
bc.x <- rbind(bc.c.pos, bc.c.neg)
bc.x <- bc.x[order(bc.x$id.quer),]

# Remove the few reads that made it into both lists
x <- !bc.x$id.quer %in% bc.x$id.quer[duplicated(bc.x$id.quer)]
bc.x <- bc.x[x,]

###Filter BLAST results for reads that occur in both in CRISPR and Integrated substrate subject searches 
bc.x <- bc.x %>%  filter(id.quer %in% bc.x.int$id.quer)
bc.x.int <- bc.x.int %>%  filter(id.quer %in% bc.x$id.quer)


###Filter BLAST results by primer
#p1: CRISPR
x <- bc.x$id.quer %in% as.character(names(p1.seqs))
bc.p1.x <-bc.x[x,] 
bc.p1.x$primer <- "p1"

#p1: Integrated Substrate
x <- bc.x.int$id.quer %in% as.character(names(p1.seqs))
bc.p1.x.int <-bc.x.int[x,] 
bc.p1.x.int$primer <- "p1"

#p2: CRISPR
x <- bc.x$id.quer %in% as.character(names(p2.seqs))
bc.p2.x <-bc.x[x,] 
bc.p2.x$primer <- "p2"

#p2: Integrated Substrate
x <- bc.x.int$id.quer %in% as.character(names(p2.seqs))
bc.p2.x.int <-bc.x.int[x,] 
bc.p2.x.int$primer <- "p2"

###Figure out which reads have ambiguos insertion sites
#p1
bc.p1.x.overlap <- data.frame(bc.p1.x.int$id.quer)

colnames(bc.p1.x.overlap) <- c("id.quer")
temp1 <- bc.p1.x.int %>% select(id.quer, start.quer)
temp2 <- bc.p1.x %>% select(id.quer, end.quer)

bc.p1.x.overlap <- bc.p1.x.overlap %>% left_join(temp1, by = "id.quer")
bc.p1.x.overlap <- bc.p1.x.overlap %>% left_join(temp2, by = "id.quer")

colnames(bc.p1.x.overlap) <- c("read", "start", "end")
bc.p1.x.overlap$overlap.len <- bc.p1.x.overlap$end - bc.p1.x.overlap$start + 1

bc.p1.x.no.overlap <- bc.p1.x.overlap %>% filter(overlap.len < 1)
bc.p1.x.overlap <- bc.p1.x.overlap %>% filter(overlap.len > 0)


#p2
bc.p2.x.overlap <- data.frame(bc.p2.x.int$id.quer)

colnames(bc.p2.x.overlap) <- c("id.quer")
temp1 <- bc.p2.x.int %>% select(id.quer, start.quer)
temp2 <- bc.p2.x %>% select(id.quer, end.quer)

bc.p2.x.overlap <- bc.p2.x.overlap %>% left_join(temp1, by = "id.quer")
bc.p2.x.overlap <- bc.p2.x.overlap %>% left_join(temp2, by = "id.quer")

colnames(bc.p2.x.overlap) <- c("read", "start", "end")
bc.p2.x.overlap$overlap.len <- bc.p2.x.overlap$end - bc.p2.x.overlap$start + 1

bc.p2.x.no.overlap <- bc.p2.x.overlap %>% filter(overlap.len < 1)
bc.p2.x.overlap <- bc.p2.x.overlap %>% filter(overlap.len > 0)


###Check for trimming
#Combine non-ambiguos reads for p1 and p2
temp1 <- bc.p1.x.int[bc.p1.x.int$id.quer %in% bc.p1.x.no.overlap$read, ]
temp2 <- bc.p2.x.int[bc.p2.x.int$id.quer %in% bc.p2.x.no.overlap$read, ]
trimming <- rbind (temp1, temp2)


df.trim <- data.frame(seq(1, 40))

temp1 <- data.frame(table(trimming$start.sub[trimming$strand.sub == "plus"]))
temp2 <- data.frame(table(trimming$start.sub[trimming$strand.sub == "minus"]))
temp1 <- rbind(temp1, temp2) ; rm(temp2)
colnames(temp1) <- c("length", "quant")

colnames(df.trim) <- c("length")
df.trim$quant <- 0
x <- df.trim$length %in% temp1$length

df.trim$quant[x] <- temp1$quant
df.trim$percent <- round(100 * df.trim$quant / sum(df.trim$quant), 2)
df.trim$percent[df.trim$length <20] <- df.trim$percent[df.trim$length <20] *-1

#Graph Trimming Events
tiz <- paste0("Trimming of Spacers with Unambiguos Insertion Sites")
tiz.sub <- paste0("(n = ", format(nrow(trimming), big.mark = ","), " Spacers)")

g1 <- ggplot() + geom_col(data = df.trim, aes(length, percent)) +
  theme_bw() +
labs(title = tiz, subtitle = tiz.sub, y = "Percent of Half-Site Spacers", x = "Position") +
   geom_text(data = df.trim,aes(x = length, y = percent, label = paste0(percent, "%")))

print(g1)

###Graph Insertion Sites for unambiguos integration events
crispr.list <- c("WT", "minus5bpi", "plus5bpi", "plus10bpi", "plus10b-Api", "IRp", "IRd",  "Del-IHFd", "WT", "WT")
genz <- read_delim("Genetic_feature_locations.txt", delim = "\t")
genz <- genz[BC,]
lengz <- genz$S2_e

x <- bc.p2.x$id.quer %in% bc.p2.x.no.overlap$read
temp.x <- bc.p2.x[x,]

df.pos <- data.frame(table(temp.x$end.sub))
df.pos$percent <- round(100 * df.pos$Freq/ sum(df.pos$Freq),2)
df.pos$Var1 <- as.numeric(as.character(df.pos$Var1))

temp1 <- data.frame(seq(1, lengz)) ; colnames(temp1) <- "Var1"
temp1 <- temp1 %>% left_join(df.pos, by = "Var1")
temp1$Freq[is.na(temp1$Freq)] <- 0
temp1$percent[is.na(temp1$percent)] <- 0
df.pos <- temp1

x <- bc.p1.x$id.quer %in% bc.p1.x.no.overlap$read
temp.x <- bc.p1.x[x,]

df.neg <- data.frame(table(temp.x$end.sub))
df.neg$percent <- round(-100 * df.neg$Freq/ sum(df.neg$Freq),2)
df.neg$Var1 <- as.numeric(as.character(df.neg$Var1))

temp1 <- data.frame(seq(1, lengz)) ; colnames(temp1) <- "Var1"
temp1 <- temp1 %>% left_join(df.neg, by = "Var1")
temp1$Freq[is.na(temp1$Freq)] <- 0
temp1$percent[is.na(temp1$percent)] <- 0
df.neg <- temp1


##Graph
g2 <- ggplot() + 
  geom_segment(aes(x = 1, y = 0, xend = lengz, yend = 0)) +
  geom_segment(aes(x = n5_s, y = 0, xend = n5_e, yend = 0, colour = "n"), alpha = 0.6, size = 2, data = genz) +
  geom_segment(aes(x = IRd_s, y = 0, xend = IRd_e, yend = 0, colour = "IR"), alpha = 0.6, size = 2, data = genz) +
  geom_segment(aes(x = n4_s, y = 0, xend = n4_e, yend = 0, colour = "n"), alpha = 0.6, size = 2, data = genz) +
  geom_segment(aes(x = IHFd_s, y = 0, xend = IHFd_e, yend = 0, colour = "IHF"), alpha = 0.6, size = 2, data = genz) +
  geom_segment(aes(x = n3_s, y = 0, xend = n3_e, yend = 0, colour = "n"), alpha = 0.6, size = 2, data = genz) +
  geom_segment(aes(x = IRp_s, y = 0, xend = IRp_e, yend = 0, colour = "IR"), alpha = 0.6, size = 2, data = genz) +
  geom_segment(aes(x = n2_s, y = 0, xend = n2_e, yend = 0, colour = "n"), alpha = 0.6, size = 2, data = genz) +
  geom_segment(aes(x = IHFp_s, y = 0, xend = IHFp_e, yend = 0, colour = "IHF"), alpha = 0.6, size = 2, data = genz) +
  geom_segment(aes(x = n1_s, y = 0, xend = n1_e, yend = 0, colour = "n"), alpha = 0.6, size = 2, data = genz) +
  geom_segment(aes(x = R1_s, y = 0, xend = R1_e, yend = 0, colour = "R"), alpha = 0.6, size = 2, data = genz) +
  geom_segment(aes(x = S1_s, y = 0, xend = S1_e, yend = 0, colour = "S1"), alpha = 0.6, size = 2, data = genz) +
  geom_segment(aes(x = R2_s, y = 0, xend = R2_e, yend = 0, colour = "R"), alpha = 0.6, size = 2, data = genz) +
  geom_segment(aes(x = S2_s, y = 0, xend = S2_e, yend = 0, colour = "S2"), alpha = 0.6, size = 2, data = genz) +
  geom_segment(data = df.pos, aes(x = Var1, y = 0, xend= Var1, yend = -1*percent), color = "navy", alpha = .8) +
  geom_segment(data = df.neg, aes(x = Var1, y = 0, xend= Var1, yend = -1*percent), color = "navy", alpha = .8) +
  labs(title = paste0("Unambigous Integrations into ", crispr.list[BC], " CRISPRs"), y = "Percent of Integration Evens", x = "Position")

print(g2)

########################### Now examine reads with Ambiguos instertion sites
###Working with ambiguos reads
#minus
temp <- bc.p1.x.overlap

temp <- temp[order(temp$read), ]
x <- bc.p1.x$id.quer %in% temp$read
temp$overlap.start <- bc.p1.x$end.sub[x]
temp$overlap.end <- bc.p1.x$end.sub[x]+temp$overlap.len


temp$s.e <- paste0(temp$overlap.start, "-", temp$overlap.end)
df.neg.amb <- data.frame(table(temp$s.e))

df.neg.amb$start <- as.numeric(as.character(word(df.neg.amb$Var1, start = 1, end =1, sep = "-")))
df.neg.amb$end <- as.numeric(as.character(word(df.neg.amb$Var1, start = 2, end =2, sep = "-")))

#plus
temp <- bc.p2.x.overlap

temp <- temp[order(temp$read), ]
x <- bc.p2.x$id.quer %in% temp$read
temp$overlap.start <- bc.p2.x$end.sub[x]-temp$overlap.len
temp$overlap.end <- bc.p2.x$end.sub[x]


temp$s.e <- paste0(temp$overlap.start, "-", temp$overlap.end)
df.pos.amb <- data.frame(table(temp$s.e))

df.pos.amb$start <- as.numeric(as.character(word(df.pos.amb$Var1, start = 1, end =1, sep = "-")))
df.pos.amb$end <- as.numeric(as.character(word(df.pos.amb$Var1, start = 2, end =2, sep = "-")))

#Clean up
df.neg.amb$percent <- round(-100 * df.neg.amb$Freq / (nrow(bc.p1.x) + nrow(bc.p2.x)) , 2)
df.pos.amb$percent <- round(100 * df.pos.amb$Freq / (nrow(bc.p1.x) + nrow(bc.p2.x)) , 2)

df.neg.cer <- df.neg ; df.pos.cer <- df.pos
df.neg.cer$percent <- round(-100 * df.neg.cer$Freq / (nrow(bc.p1.x) + nrow(bc.p2.x)) , 2)
df.pos.cer$percent <- round(100 * df.pos.cer$Freq / (nrow(bc.p1.x) + nrow(bc.p2.x)) , 2)


df.pos.amb$mid <- (df.pos.amb$end + df.pos.amb$start ) /2
df.pos.amb$wid <- df.pos.amb$end - df.pos.amb$start + 1

df.neg.amb$mid <- (df.neg.amb$end + df.neg.amb$start ) /2
df.neg.amb$wid <- df.neg.amb$end - df.neg.amb$start + 1


#Graph (Note: Alignment to minus is actually p1, meaning insertion came from a nucleophilic attack on positive strand)
g3 <- ggplot() + 
  ylim(-10, 40) +
  theme_bw() +
  geom_segment(aes(x = 1, y = 0, xend = lengz, yend = 0)) +
  geom_segment(aes(x = n5_s, y = 0, xend = n5_e, yend = 0, colour = "n"), alpha = 0.6, size = 2, data = genz) +
  geom_segment(aes(x = IRd_s, y = 0, xend = IRd_e, yend = 0, colour = "IR"), alpha = 0.6, size = 2, data = genz) +
  geom_segment(aes(x = n4_s, y = 0, xend = n4_e, yend = 0, colour = "n"), alpha = 0.6, size = 2, data = genz) +
  geom_segment(aes(x = IHFd_s, y = 0, xend = IHFd_e, yend = 0, colour = "IHF"), alpha = 0.6, size = 2, data = genz) +
  geom_segment(aes(x = n3_s, y = 0, xend = n3_e, yend = 0, colour = "n"), alpha = 0.6, size = 2, data = genz) +
  geom_segment(aes(x = IRp_s, y = 0, xend = IRp_e, yend = 0, colour = "IR"), alpha = 0.6, size = 2, data = genz) +
  geom_segment(aes(x = n2_s, y = 0, xend = n2_e, yend = 0, colour = "n"), alpha = 0.6, size = 2, data = genz) +
  geom_segment(aes(x = IHFp_s, y = 0, xend = IHFp_e, yend = 0, colour = "IHF"), alpha = 0.6, size = 2, data = genz) +
  geom_segment(aes(x = n1_s, y = 0, xend = n1_e, yend = 0, colour = "n"), alpha = 0.6, size = 2, data = genz) +
  geom_segment(aes(x = R1_s, y = 0, xend = R1_e, yend = 0, colour = "R"), alpha = 0.6, size = 2, data = genz) +
  geom_segment(aes(x = S1_s, y = 0, xend = S1_e, yend = 0, colour = "S1"), alpha = 0.6, size = 2, data = genz) +
  geom_segment(aes(x = R2_s, y = 0, xend = R2_e, yend = 0, colour = "R"), alpha = 0.6, size = 2, data = genz) +
  geom_segment(aes(x = S2_s, y = 0, xend = S2_e, yend = 0, colour = "S2"), alpha = 0.6, size = 2, data = genz) +
  
  geom_col(dat = df.pos.amb, aes(x = mid, y = -1 * percent), fill = "salmon", alpha = .4, width = df.pos.amb$wid+1, position = "identity") +
  geom_col(dat = df.neg.amb, aes(x = mid, y = -1 * percent), fill = "orange", alpha = .4, width = df.neg.amb$wid+1, position = "identity") +
  
  geom_segment(data = df.pos.cer, aes(x = Var1, y = 0, xend= Var1, yend = -1 * percent), color = "navy") +
  geom_segment(data = df.neg.cer, aes(x = Var1, y = 0, xend= Var1, yend = -1 * percent), color = "midnight blue") +
  
  
  labs(title = paste0("Integrations into ", crispr.list[BC], " CRISPRs"), y = "Percent of Integration Evens", x = "Position")

print(g3)

#########################  Table with Metrics #########################
### reads from p1 ###
##First boundary (L-R1: preferred in polarized acquisition)
#certain
df.metric <- data.frame(-1*df.neg.cer$percent[df.neg.cer$Var1 == genz$R1_s], "Leader-Repeat1", "plus", "certain")
colnames(df.metric) <- c("percent", "location", "strand.attacked", "insertion.class")
df.metric$location <- as.character(df.metric$location)
df.metric$strand.attacked <- as.character(df.metric$strand.attacked)
df.metric$insertion.class <- as.character(df.metric$insertion.class)

#ambiguos
x <- df.neg.amb %>%  filter(start <= genz$R1_s, end >= genz$R1_s)
df.metric[2,] <- c(sum(-1 * x$percent), "Leader-Repeat1", "plus", "ambiguos")
df.metric$percent <- as.numeric(as.character(df.metric$percent))

#combine
df.metric[3,] <- c(sum(df.metric$percent[1:2]), "Leader-Repeat1", "plus", "TOTAL")

##Second boundary (S1-R2)
#certain
df.metric[4,] <- c(-1*df.neg.cer$percent[df.neg.cer$Var1 == genz$R2_s], "Spacer1-Repeat2", "plus", "certain")

#ambiguos
x <- df.neg.amb %>%  filter(start <= genz$R2_s, end >= genz$R2_s)
df.metric[5,] <- c(sum(-1 * x$percent), "Spacer1-Repeat2", "plus", "ambiguos")
df.metric$percent <- as.numeric(as.character(df.metric$percent))

#combine
df.metric[6,] <- c(sum(df.metric$percent[4:5]), "Spacer1-Repeat2", "plus", "TOTAL")



### reads from p2 ###
##First boundary (R1-S1: preferred in polarized acquisition)
#certain
df.metric[7,] <- c(df.pos.cer$percent[df.pos.cer$Var1 == genz$R1_e], "Repeat1-Spacer1", "minus", "certain")

#ambiguos
x <- df.pos.amb %>%  filter(start <= genz$R1_e, end >= genz$R1_e)
df.metric[8,] <- c(sum(x$percent), "Repeat1-Spacer1", "minus", "ambiguos")
df.metric$percent <- as.numeric(as.character(df.metric$percent))

#combine
df.metric[9,] <- c(sum(df.metric$percent[7:8]), "Repeat1-Spacer1", "minus", "TOTAL")

##Second boundary (R2-S2)
#certain
df.metric[10,] <- c(df.pos.cer$percent[df.pos.cer$Var1 == genz$R2_e], "Repeat2-Spacer2", "minus", "certain")

#ambiguos
x <- df.pos.amb %>%  filter(start <= genz$R2_e, end >= genz$R2_e)
df.metric[11,] <- c(sum(x$percent), "Repeat2-Spacer2", "minus", "ambiguos")
df.metric$percent <- as.numeric(as.character(df.metric$percent))

#combine
df.metric[12,] <- c(sum(df.metric$percent[10:11]), "Repeat2-Spacer2", "minus", "TOTAL")
df.metric$percent <- as.numeric(as.character(df.metric$percent))


###Off-site if we assume ambiguos were correct ( == confirmed off-sites)
x <- df.metric %>% filter(insertion.class == "TOTAL")
df.metric[13,] <- c(round(100-sum(x$percent),2), "off-site", "either", "TOTAL")

#order table 
df.metric <- df.metric[order(df.metric$insertion.class),]

df.metric


# Create directory for table output
dir.create("Tables")


#Write out table
write.table(df.metric, file = paste0("Tables/BC", BC, "_integration_locations.txt"), 
            row.names = FALSE, 
            quote = FALSE, 
            sep = "\t")



##Metrics
##Count of reads collected from each primer
x <- length(p1.seqs) + length(p2.seqs) #+ length(p3A.seqs) + length(p3B.seqs)
print(paste0("Primers matched to ", format(x, big.mark = ","), " of ", format(length(seq), big.mark = ","), 
       " total reads [", round(100*x/length(seq),2), "%]"))
print(paste0("p1 = ", format(length(p1.seqs), big.mark = ",") ))
print(paste0("p2 = ", format(length(p2.seqs), big.mark = ",") ))



#Percent of total reads that came from ambiguos vs. certain integrations
print(paste0("Ambiguos integration events = ", format(nrow(bc.p1.x.overlap) + nrow(bc.p2.x.overlap), big.mark = ","),
       " [", round(100*(nrow(bc.p1.x.overlap) + nrow(bc.p2.x.overlap)) / nrow(bc.x),2), "% of mapped integrations]" ))
print(paste0("Certin integration events = ", format(nrow(bc.p1.x.no.overlap) + nrow(bc.p2.x.no.overlap), big.mark = ","),
       " [", round(100*(nrow(bc.p1.x.no.overlap) + nrow(bc.p2.x.no.overlap)) / nrow(bc.x),2), "% of mapped integrations]" ))



#Number of total reads that mapped to both CRISPR and integrated substrate
print(paste0("Number of reads with Blast result for CRISPR and integrated substrate = ", 
       format(nrow(bc.x), big.mark = ",") , " [" , round(100*nrow(bc.x)/length(seq), 2) , "% of total reads]"))


