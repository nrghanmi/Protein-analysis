
library(ggplot2)
library(writexl)


rm(list=ls())
# ---please set your file path----
setwd("/")

#---- map-mut site ---------
ResultsFrq = read.csv("result_freq.csv")

#------ AlphaMissense site -----
AlphaMissense <- readxl::read_excel(
  path = "AlphaMissense-Search-P42338 3.xls",
)

ResultsFrq$`protein variant` <- paste("p.",ResultsFrq$Wild.Type,ResultsFrq$Residue.Position, ResultsFrq$Mutant.Residue,  sep = "")

# Display the modified data frame

ResultsFrq$ID <- toupper(ResultsFrq$`protein variant`)
AlphaMissense$ID <- toupper(AlphaMissense$`protein variant`)

merged_df <- merge(ResultsFrq, AlphaMissense, by = "ID")
merged_df <- merged_df[,c(1:9,13:16,18:23)]

merged_df$`pathogenicity class` <- factor(merged_df$`pathogenicity class`, levels = c("likely_benign", "likely_pathogenic", "ambiguous"), labels = c("Benign", "Pathogenic", "Ambiguous"))

#write_xlsx(as.data.frame(merged_df),"finalAnalysis.xlsx")


state_counts <- table(merged_df$`pathogenicity class`)

# Convert state counts to a data frame
data <- as.data.frame(state_counts)
names(data) <- c("category", "count")

# Compute percentages
data$fraction <- data$count / sum(data$count)

# Compute the cumulative percentages (top of each rectangle)
data$ymax <- cumsum(data$fraction)

# Compute the bottom of each rectangle
data$ymin <- c(0, head(data$ymax, n=-1))

# Compute label position
data$labelPosition <- (data$ymax + data$ymin) / 2

# Compute a good label
data$label <- paste0(data$category,"\n" , round((data$count/nrow(merged_df)*100),2), "%" )

# Make the plot
g<-ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
  geom_rect() +
  geom_label( x=3.5, aes(y=labelPosition, label=label), size=3) +
  scale_fill_brewer(palette=4) +
  coord_polar(theta="y") +
  xlim(c(2, 4)) +
  theme_void() +
  theme(legend.position = "none")+   labs(title = "Pathogenicity Class")

pdf("PathogenicityClass.pdf")
print(g)
dev.off()

#------ M location -----


MutationsLocation_counts <- table(merged_df$Mutations.Location)

# Convert state counts to a data frame
data_mutationsLocation <- as.data.frame(MutationsLocation_counts)
names(data_mutationsLocation) <- c("category", "count")

# Compute percentages
data_mutationsLocation$fraction <- data_mutationsLocation$count / sum(data_mutationsLocation$count)

# Compute the cumulative percentages (top of each rectangle)
data_mutationsLocation$ymax <- cumsum(data_mutationsLocation$fraction)

# Compute the bottom of each rectangle
data_mutationsLocation$ymin <- c(0, head(data_mutationsLocation$ymax, n=-1))

# Compute label position
data_mutationsLocation$labelPosition <- (data_mutationsLocation$ymax + data_mutationsLocation$ymin) / 2

# Compute a good label
data_mutationsLocation$label <- paste0(data_mutationsLocation$category,"\n" , round((data_mutationsLocation$count/nrow(merged_df)*100),2), "%" )

# Make the plot
g<-ggplot(data_mutationsLocation, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
  geom_rect() +
  geom_label( x=3.5, aes(y=labelPosition, label=label), size=3) +
  scale_fill_brewer(palette=4) +
  coord_polar(theta="y") +
  xlim(c(2, 4)) +
  theme_void() +
  theme(legend.position = "none")+   labs(title = "Mutations Location")

pdf("MutationsLocation.pdf")
print(g)
dev.off()

#------- based on location Benign  ----------

subset_Benign <- subset(merged_df, `pathogenicity class` == "Benign")


state_counts_benign <- table(subset_Benign$Mutations.Location)

# Convert state counts to a data frame
data_benign <- as.data.frame(state_counts_benign)
names(data_benign) <- c("location", "count")

# Compute percentages
data_benign$fraction <- data_benign$count / sum(data_benign$count)

# Compute the cumulative percentages (top of each rectangle)
data_benign$ymax <- cumsum(data_benign$fraction)

# Compute the bottom of each rectangle
data_benign$ymin <- c(0, head(data_benign$ymax, n=-1))

# Compute label position
data_benign$labelPosition <- (data_benign$ymax + data_benign$ymin) / 2

# Compute a good label
data_benign$label <- paste0(data_benign$location,"\n" , round((data_benign$count/nrow(subset_Benign))*100,2), "%" )

# Make the plot
g<-ggplot(data_benign, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=location)) +
  geom_rect() +
  geom_label( x=3.5, aes(y=labelPosition, label=label), size=3) +
  scale_fill_brewer(palette=4) +
  coord_polar(theta="y") +
  xlim(c(2, 4)) +
  theme_void() +
  theme(legend.position = "none")+   labs(title = "Mutations Location - Benign")

pdf("MutationsLocationBenign.pdf")
print(g)
dev.off()
#----------- Pathogenic -----

subset_Pathogenic <- subset(merged_df, `pathogenicity class` == "Pathogenic")


state_counts_pathogenic <- table(subset_Pathogenic$Mutations.Location)

# Convert state counts to a data frame
data_pathogenic <- as.data.frame(state_counts_pathogenic)
names(data_pathogenic) <- c("location", "count")

# Compute percentages
data_pathogenic$fraction <- data_pathogenic$count / sum(data_pathogenic$count)

# Compute the cumulative percentages (top of each rectangle)
data_pathogenic$ymax <- cumsum(data_pathogenic$fraction)

# Compute the bottom of each rectangle
data_pathogenic$ymin <- c(0, head(data_pathogenic$ymax, n=-1))

# Compute label position
data_pathogenic$labelPosition <- (data_pathogenic$ymax + data_pathogenic$ymin) / 2

# Compute a good label
data_pathogenic$label <- paste0(data_pathogenic$location,"\n" , round((data_pathogenic$count/nrow(subset_Pathogenic))*100,2), "%" )

# Make the plot
g<-ggplot(data_pathogenic, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=location)) +
  geom_rect() +
  geom_label( x=3.5, aes(y=labelPosition, label=label), size=3) +
  scale_fill_brewer(palette=4) +
  coord_polar(theta="y") +
  xlim(c(2, 4)) +
  theme_void() +
  theme(legend.position = "none")+   labs(title = "Mutations Location - Pathogenic")


pdf("MutationsLocationPathogenic.pdf")
print(g)
dev.off()

#----------- Ambiguous -----

subset_Ambiguous <- subset(merged_df, `pathogenicity class` == "Ambiguous")


state_counts_ambiguous <- table(subset_Ambiguous$Mutations.Location)

# Convert state counts to a data frame
data_ambiguous  <- as.data.frame(state_counts_ambiguous)
names(data_ambiguous) <- c("location", "count")

# Compute percentages
data_ambiguous$fraction <- data_ambiguous$count / sum(data_ambiguous$count)

# Compute the cumulative percentages (top of each rectangle)
data_ambiguous$ymax <- cumsum(data_ambiguous$fraction)

# Compute the bottom of each rectangle
data_ambiguous$ymin <- c(0, head(data_ambiguous$ymax, n=-1))

# Compute label position
data_ambiguous$labelPosition <- (data_ambiguous$ymax + data_ambiguous$ymin) / 2

# Compute a good label
data_ambiguous$label <- paste0(data_ambiguous$location,"\n" , round((data_ambiguous$count/nrow(subset_Ambiguous))*100,2), "%" )

# Make the plot
g<-ggplot(data_ambiguous, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=location)) +
  geom_rect() +
  geom_label( x=3.5, aes(y=labelPosition, label=label), size=3) +
  scale_fill_brewer(palette=4) +
  coord_polar(theta="y") +
  xlim(c(2, 4)) +
  theme_void() +
  theme(legend.position = "none")+   labs(title = "Mutations Location - Ambiguous")

pdf("MutationsLocationAmbiguous.pdf")
print(g)
dev.off()