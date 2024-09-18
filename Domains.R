
library(ggplot2)
library(writexl)
library(ggplot2)
library(ggrepel)


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

filtered_df <- subset(merged_df, frequency >4)

# View the filtered data frame
print(filtered_df)

mut.df <- data.frame("AA" = filtered_df$Residue.Position, 
                     "Mut" = filtered_df$Mutant.Residue, 
                     "Type" = filtered_df$`pathogenicity class`, 
                     "Freq" = filtered_df$frequency)

domain.df <- data.frame("Feature" = c("Start", "PIK3-ABD", "PIK3-RBD", "C1", "Nuclear", "PIK3-Catalytic", "End"), 
                        "Type" = c("str", "dom", "dom","dom","dom","dom", "str"),
                        "Start" = c(1, 26, 194, 327,410, 772, 1070), 
                        "End" = c(1, 115, 285,496,418,1053, 1070))

str.fill <- "#E1E1E1"
str.col <- "#16161D"

dom.fill <- c("PIK3-ABD" = "#CD0BBC", "PIK3-RBD" = "#AA4BAB", "C1" = "tan", "Nuclear" = "#F5C710", "PIK3-Catalytic"="#28E2E5")
dom.col <- c("#16161D")



gp <- ggplot() +
  geom_rect(data = subset(head(domain.df,n=1), Type == "str"),
            mapping = aes(xmin = 1, xmax = 1070, ymin = 0.3, ymax = 0.7),
            fill = str.fill,
            colour = str.col)

gp <- gp + scale_y_continuous(limits = c(0, 25), breaks = 0:25)

gp <- gp + geom_segment(data = mut.df, 
                        mapping = aes(x = AA, xend = AA, y = 0.7, yend = Freq)) +
  geom_point(data = mut.df,
             mapping = aes(x = AA, y = Freq, fill = Type),
             shape = 21,
             size = 2) +
  geom_text_repel(data = mut.df,
                  mapping = aes(x = AA, y = Freq, label = Mut),
                  bg.colour = "white",
                  seed = 12345,
                  nudge_y = 0.25)



gp <- gp + geom_rect(data = subset(domain.df, Type == "dom"),
                     mapping = aes(xmin = Start, xmax = End, ymin = 0.2, ymax = 1.5, fill = Feature, group = Feature),
                     fill = dom.fill[subset(domain.df, Type == "dom")$Feature],
                     colour = dom.col) +
                     geom_text_repel(data=subset(domain.df, Type == "dom"), aes(x=Start+5, y=1.7, label=Feature), size=4,min.segment.length = unit(0, 'lines'),
                                     nudge_y = 2, fontface = 'bold')
#                     geom_text(data=subset(domain.df, Type == "dom"), aes(x=(Start+End)/2, y=1.7, label=Feature), size=4,check_overlap = T)
gp
gp <- gp +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(linetype = "dotted"),text = element_text(size=20)) +
  labs(x = "Range", y = "Freq.", fill = "Pathogenicity class")

gp

pdf("Domains>4.pdf",width=20)
print(gp)
dev.off()