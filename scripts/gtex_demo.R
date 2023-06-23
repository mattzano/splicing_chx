#transcr <- read.delim(file="~/Downloads/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct.gz", skip=2)
#head(transcr)

# install CePa package from CRAN
install.packages('CePa')
# Attached require library
library(CePa)
# Load in gct file
Data <- read.gct("~/Downloads/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct.gz") 

header <- read.table("~/Downloads/header_gtex.txt", skip = 2)
file <- read.table("~/Downloads/output_GTEx.txt")

names(file) <- header
file2 <- file %>% column_to_rownames("transcript_id")
file3 <- file2[,c(2:17383)]

transcr <- as.data.frame(t(file3))


meta_gtex <- read.table("~/Desktop/gtex.txt", sep = "\t", header = T)
subjd <- read.table("~/Desktop/SUBJID .txt", sep = "\t", header = T)

transcr2 <- transcr %>%
  rownames_to_column("SAMPID") %>%
  left_join(meta_gtex)

transcr3 <- transcr2 %>%
  separate(col = SAMPID, into = c("foo", "bar"), sep = "-", remove = F) %>%
  mutate(SUBJID = paste0(foo, "-", bar)) %>%
  #mutate(SUBJID = gsub("^-*","",SAMPID)) %>%
  #mutate(subjd = strsplit(SAMPID, "-"))
  #mutate(SUBJID2 = gsub("^.*-*-","",SAMPID)) #%>%
  left_join(subjd, by = "SUBJID")

transcr3$inclusion = rowSums(transcr3[,c(4,7,8,11)])
transcr3$exclusion = rowSums(transcr3[,c(5,6,9,12)])

transcr3 %>%
  dplyr::filter(!is.na(SMTS)) %>%
  dplyr::filter(SMTS != "") %>%
  ggplot(aes(x = inclusion, y = exclusion)) +
  geom_smooth(method='lm', se = F) +
  geom_point() +
  facet_wrap(facets = vars(SMTS)) +
  theme_classic()

#head(transcr3)




