#################################################################################################
#  Encoding "article" and "author" variable in to factors.

library(readxl);
library(WriteXLS);
library(dplyr);
library(tidyr);


#################################################################################################
# Read in the second tab of the sheet because it contains all articles.

tab <- read_excel("Belman-Wolfson Minimum wage data for Andrews and Kasy.xls", 
                  sheet = "Data for Wolfson-Belman SSRN WP");

# Encode "article" to a factor.
article_factor <- match(tab$article, unique(tab$article));
tab <- cbind.data.frame(tab, article_factor);
length(unique(tab$article))


# Generate a unique list of article, and keep authors and year information.
list <- tab %>% 
  distinct(article, authors) %>%
  select(article, article_factor, authors, year) 

wolfson_belman_ssrn <- list

WriteXLS( wolfson_belman_ssrn, ExcelFileName = "wolfson_belman_ssrn unique article.xls"); 

#################################################################################################
# Manually code up "article_id" variable that contains citation of the article (authors, year)
# as the unique identifier for the article.

# Manually fill in the article's publication status: 
# published: = 1 if the article was published in a peer-reviewed journal.
# published_journal: name of the journal published.
# publication_year_current: the most current year in which the article is published.
# working_paper: = 1 if the article is a working paper, or a non peer-reviewed report.

# Read in the worksheet with the above publication status and merge the status back to the
# original worksheet (with all specifications from each article).

pub_tab <- as.data.frame(read_excel("Belman-Wolfson Minimum wage study publication status.xls"));

wolfson_belman_ssrn_pub_status <- right_join(x=pub_tab, y=tab, 
                                            by=c("article", "article_factor","year","authors"));

try(if(nrow(wolfson_belman_ssrn_pub_status) != nrow(tab)) stop("we lost some rows!"))

# Read in the first tab of the worksheet with more results from the specification
# and merge the publication status back to the original worksheet.


tab1 <- read_excel("Belman-Wolfson Minimum wage data for Andrews and Kasy.xls", 
                  sheet = "Data for Belman-Wolfson (2014)");

# The article year differs from the other tab.  So remove the year from tab 2 before merging
pub_tab_key <- select(pub_tab,article,authors, article_factor);
belman_wolfson_2014_pub_status <- right_join(x=pub_tab_key, y=tab1, 
                                            by=c("article","authors"));

try(if(is.element(NA, belman_wolfson_2014_pub_status$article_factor)) stop("article not matched!"))

# An article changed name from "The Effect of MWs on Wages & Employment" in 2010 to
# "The Effect of MWs on Labor Market Outcomes" in 2012, which is article # 26
# We manually match them

belman_wolfson_2014_pub_status$article_factor[is.na(belman_wolfson_2014_pub_status$article_factor)]<- 26;
pub_tab_full <- select(pub_tab,-article,-authors,-year);

belman_wolfson_2014_pub_status <- left_join(x=belman_wolfson_2014_pub_status, y=pub_tab_full,
                                             by=c("article_factor"));

try(if(nrow(belman_wolfson_2014_pub_status) != nrow(tab1)) stop("we lost some rows!"))

#################################################################################################
WriteXLS(c("belman_wolfson_2014_pub_status", "wolfson_belman_ssrn_pub_status"), 
         ExcelFileName = "Belman-Wolfson Minimum wage encoded.xls",
         SheetNames = c("Data for Belman-Wolfson (2014)","Data for Wolfson-Belman SSRN WP")); 


