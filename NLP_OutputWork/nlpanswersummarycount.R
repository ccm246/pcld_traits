# looking at NLP question output success

library(dplyr)
library(ggplot2)
library(tidyr)

a <- read.csv('answers_with_summary.csv', header=T)

head(a)

#just looking at knowns
a %>%
  filter(known == 1) %>%            
  count(question.last.word, sort = TRUE)

a %>%
  filter(known == 1) %>%
  count(question.last.word) %>%
  ggplot(aes(x = reorder(question.last.word, n), y = n)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(x = "Last Word in Question", y = "Number of Known Answers", title = "Known Answers by Question Type") +
  theme_minimal()

length(unique(a$question.last.word))
unique(a$question.last.word)
