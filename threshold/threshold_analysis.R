dia_w <- read.delim("~/PHD/Leven/threshold/threshold_output_diareg_var_word.txt")
dia_l <- read.delim("~/PHD/Leven/threshold/threshold_output_diareg_var_loc.txt")
gtrp_w <- read.delim("~/PHD/Leven/threshold/threshold_output_gtrp_var_word.txt")
gtrp_l <- read.delim("~/PHD/Leven/threshold/threshold_output_gtrp_var_loc.txt")

dia_w_20 <- read.delim("~/PHD/Leven/threshold/threshold_output_diareg_var_word_sub_nl.txt")
dia_l_20 <- read.delim("~/PHD/Leven/threshold/threshold_output_diareg_var_loc_sub_nl.txt")
gtrp_w_20 <- read.delim("~/PHD/Leven/threshold/threshold_output_gtrp_var_word_sub_nl.txt")
gtrp_l_20 <- read.delim("~/PHD/Leven/threshold/threshold_output_gtrp_var_loc_sub_nl.txt")




library(mgcv)
library(itsadug)
library(dplyr)

add_words = gam(corr ~ s(loc_ssize, word_ssize), data = dia_w)
fvisgam(add_words, view = c('word_ssize', 'loc_ssize'), main = "Adding random sets of words (DiaReg)")
table(dia_w$loc_ssize, dia_w$word_ssize)
boxplot(dia_w$corr ~ dia_w$word_ssize)

add_locations = gam(corr ~ s(loc_ssize, word_ssize), data = dia_l)
fvisgam(add_locations, view = c('loc_ssize', 'word_ssize'), main = "Adding random sets of locations (DiaReg)")
table(dia_l$loc_ssize, dia_l$word_ssize)
boxplot(dia_l$corr ~ dia_l$loc_ssize)

add_words = gam(corr ~ s(loc_ssize, word_ssize), data = gtrp_w)
fvisgam(add_words, view = c('word_ssize', 'loc_ssize'), main = "Adding random sets of words (GTRP)")
table(dia_l$loc_ssize, dia_l$word_ssize)
boxplot(dia_l$corr ~ dia_l$loc_ssize)

add_locations = gam(corr ~ s(loc_ssize, word_ssize), data = gtrp_l)
fvisgam(add_locations, view = c('loc_ssize', 'word_ssize'), main = "Adding random sets of locations (GTRP)")
table(dia_l$loc_ssize, dia_l$word_ssize)
boxplot(dia_l$corr ~ dia_l$loc_ssize)

# distributional values

par(mfrow=c(1,2))

summary_dia_l = dia_l %>%
  group_by(loc_ssize, word_ssize) %>%
  summarise(lower_bound = min(corr),
            lower_95 = quantile(corr, probs = 0.05, na.rm=TRUE),
            upper_95 = quantile(corr, probs = 0.95, na.rm=TRUE),
            upper_bound = max(corr))

m_summ = gam(lower_95 ~ s(loc_ssize, word_ssize), data = summary_dia_l)
fvisgam(m_summ, view = c('loc_ssize', 'word_ssize'), main = 'DiaReg 0.05 bound
        [adding locations]')

summary_gtrp_l = gtrp_l %>%
  group_by(loc_ssize, word_ssize) %>%
  summarise(lower_bound = min(corr),
            lower_95 = quantile(corr, probs = 0.05, na.rm=TRUE),
            upper_95 = quantile(corr, probs = 0.95, na.rm=TRUE),
            upper_bound = max(corr))

m_summ = gam(lower_95 ~ s(loc_ssize, word_ssize), data = summary_gtrp_l)
fvisgam(m_summ, view = c('loc_ssize', 'word_ssize'), main = 'GTRP 0.05 bound
        [adding locations]')

summary_dia_w = dia_w %>%
  group_by(loc_ssize, word_ssize) %>%
  summarise(lower_bound = min(corr),
            lower_95 = quantile(corr, probs = 0.05, na.rm=TRUE),
            upper_95 = quantile(corr, probs = 0.95, na.rm=TRUE),
            upper_bound = max(corr))

m_summ = gam(lower_95 ~ s(word_ssize, loc_ssize), data = summary_dia_w)
fvisgam(m_summ, view = c('word_ssize', 'loc_ssize'), main = 'DiaReg 0.05 bound
        [adding words]')


summary_gtrp_w = gtrp_w %>%
  group_by(loc_ssize, word_ssize) %>%
  summarise(lower_bound = min(corr),
            lower_95 = quantile(corr, probs = 0.05, na.rm=TRUE),
            upper_95 = quantile(corr, probs = 0.95, na.rm=TRUE),
            upper_bound = max(corr))

m_summ = gam(lower_95 ~ s(word_ssize, loc_ssize), data = summary_gtrp_w)
fvisgam(m_summ, view = c('word_ssize', 'loc_ssize'), main = 'GTRP 0.05 bound
        [adding words]')





summary_dia_l = dia_l_20 %>%
  group_by(loc_ssize, word_ssize) %>%
  summarise(lower_bound = min(corr),
            lower_95 = quantile(corr, probs = 0.05, na.rm=TRUE),
            upper_95 = quantile(corr, probs = 0.95, na.rm=TRUE),
            upper_bound = max(corr))

m_summ = gam(lower_95 ~ s(loc_ssize, word_ssize), data = summary_dia_l)
fvisgam(m_summ, view = c('loc_ssize', 'word_ssize'), main = 'DiaReg 0.05 bound
        [adding locations]', ylim = c(0, 20))

summary_gtrp_l = gtrp_l_20 %>%
  group_by(loc_ssize, word_ssize) %>%
  summarise(lower_bound = min(corr),
            lower_95 = quantile(corr, probs = 0.05, na.rm=TRUE),
            upper_95 = quantile(corr, probs = 0.95, na.rm=TRUE),
            upper_bound = max(corr), ylim = c(0, 20))

m_summ = gam(lower_95 ~ s(loc_ssize, word_ssize), data = summary_gtrp_l)
fvisgam(m_summ, view = c('loc_ssize', 'word_ssize'), main = 'GTRP 0.05 bound
        [adding locations]')

summary_dia_w = dia_w_20 %>%
  group_by(loc_ssize, word_ssize) %>%
  summarise(lower_bound = min(corr),
            lower_95 = quantile(corr, probs = 0.05, na.rm=TRUE),
            upper_95 = quantile(corr, probs = 0.95, na.rm=TRUE),
            upper_bound = max(corr))

m_summ = gam(lower_95 ~ s(word_ssize, loc_ssize), data = summary_dia_w)
fvisgam(m_summ, view = c('word_ssize', 'loc_ssize'), main = 'DiaReg 0.05 bound
        [adding words]', ylim = c(0, 20))


summary_gtrp_w = gtrp_w_20 %>%
  group_by(loc_ssize, word_ssize) %>%
  summarise(lower_bound = min(corr),
            lower_95 = quantile(corr, probs = 0.05, na.rm=TRUE),
            upper_95 = quantile(corr, probs = 0.95, na.rm=TRUE),
            upper_bound = max(corr))

m_summ = gam(lower_95 ~ s(word_ssize, loc_ssize), data = summary_gtrp_w)
fvisgam(m_summ, view = c('word_ssize', 'loc_ssize'), main = 'GTRP 0.05 bound
        [adding words]', ylim = c(0, 20))

