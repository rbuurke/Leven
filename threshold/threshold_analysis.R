dia_w <- read.delim("~/PHD/Leven/threshold/threshold_output_diareg_var_word.txt")
dia_l <- read.delim("~/PHD/Leven/threshold/threshold_output_diareg_var_loc.txt")
gtrp_w <- read.delim("~/PHD/Leven/threshold/threshold_output_gtrp_var_word.txt")
gtrp_l <- read.delim("~/PHD/Leven/threshold/threshold_output_gtrp_var_loc.txt")


library(mgcv)
library(itsadug)


add_words = gam(corr ~ s(loc_ssize, word_ssize), data = dia_w)
fvisgam(add_words, view = c('word_ssize', 'loc_ssize'))
table(dia_w$loc_ssize, dia_w$word_ssize)
boxplot(dia_w$corr ~ dia_w$word_ssize)

add_locations = gam(corr ~ s(loc_ssize, word_ssize), data = dia_l)
fvisgam(add_locations, view = c('loc_ssize', 'word_ssize'))
table(dia_l$loc_ssize, dia_l$word_ssize)
boxplot(dia_l$corr ~ dia_l$loc_ssize)

add_words = gam(corr ~ s(loc_ssize, word_ssize), data = gtrp_w)
fvisgam(add_words, view = c('word_ssize', 'loc_ssize'))
table(dia_l$loc_ssize, dia_l$word_ssize)
boxplot(dia_l$corr ~ dia_l$loc_ssize)

add_locations = gam(corr ~ s(loc_ssize, word_ssize), data = gtrp_l)
fvisgam(add_locations, view = c('loc_ssize', 'word_ssize'))
table(dia_l$loc_ssize, dia_l$word_ssize)
boxplot(dia_l$corr ~ dia_l$loc_ssize)
