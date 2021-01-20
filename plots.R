library(readr)

cor_vec_rnd <- read_delim("pdf/PHD/LEVEN_TEST/Leven/cor_vec_rnd.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
cor_vec_gtrp <- read_delim("pdf/PHD/LEVEN_TEST/Leven/cor_vec_gtrp.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
cor_vec_diareg_om <- read_delim("pdf/PHD/LEVEN_TEST/Leven/cor_vec_diareg_om.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
cor_vec_diareg_jv <- read_delim("pdf/PHD/LEVEN_TEST/Leven/cor_vec_diareg_jv.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

boxplot(cor_vec_rnd)
boxplot(cor_vec_gtrp)
boxplot(cor_vec_diareg_om)
boxplot(cor_vec_diareg_jv)

vec_list = vector()
for (sample in cor_vec_rnd) {
  vec_list = c(vec_list, var(sample))
}

plot(vec_list)

vec_list = vector()
for (sample in cor_vec_gtrp) {
  vec_list = c(vec_list, var(sample))
}

plot(vec_list)

vec_list = vector()
for (sample in cor_vec_diareg_om) {
  vec_list = c(vec_list, var(sample))
}

plot(vec_list)

vec_list = vector()
for (sample in cor_vec_diareg_jv) {
  vec_list = c(vec_list, var(sample))
}

plot(vec_list)
