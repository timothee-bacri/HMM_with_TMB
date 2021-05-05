library(readxl)
library(readr)
library(tidyverse)
library(pbapply)
library(lubridate)
library(ggplot2)

# rm(list = ls())
# files <- list.files(recursive = TRUE)
# i = 0
# hospital_data <- list()
# pb <- txtProgressBar(min = 0, max = length(files), style = 3)
# for(file in files) {
#   tryCatch(
#     hospital_data[[file]] <- read_excel(file, range = cell_cols("C")) %>%
#       bind_cols(read_excel(file, range = cell_cols("L"))) %>%
#       distinct(.keep_all = TRUE),
#     error = function(e) {message(file)}
#   )
#   if (is.null(hospital_data[[file]])) {
#     tryCatch(
#       hospital_data[[file]] <- read_delim(file,
#                                           "\t",
#                                           col_types = cols_only(IPP = "c", SDATE_ARR = col_datetime(format = "%d/%m/%Y %R")),
#                                           escape_double = FALSE,
#                                           trim_ws = TRUE) %>%
#         distinct(.keep_all = TRUE),
#       error = function(e) {message(file)}
#     )
#   }
#   if (is.null(hospital_data[[file]])) (message(paste("NULL", file)))
#   setTxtProgressBar(pb, i)
#   i <- i + 1
# }

# Dates are number of days since 30/12/1899 for Excel
# Hours are percentages of 24*60 minutes in a day

# Start week on Monday
options("lubridate.week.start" = 1)

begin <- as.Date("01-01-2010", format = "%d-%m-%Y")
end <- as.Date("31-12-2019", format = "%d-%m-%Y")
nb_days <- as.integer(end-begin+1)
days <- seq(begin, end, by = "day")

# all_minutes <- tibble(YEAR = rep(year(days), each = 24*60),
#                       MONTH = rep(month(days), each = 24*60),
#                       DAY = rep(mday(days), each = 24*60),
#                       HOUR = rep(rep(0:23, each = 60), nb_days),
#                       MINUTE = rep(0:59, 24*nb_days)) %>%
#   mutate_if(is.numeric, as.integer)

all_hours <- tibble(YEAR = rep(year(days), each = 24),
                    MONTH = rep(month(days), each = 24),
                    DAY = rep(mday(days), each = 24),
                    HOUR = rep(0:23, nb_days)) %>%
  mutate_if(is.numeric, as.integer)

# all_days <- tibble(YEAR = year(days),
#                    MONTH = month(days),
#                    DAY = mday(days)) %>%
#   mutate_if(is.numeric, as.integer)
# 
# all_wdays <- tibble(YEAR = year(days),
#                     WDAY = wday(days)) %>%
#   distinct() %>%
#   mutate_if(is.numeric, as.integer)

# lhour <- as_tibble(read.csv("grouped_hospital_data.csv", stringsAsFactors = FALSE)) %>%
#   mutate(DATE = ymd_hms(DATE))

file <- "O:/Downloads/Book2.xlsx"
sheets <- excel_sheets(file)[- 1]

t <- pblapply(sheets, function(x) {
  read_excel(file, sheet = x, range = cell_cols("C:O"), col_types = rep("text", 13)) %>%
    select(IPP, DATES, HEURE_ENTREE)
})

l <- t %>%
  bind_rows() %>%
  mutate(DATE = ymd(DATES),
         HEURE = paste0(floor(as.numeric(HEURE_ENTREE) * 24),
                        ":",
                        round((as.numeric(HEURE_ENTREE) * 24 * 60) %% 60)),
         CDATE = ymd_hm(paste(DATE, HEURE))) %>%
  select(IPP, CDATE) %>%
  distinct() %>%
  mutate(YEAR = year(CDATE),
         MONTH = month(CDATE),
         DAY = day(CDATE),
         WDAY = wday(CDATE),
         WEEK = week(CDATE),
         HOUR = hour(CDATE),
         MINUTE = minute(CDATE)) %>%
  select(CDATE, YEAR, MONTH, DAY, WEEK, WDAY, HOUR, MINUTE) %>%
  mutate_if(is.numeric, as.integer)

# lday <- l %>%
#   group_by(YEAR, MONTH, DAY) %>%
#   summarise(PATIENTS = n()) %>%
#   right_join(all_days) %>%
#   replace_na(list(PATIENTS = integer(1))) %>%
#   mutate(DATE = ymd(paste(YEAR, MONTH, DAY)),
#          WDAY = as.integer(wday(DATE))) %>%
#   select(DATE, YEAR, MONTH, DAY, WDAY, PATIENTS)
# ggplot(lday, aes(x = PATIENTS)) + geom_histogram(binwidth = 0.5) + ggtitle("Per day")

lhour <- l %>%
  group_by(YEAR, MONTH, DAY, HOUR) %>%
  summarise(PATIENTS = n()) %>%
  right_join(all_hours) %>%
  replace_na(list(PATIENTS = integer(1))) %>%
  mutate(DATE = ymd_h(paste(YEAR, MONTH, DAY, HOUR)),
         WDAY = as.integer(wday(DATE))) %>%
  select(DATE, YEAR, MONTH, DAY, WDAY, HOUR, PATIENTS)

write.csv(lhour, "data/grouped_hospital_data.csv", row.names = FALSE)
# lhour <- as_tibble(read.csv("data/grouped_hospital_data.csv",
#                             stringsAsFactors = FALSE)) %>%
#   mutate(DATE = ymd_hms(DATE))

# lhour %>%
#   select(HOUR, PATIENTS) %>%
#   mutate(HOUR = as.factor(HOUR)) %>%
#   ggplot(aes(x = HOUR, y = PATIENTS)) +
#   geom_boxplot() +
#   ggtitle("Per hour")


# ggplot(lhour, aes(x = PATIENTS)) + geom_histogram(binwidth = 0.5) + ggtitle("Per hour")
# 
# lday %>%
#   group_by(DAY) %>%
#   summarize(avg = mean(PATIENTS)) %>%
#   ggplot(aes(x = DAY, y = avg)) + geom_line() + geom_point() + ggtitle("Per day avg")
# 
# lday %>%
#   group_by(MONTH) %>%
#   summarize(avg = mean(PATIENTS)) %>%
#   ggplot(aes(x = MONTH, y = avg)) + geom_line() + geom_point() + ggtitle("Per month avg")
# 
# lhour %>%
#   group_by(HOUR) %>%
#   summarize(avg = mean(PATIENTS)) %>%
#   ggplot(aes(x = HOUR, y = avg)) + geom_line() + geom_point() + ggtitle("Per hour")
# 
# lhour %>%
#   group_by(WDAY) %>%
#   summarize(avg = mean(PATIENTS)) %>%
#   ggplot(aes(x = WDAY, y = avg)) + geom_line() + geom_point() + ggtitle("Per wday")
# 
# lwday <- l %>%
#   group_by(YEAR, WEEK, WDAY) %>%
#   summarise(PATIENTS = n())
# ggplot(lwday, aes(x = WDAY, y = PATIENTS, group = WDAY)) + geom_boxplot() + ggtitle("Per wday")

# day : 3 652
# hour : 87 648 groups
# minute : 808 596 groups

# CHECK NA
# sapply(t, function(e){sum(is.na(e[, 3]))})
# sum(is.na(l[, 3]))
# 
# ggplot(lhour, aes(x = DATE, y = PATIENTS)) + geom_line()
# 
# tmblist <- list()
# M <- 1:8
# for (m in M) {
#   lambda <- seq(5, 25, length.out = m)
#   gamma <- matrix(0.2 / (m - 1), nrow = m, ncol = m)
#   diag(gamma) <- 0.8
#   # DM.estimate(x = lhour$patients,
#   #             m = m,
#   #             lambda0 = lambda,
#   #             gamma0 = gamma,
#   #             method = METHOD)
#   
#   tgamma = Gamma_n2w(m, gamma)
#   tlambda = log(lambda)
#   parameters = list(tlambda = tlambda, tgamma = tgamma)
#   TMB_data = list(x = lhour$PATIENTS, m = m)
#   tmblist[[m]] <- TMB.estimate(TMB_data = TMB_data,
#                                parameters = parameters,
#                                gradient = TRUE,
#                                hessian = TRUE,
#                                method = METHOD,
#                                std_error = TRUE)
# }
# plot(M, unlist(lapply(tmblist, `[[`, "mllk")), type = "o")
# plot(M, unlist(lapply(tmblist, `[[`, "AIC")), type = "o")
# plot(M, unlist(lapply(tmblist, `[[`, "BIC")), type = "o")
# 
# 1-tmblist[[2]]$BIC/tmblist[[1]]$BIC
# 1-tmblist[[3]]$BIC/tmblist[[2]]$BIC
# 1-tmblist[[4]]$BIC/tmblist[[3]]$BIC