# simu1time <- as.POSIXct("2021-05-22 01:10:58 CEST")
# now <- Sys.time()
# idx_coverage <- 190
# idx_sample <- 1

# end_benchmark <- as.numeric(simu1time) + 1000
# time <- as.numeric(now) - end_benchmark
# time <- time / ((idx_coverage-1)*500+idx_sample) * 500*200
# endtime <- end_benchmark + time
# endtime <- as.POSIXct(endtime, origin = "1970-01-01 00:00.00 UTC")
# print(paste("end estimated at:", endtime))

# start <- lambtime
start <- simu1time
now <- Sys.time()

print(idx_coverage)
print(start)

# SIMU1
# benchmark_duration <- sum(method_comparison_df_simu1$time)
# bootstrap_duration <- 29 # Measured

# SIMU2
benchmark_duration <- sum(method_comparison_df_simu2$time)
bootstrap_duration <- 176

start_coverage <- start + benchmark_duration + bootstrap_duration
current_duration_coverage <- now - start_coverage
# duration <=> idx_coverage
# x <=> COVERAGE_SAMPLES
# x = duration / idx_coverage * COVERAGE_SAMPLES
duration_coverage <- current_duration_coverage / (idx_coverage-0.5) * COVERAGE_SAMPLES
end_coverage <- start_coverage + duration_coverage
end_coverage <- as.POSIXct(end_coverage, origin = "1970-01-01 00:00.00 UTC")
print(end_coverage)