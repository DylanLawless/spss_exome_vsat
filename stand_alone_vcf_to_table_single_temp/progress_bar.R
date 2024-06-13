# Progress bar ----
# calculate percentage of completion
percent_complete <- round(f / length(file_list) * 100, 2)

# calculate number of symbols to print
num_symbols <- round(percent_complete/2)

# print progress bar
cat("\nFinished processing sample: ", f, "of ", length(file_list), "\n",
	 "\r[", paste(rep("=", num_symbols), collapse = ""), ">",
	 paste(rep(" ", 50 - num_symbols), collapse = ""), "]",
	 percent_complete, "%\n")
