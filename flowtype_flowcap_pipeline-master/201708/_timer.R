time_to_run = as.POSIXct("2017-05-11 7:00:00")
while(TRUE) {
    Sys.sleep(1)
    if(Sys.time == time_to_run) {
        source("code/06b_dist_classify.R")
	source("code/07_dist_compare.R")
    }
}
