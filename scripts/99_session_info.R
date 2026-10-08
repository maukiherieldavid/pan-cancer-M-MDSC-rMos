# Capture the software environment after the analysis is run.
source(file.path("scripts", "00_setup.R"))
writeLines(capture.output(sessionInfo()), file.path(RESULTS_DIR, "sessionInfo.txt"))
