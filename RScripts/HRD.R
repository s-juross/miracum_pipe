library(sequenza)
library(scarHRD)

args <- commandArgs()
ID <- args[6]
path <- args[7]
alternative_file <- args[8]

print(args)
print(ID)
print(path)
setwd(path)
getwd()

input_file <- paste0(path, "/", ID, ".small.seqz.gz")
print(input_file)
# Determine the HRD-Score
tryCatch({
		result <- scar_score(input_file, reference = "grch37", seqz = TRUE)
		write.table(x = result, file = paste0(path, "/", ID, "_HRD.txt"))
	}, error = function(e) {
		cat("\nBeim Berechnen des HRD-Scores ist ein Fehler aufgetreten.\n")
	}
)
# CNV analysis of SEQUENZA to determine tumor purity and ploidy
tryCatch({
	extr <- sequenza.extract(input_file, verbose = FALSE)
	cp <- sequenza.fit(extr)
	sequenza.results(
	  sequenza.extract = extr, cp.table = cp, sample.id = ID,
	  out.dir = paste0(path, "/", ID, "_sequenza")
	)
}, error = function(e) {
	cat("\nBeim Berechnen der Variabeln purity und ploidy ist ein Fehler",
		"aufgetreten.\n")
	# hier nun Ersatzfile anlegen für alternative Lösungen anlegen.
	# Wird für Report gebraucht.	
	file_conn <- file(alternative_file, "w")
	writeLines("\"cellularity\"\t\"ploidy\"\t\"SLPP\"", file_conn)
	close(file_conn)
})
