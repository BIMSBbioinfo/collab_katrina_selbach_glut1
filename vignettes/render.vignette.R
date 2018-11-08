library(rmarkdown)

args = commandArgs(trailingOnly=TRUE)

pathToVignette <- args[1] # provide path to the vignette you'd like to render.
dataFolder <- args[2] #provide path to the `data` folder in the downloaded repository ()
outDir <- getwd()

rmarkdown::render(input = pathToVignette, 
		output_dir = outDir, 
        	intermediates_dir = outDir, 
		output_format = rmarkdown::html_document(
						code_folding = 'show',
						toc = TRUE, 
            					toc_float = TRUE, 
						theme = "spacelab", 
						number_sections = TRUE),
		params = list(workdir = outDir, datadir = dataFolder), 
		quiet = FALSE)
