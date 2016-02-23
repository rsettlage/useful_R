args <- commandArgs(trailingOnly = F)
myargument <- args[length(args)]
filemask <- sub("-","",myargument)

print(filemask)

files<-dir(pattern=filemask)
for(i in 1:length(files)){
	print(files[i])
}