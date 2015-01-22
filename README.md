# RDataInfrastructure

This project provides capabilities to automatically download projects from data repositories (GEO, etc).  The data is then stored in a standard formats depending on its type.

Installation
--------------

You must first get a token from [github](https://github.com/settings/applications). Click generate new token and then save it to a text file (you only get to copy it once). Read that text file in and install repo.  You mush have access to this repo granted to you in the first place since it is private.

```{R}
library(devtools)

auth=read.delim2("githubToken.txt", header=FALSE)

install_github("Selventa/RDataInfrastructure", auth_token = as.character(unlist(auth)))
```
