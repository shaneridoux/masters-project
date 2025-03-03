# 250301
# Shane Ridoux
# Text Me Updates Script

library(RCurl)
textme <- function(api, project, channel, event, description){
headers = c(
  "Content-Type" = "application/json",
  "Authorization" = paste0("Bearer ",api)
)
params = paste0("{
  \"project\": \"",project,"\",
  \"channel\": \"",channel,"\",
  \"event\": \"",event," is complete\",
  \"description\": \"",description,"\",
  \"icon\": \"🔥\",
  \"notify\": true
}")
res <- postForm("https://api.logsnag.com/v1/log", .opts=list(postfields = params, httpheader = headers, followlocation = TRUE), style = "httppost")
cat(res)
}
# textme(api = "a74417411cf1cb3145ee97bc673ea190", project = "test", channel = "crawler", event = "Test Text Me", description = "This is a test for the textme function.")
