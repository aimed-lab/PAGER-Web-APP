###############
### run app ###


#source("/home/zongliang/DEseq2/lmfit_code/ui.R")
#source("/home/zongliang/DEseq2/lmfit_code/server.R")
#shinyApp(ui = ui, server = server)
#shiny::runApp()
### end run app ###
###############

#########################
### rshiny deplotment ###
#install.packages('rsconnect')
#library(rsconnect)
#rsconnect::setAccountInfo(name='zongyue',
#                          token='E5EE4BA28D8A8564E4D3F19E7DF2F068',
#                          secret='6I0mvfpFgum4LTlC9r+Vw7/t72a1I1H6JK2Ed34I')
#
#rsconnect::deployApp('/home/zongliang/DEseq2/lmfit_code')
## unable-to-deploy-app-unable-to-determine-package-source/
#options(repos = BiocManager::repositories())
#shiny::runApp()
#deployApp()


#https://stackoverflow.com/questions/26222052/how-can-i-specify-which-shiny-account-to-use-when-deploying
options(repos = BiocManager::repositories())
rsconnect::setAccountInfo(name='aimed-lab',
                          token='5CCB428089433B7CD811E43AFB753120',
                          secret='+Z13D47ZqWAMSULmtg/iDypnb11uiAwt+nyottTi')
rsconnect::deployApp('/home/zongliang/DEseq2/FACT')
setwd("/home/zongliang/DEseq2/FACT/") # or set the full path below
rsconnect::deployApp("/home/zongliang/DEseq2/FACT", #yAppNameOffline
                     appName = "FACT", #yAppNameOnline
                     account = "aimed-lab") 


shiny::sign_out_from_shiny(
  session = shiny::getDefaultReactiveDomain(),
  redirect_page = "?page=sign_in"
)
#########################