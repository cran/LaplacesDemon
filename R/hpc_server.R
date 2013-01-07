###########################################################################
# server_Listening                                                        #
#                                                                         #
# The server_Listening function is not intended to be called directly by  #
# the user. It is an internal-only function that is intended to prevent   #
# cluster problems while using the INCA algorithm through the             #
# LaplacesDemon.hpc function.                                             #
###########################################################################

server_Listening <- function(n=2, port=19009)
     {
     slist <- vector('list', n)
     for (i in 1:n) {
          slist[[i]] <- socketConnection("localhost", port, server=TRUE,
               open="r+")
          cat("\nClient", i, "Connected")}
     tmp <- NULL
     stop_server <- FALSE
     cat("\nStart listening...")
     repeat
          {
          ready <- which(socketSelect(slist, TRUE))
          for (i in ready) {
               #print(paste("Socket", i, "ready to write"))
               con <- slist[[i]]
               #print("Write message...")
               serialize(tmp, con)
               #print("Read message...")
               buf <- try(unserialize(con), silent=TRUE)
               if(is.matrix(buf)) {
                    tmp <- buf
                    }
               else
                    {
                    stop_server <- TRUE
                    break
                    }
               }
          if(stop_server == TRUE) break
          }
     for (i in 1:n) {
          close(slist[[i]])
          cat("\nClose connection", i)
          }
     }
