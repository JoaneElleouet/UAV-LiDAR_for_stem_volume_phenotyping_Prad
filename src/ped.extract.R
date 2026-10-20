ped.extract<-function(x, reference){ #would be great to have a column with depth for individual trees
  rpbc.ped <- reference
  i=0
  ped.length<-c(0) #initialise pedigree length
  keep_going = TRUE

  while(keep_going)
  {
    i<-i+1
    temp<-unique(c(rpbc.ped[which(rpbc.ped$Tree%in%x),'Tree'], #get trees and their parents in 1 vector
                   rpbc.ped[which(rpbc.ped$Tree%in%x),'Mom'],
                   rpbc.ped[which(rpbc.ped$Tree%in%x),'Dad']))
    ped.length[i]<-length(temp) #record pedigree length (how many trees are picked up in this round)
    x<-temp

    if(i>1){ #if not first round
      if(ped.length[i]==ped.length[i-1]){ #if same number is picked up as previous round then finish
        keep_going <- FALSE
      }
    }

  }

  temp.ped<-rpbc.ped[which(rpbc.ped$Tree%in%temp),]
  return(list(pedigree = temp.ped,
              ped_depth = i)) #returning i to see depth of pedigree - can drop later
}
