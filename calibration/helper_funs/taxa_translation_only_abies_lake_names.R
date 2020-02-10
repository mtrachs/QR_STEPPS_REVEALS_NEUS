taxa <- read.csv(paste(data.loc,l.name,'taxon_names.csv',sep=''))
taxa <- taxa$x

species.english <- c("Ash","Beech","Birch","Hemlock","Hickory","Maple",
                    "Oak","Pine","Spruce","Tamarack","Chestnut","Other conifer","Other hardwood")

#These taxa were used as other conifers and other hardwood in Dawson et al.
other_conifer <- c('Abies')#,'Juniperus','Cedrus','Cupressaceae','Taxus','Thuja','Toxicodendron')

other_hardwood <- c('Alnus','Celtis','Morus','Platanus','Prunus','Rhamnus','Rhus','Robinia','Sambucus',
                    'Shepherdia','Corylus','Juglans','Magnolia','Moraceae','Morus','Nyssa','Ostrya',#'Carpinus', remove Carpinus as it only occurs as Ostrya.Carpinus
                    'Rhamnaceae','Salix','Tilia',"Ulmus","Populus")# Ericaceae

#-------------------------------------------------------------------------------------------------------------------
species.latin <- c("Fraxinus","Fagus","Betula","Tsuga","Carya","Acer",
                  "Quercus","Pinus","Picea","Larix","Castanea",'Conifer','Deciduae')


taxa.translate <- as.vector(matrix(nrow=length(taxa)))
taxa.english <- as.vector(matrix(nrow=length(taxa)))

  for(x in species.latin[!species.latin%in%c('Conifer','Deciduae')]) {
    col_idx <- grep(x,taxa)
    taxa.translate[col_idx] <- x
    taxa.english[col_idx] <- species.english[x==species.latin]
  }

  for (x in other_conifer) {
    col_idx <- grep(x,taxa)
    taxa.translate[col_idx] <- 'Conifer'
    taxa.english[col_idx] <- 'Other conifer'
  }

for (x in other_hardwood) {
  col_idx <- grep(x,taxa)
  taxa.translate[col_idx] <- 'Deciduae'
  taxa.english[col_idx] <- 'Other hardwood'
}

write.csv(data.frame(target = as.character(taxa),latin=taxa.translate,match=taxa.english),
          paste(data.loc,l.name,'taxon_names_translated.csv',sep=''),row.names = FALSE)

          