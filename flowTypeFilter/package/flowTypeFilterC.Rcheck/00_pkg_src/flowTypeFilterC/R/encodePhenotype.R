#Method to encoding a phenotype string into something human-readable:
if (!isGeneric("encodePhenotype")) 
{
  setGeneric("encodePhenotype", function(pheno.string, marker.names) standardGeneric("encodePhenotype"))
}

setMethod("encodePhenotype", signature(pheno.string="character", marker.names="character"),
          function(pheno.string, marker.names)
          {
            
            #MarkerNames<- ThePhenotypes@MarkerNames[ThePhenotypes@PropMarkers]
            #PartitionsPerMarker<- ThePhenotypes@PartitionsPerMarker
            
            pasteOneMarker <- function(markerName)
            {
              markerMatch <- regexpr(paste(markerName, '[+|-]', sep=''), pheno.string)
              if(markerMatch[1] == -1)
                return('0')
              markerStr <- substr(pheno.string, markerMatch[1], markerMatch[1]+attr(markerMatch, 'match.length')-1)
              modStr <- sub(markerName, '', markerStr)
              if(modStr=='-')
                return('1')
              else
                return(as.character(nchar(modStr)+1))
            }
            
            paste(sapply(marker.names, pasteOneMarker), collapse='')
                                
            
          })


## Example of applying to a whole Phenotypes object:
#sapply(res@PhenoCodes, function(x,y){encodePhenotype(y,x)}, res)
