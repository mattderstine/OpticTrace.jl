
#a standalone script to load the Zemax browser with all the glass catalogs loaded
using OpticTrace # load all the lens design stuff
loadRICatalog!( "glass/schott")
loadRICatalog!("glass/hikari")
loadRICatalog!("glass/hoya")
loadRICatalog!("glass/ohara")
loadAGFCatalog!("/Users/matt/Development/Software/Zemax/Glasscat")
zemaxBrowser(".")