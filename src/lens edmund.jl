#=

    lens edmund.jl
    definitions of edmund optics parts

=#

export lens_EO68001, lens_EO67652,lens_EO67548

#riN_SF11 = getRefractiveIndexFunc(dirBaseRefractiveIndex, "glass/schott/N-SF11.yml")
#riN_BK7 = getRefractiveIndexFunc(dirBaseRefractiveIndex, "glass/schott/N-BK7.yml")






#=
VERS 171115
DBDT 0 38398 1.000E+10 0.000E+00
DBDT 1 0 1 0 0 0 0
MODE SEQ
NAME 38398
NOTE 0 12.5mm Dia. x 15.0mm FL, NIR II Coated, Plano-Convex Lens
PFIL 0 0 0
LANG 0
UNIT MM X W X CM MR CPMM
ENPD 12.5
ENVD 20 1 0
GFAC 0 0
GCAT SCHOTT OHARA MISC CORNING INFRARED 
SDMA 0 1 0
OMMA 1 1
FTYP 0 0 1 1 0 0 0 1
ROPD 2
HYPR 0
PICB 1
XFLN 0
YFLN 0
FWGN 1
WAVM 1 0.58760000000000001 1
PWAV 1
GLRS 1 0
SURF 0
  COMM 38398
  FIMP 

  CURV 0.0
  DISZ INFINITY
  MEMA 0 0 0 0 1 ""
SURF 1
  STOP
  FIMP 

  CURV 8.496176720475789867E-02
  COAT EO_NIRII_785
  DISZ 3.1000000000000001
  GLAS N-SF11
  DIAM 5.75 1 0 0 1 ""
  OEMA 0.5 0 0 0 0 ""
  MEMA 6.25 0 0 0 1 ""
  CLAP 0 5.75 0
SURF 2
  FIMP 

  CURV 0.0
  COAT EO_NIRII_785
  DISZ 13.262120360019384
  MAZH 0 0
  DIAM 5.75 1 0 0 1 ""
  OEMA 0.5 0 0 0 0 ""
  MEMA 6.25 0 0 0 1 ""
  FLAP 0 5.75 0
SURF 3
  FIMP 

  CURV 0.0
  DISZ 0
  MEMA 0 0 0 0 1 ""

=#


function lens_EO38398(base, dir, lambda; order = "forward", lensname = "EO38-398")
  riN_SF11 = getRefractiveIndexFunc(dirBaseRefractiveIndex, "glass/schott/N-SF11.yml")
  lensSinglet(base, dir, 8.496176720475789867E-02, 0.0, 3.1, lambda, riN_SF11, 5.75; order = order, lensname = lensname)
end

#=

VERS 171115
DBDT 0 68001 1.000E+10 0.000E+00
DBDT 1 0 1 0 0 0 0
MODE SEQ
NAME 68001
NOTE 0 25mm Diameter x -50 FL, NIR II Coated, Plano-Concave Lens
PFIL 0 0 0
LANG 0
UNIT MM X W X CM MR CPMM
ENPD 25
ENVD 20 1 0
GFAC 0 0
GCAT SCHOTT OHARA MISC CORNING INFRARED 
SDMA 0 1 0
OMMA 1 1
FTYP 0 0 1 1 0 0 0 1
ROPD 2
HYPR 0
PICB 1
XFLN 0
YFLN 0
FWGN 1
WAVM 1 0.58760000000000001 1
PWAV 1
GLRS 1 0
SURF 0
  COMM 68001
  FIMP 

  CURV 0.0
  DISZ INFINITY
  MEMA 0 0 0 0 1 ""
SURF 1
  STOP
  FIMP 

  CURV -3.869969040247679681E-02
  COAT EO_NIRII_517
  DISZ 3.5
  GLAS N-BK7
  DIAM 12 1 0 0 1 ""
  OEMA 0.5 0 0 0 0 ""
  MEMA 12.5 0 0 0 1 ""
  CLAP 0 12 0
SURF 2
  FIMP 

  CURV 0.0
  COAT EO_NIRII_517
  DISZ -52.307642959813563
  MAZH 0 0
  DIAM 12 1 0 0 1 ""
  OEMA 0.5 0 0 0 0 ""
  MEMA 12.5 0 0 0 1 ""
  FLAP 0 12 0
SURF 3
  FIMP 

  CURV 0.0
  DISZ 0
  MEMA 0 0 0 0 1 ""
=#

function lens_EO68001(base, dir, lambda; order="forward", lensname = "EO68-001")
  riN_BK7 = getRefractiveIndexFunc(dirBaseRefractiveIndex, "glass/schott/N-BK7.yml")  
  lensSinglet(base, dir, -3.869969040247679681E-02, 0.0, 3.5, lambda, riN_BK7, 12.5; order = order, lensname = lensname)
    
end

#=
VERS 171115
DBDT 0 67548 1.000E+10 0.000E+00
DBDT 1 0 1 0 0 0 0
MODE SEQ
NAME 67548
NOTE 0 25.0mm Dia. x 75.0mm FL, NIR II Coated, Plano-Convex Lens
PFIL 0 0 0
LANG 0
UNIT MM X W X CM MR CPMM
ENPD 25
ENVD 20 1 0
GFAC 0 0
GCAT SCHOTT OHARA MISC CORNING INFRARED 
SDMA 0 1 0
OMMA 1 1
FTYP 0 0 1 1 0 0 0 1
ROPD 2
HYPR 0
PICB 1
XFLN 0
YFLN 0
FWGN 1
WAVM 1 0.58760000000000001 1
PWAV 1
GLRS 1 0
SURF 0
  COMM 67548
  FIMP 

  CURV 0.0
  DISZ INFINITY
  MEMA 0 0 0 0 1 ""
SURF 1
  STOP
  FIMP 

  CURV 2.579979360165119903E-02
  COAT EO_NIRII_517
  DISZ 4.5
  GLAS N-BK7
  DIAM 12 1 0 0 1 ""
  OEMA 0.5 0 0 0 0 ""
  MEMA 12.5 0 0 0 1 ""
  CLAP 0 12 0
SURF 2
  FIMP 

  CURV 0.0
  COAT EO_NIRII_517
  DISZ 72.033451490638825
  MAZH 0 0
  DIAM 12 1 0 0 1 ""
  OEMA 0.5 0 0 0 0 ""
  MEMA 12.5 0 0 0 1 ""
  FLAP 0 12 0
SURF 3
  FIMP 

  CURV 0.0
  DISZ 0
  MEMA 0 0 0 0 1 ""
=#


function lens_EO67548(base, dir, lambda; order = "forward", lensname = "EO67-548")
  riN_BK7 = getRefractiveIndexFunc(dirBaseRefractiveIndex, "glass/schott/N-BK7.yml")
  lensSinglet(base, dir, 2.579979360165119903E-02, 0.0, 4.5, lambda, riN_BK7, 12.5; order = order, lensname = lensname)
    
end


#=
VERS 171115
DBDT 0 67652 1.000E+10 0.000E+00
DBDT 1 0 1 0 0 0 0
MODE SEQ
NAME 67652
NOTE 0 25mm Dia. x 75mm FL, NIR II Coated, Double-Convex Lens
PFIL 0 0 0
LANG 0
UNIT MM X W X CM MR CPMM
ENPD 25
ENVD 20 1 0
GFAC 0 0
GCAT SCHOTT OHARA MISC CORNING INFRARED 
SDMA 0 1 0
OMMA 1 1
FTYP 0 0 1 1 0 0 0 1
ROPD 2
HYPR 0
PICB 1
XFLN 0
YFLN 0
FWGN 1
WAVM 1 0.58760000000000001 1
PWAV 1
GLRS 1 0
SURF 0
  COMM 67652
  FIMP 

  CURV 0.0
  DISZ INFINITY
  MEMA 0 0 0 0 1 ""
SURF 1
  STOP
  FIMP 

  CURV 1.304461257500649958E-02
  COAT EO_NIRII_517
  DISZ 3.5
  GLAS N-BK7
  DIAM 12 1 0 0 1 ""
  OEMA 0.5 0 0 0 0 ""
  MEMA 12.5 0 0 0 1 ""
  CLAP 0 12 0
SURF 2
  FIMP 

  CURV -1.304461257500649958E-02
  COAT EO_NIRII_517
  DISZ 73.586785848153951
  MAZH 0 0
  DIAM 12 1 0 0 1 ""
  OEMA 0.5 0 0 0 0 ""
  MEMA 12.5 0 0 0 1 ""
  FLAP 0 12 0
SURF 3
  FIMP 

  CURV 0.0
  DISZ 0
  MEMA 0 0 0 0 1 ""
=#

function lens_EO67652(base, dir, lambda; order="forward", lensname = "EO67-652")
  riN_BK7 = getRefractiveIndexFunc(dirBaseRefractiveIndex, "glass/schott/N-BK7.yml")

  lensSinglet(base, dir, 1.304461257500649958E-02, -1.304461257500649958E-02, 3.5, lambda, riN_BK7, 12.5; order = order, lensname = lensname)

end
