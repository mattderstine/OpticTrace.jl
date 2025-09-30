#=

    lens thorlabs.jl
    definitions of thorlabs parts

=#

export lensAC508180AB, lensAC127050A, lens_ACL12708U, lens_TLF220APC, lens_TLF357775_405
export lens_TLAC254_060





function lensAC508180AB(base, dir, lambda; order = "forward", lensname = "AC508180AB")
    #compute the refractive indexes


    curvAC508180_1= 6.923078165548776121E-03
    glassAC508180_1 = "N-LAK22"
    thickAC508180_1 = 9.5
    semiDiamAC508180 = 25.4
    curvAC508180_2= -8.663404145164059142E-03
    glassAC508180_2 = "N-SF6"
    thickAC508180_2 = 4.
    curvAC508180_3= -3.046912099599362322E-03
    thickAC508180 = thickAC508180_1 + thickAC508180_2
    bflAC508180= 173.09544577349999

    riN_LAK22 = getRefractiveIndexFunc(dirBaseRefractiveIndex, "glass/schott/N-LAK22.yml")
    riN_SF6 = getRefractiveIndexFunc(dirBaseRefractiveIndex, "glass/schott/N-SF6.yml")
    riLens1 = riN_LAK22(lambda)
    riLens2 = riN_SF6(lambda)

    if (order == "forward" )
        base1 = base + thickAC508180_1 .* dir
        base2 = base1 + thickAC508180_2 .* dir

        lens = [
        refractSphere("$(lensname)_1", base, dir, refIndexDefault,riLens1,
            curvAC508180_1, semiDiamAC508180, "testcoat"),
        refractSphere("$(lensname)_2", base1, dir, riLens1,riLens2,
            curvAC508180_2, semiDiamAC508180, "testcoat"),
        refractSphere("$(lensname)_3", base2, dir, riLens2, refIndexDefault,
            curvAC508180_3, semiDiamAC508180, "testcoat")
        ]
    elseif (order == "reverse" )
        base1 = base + thickAC508180_2 .* dir
        base2 = base1 + thickAC508180_1 .* dir

        lens = [
        refractSphere("$(lensname)_3", base, dir, refIndexDefault, riLens2,
            -curvAC508180_3, semiDiamAC508180, "testcoat"),
        refractSphere("$(lensname)_2", base1, dir, riLens2, riLens1,
            -curvAC508180_2, semiDiamAC508180, "testcoat"),
        refractSphere("$(lensname)_1", base2, dir, riLens1, refIndexDefault,
            -curvAC508180_1, semiDiamAC508180, "testcoat")
        ]
    else
        error("lensAC508180AB: order must be 'forward' or 'reverse'")
    end
    lens
end

const curvAC127050A_1=  3.654970760233920000E-002
const glassAC127050A_1 = "N-BK7"
const thickAC127050A_1 = 3.5
const semiDiamAC127050A = 6.35
const curvAC127050A_2= -4.436557231588290200E-002
const glassAC127050A_2 = "SF2"
const thickAC127050A_2 = 1.5
const curvAC127050A_3= -1.088968746596970100E-002
thickAC127050A = thickAC127050A_1 + thickAC127050A_2
const bflAC127050A= 4.721088792319E+1



function lensAC127050A(base, dir, lambda; order = "forward", lensname = "AC127050A")
    #compute the refractive indexes
    riN_SF2 = getRefractiveIndexFunc(dirBaseRefractiveIndex, "glass/schott/N-SF2.yml")
    riN_BK7 = getRefractiveIndexFunc(dirBaseRefractiveIndex, "glass/schott/N-BK7.yml")
    
    riLens1 = riN_BK7(lambda)
    riLens2 = riN_SF2(lambda)

    if (order == "forward" )
        base1 = base + thickAC127050A_1 .* dir
        base2 = base1 + thickAC127050A_2 .* dir

        lens = [
        refractSphere("$(lensname)_1", base, dir, refIndexDefault,riLens1,
            curvAC127050A_1, semiDiamAC127050A, "testcoat"),
        refractSphere("$(lensname)_2", base1, dir, riLens1,riLens2,
            curvAC127050A_2, semiDiamAC127050A, "testcoat"),
        refractSphere("$(lensname)_3", base2, dir, riLens2, refIndexDefault,
            curvAC127050A_3, semiDiamAC127050A, "testcoat")
        ]
    elseif (order == "reverse" )
        base1 = base + thickAC127050A_2 .* dir
        base2 = base1 + thickAC127050A_1 .* dir

        lens = [
        refractSphere("$(lensname)_3", base, dir, refIndexDefault, riLens2,
            -curvAC127050A_3, semiDiamAC127050A, "testcoat"),
        refractSphere("$(lensname)_2", base1, dir, riLens2, riLens1,
            -curvAC127050A_2, semiDiamAC127050A, "testcoat"),
        refractSphere("$(lensname)_1", base2, dir, riLens1, refIndexDefault,
            -curvAC127050A_1, semiDiamAC127050A, "testcoat")
        ]
    else
        error("lensAC127050A: order must be 'forward' or 'reverse'")
    end
    lens
end

const curvAC127019AB_1=  6.998984237749111825E-02
const glassAC127019AB_1 = "N-LAK10"
const thickAC127019AB_1 = 4
const semiDiamAC127019AB = 6.35
const curvAC127019AB_2= -7.266197974165333751E-02
const glassAC127019AB_2 = "N-SF57"
const thickAC127019AB_2 = 0.99980241223690003
const curvAC127019AB_3= -1.459515106998706020E-02
# thickAC127019AB = thickAC127019AB_1 + thickAC127019AB_2
const bflAC127019AB= 15.8738491298




function lensAC127019AB(base, dir, lambda; order = "forward", lensname = "AC127019AB")
    #compute the refractive indexes
    riN_SF57 = getRefractiveIndexFunc(dirBaseRefractiveIndex, "glass/schott/N-SF57.yml")
    riN_LAK10 = getRefractiveIndexFunc(dirBaseRefractiveIndex, "glass/schott/N-LAK10.yml")
    riLens1 = riN_LAK10(lambda)
    riLens2 = riN_SF57(lambda)

    if (order == "forward" )
        base1 = base + thickAC127019AB_1 .* dir
        base2 = base1 + thickAC127019AB_2 .* dir

        lens = [
        refractSphere("$(lensname)_1", base, dir, refIndexDefault,riLens1,
            curvAC127019AB_1, semiDiamAC127019AB, "testcoat"),
        refractSphere("$(lensname)_2", base1, dir, riLens1,riLens2,
            curvAC127019AB_2, semiDiamAC127019AB, "testcoat"),
        refractSphere("$(lensname)_3", base2, dir, riLens2, refIndexDefault,
            curvAC127019AB_3, semiDiamAC127019AB, "testcoat")
        ]
    else
        base1 = base + thickAC127019AB_2 .* dir
        base2 = base1 + thickAC127019AB_1 .* dir

        lens = [
        refractSphere("$(lensname)_3", base, dir, refIndexDefault, riLens2,
            -curvAC127019AB_3, semiDiamAC127019AB, "testcoat"),
        refractSphere("$(lensname)_2", base1, dir, riLens2, riLens1,
            -curvAC127019AB_2, semiDiamAC127019AB, "testcoat"),
        refractSphere("$(lensname)_1", base2, dir, riLens1, refIndexDefault,
            -curvAC127019AB_1, semiDiamAC127019AB, "testcoat")
        ]
    end
    lens
end

const semiDiamACL12708U = 12.7/2
const thickACL12708U = 7.5
const bflACL12708U = 3.7
const osurf2C_ACL12708U = -6.390020150881667300E-002
const osurf1C_ACL12708U = 2.103879571162499400E-001
const osurf1ϵ_ACL12708U = conicToϵ(-1.205070933941 )
const osurf1ASP_ACL12708U = [ 5.3324182908E-4, 1.116288700581E-5,  -3.745566605609E-7, -7.634201667258E-9, 1.360220018284E-10]

const semiDiamACL108U = 10.0/2
const thickACL108U = 5.8
const bflACL108U = 4.
const osurf2C_ACL108U = 0.0
const osurf1C_ACL108U = 2.389486260454002600E-001
const osurf1ϵ_ACL108U = conicToϵ( -6.027E-1) 
const osurf1ASP_ACL108U = [ 2.21E-4]

const semiDiamACL1512U= 15.0/2
const thickACL1512U = 8.
const bflACL1512U = 7.
const osurf2C_ACL1512U = 0.0
const osurf1C_ACL1512U = 1.593117731400350600E-001
const osurf1ϵ_ACL1512U = conicToϵ( -6.139E-1 )
const osurf1ASP_ACL1512U = [0., 0., 0., 6.8E-5]

const semiDiamACL1815U= 18.0/2
const thickACL1815U = 8.2
const bflACL1815U = 10.
const osurf2C_ACL1815U = 0.0
const osurf1C_ACL1815U = 1.279099513942184900E-001
const osurf1ϵ_ACL1815U = conicToϵ( -1.8165) 
const osurf1ASP_ACL1815U = [ 2.933E-4]

const semiDiamACL2018U= 20.0/2
const thickACL2018U = 8.
const bflACL2018U = 13.
const osurf2C_ACL2018U = 0.0
const osurf1C_ACL2018U = 1.062089767827176800E-001
const osurf1ϵ_ACL2018U = conicToϵ( -6.3916E-1 )
const osurf1ASP_ACL2018U = [1.7E-5]

const semiDiamACL25416U= 25.4/2
const thickACL25416U = 14.
const bflACL25416U = 7.3
const osurf2C_ACL25416U = -1.428581996804618000E-002
const osurf1C_ACL25416U = 1.134018728240511100E-001
const osurf1ϵ_ACL25416U = conicToϵ( -9.991715230203E-1 )
const osurf1ASP_ACL25416U = [ 8.68216737142E-5, 0.,  6.376012336251E-8, 0., 2.407308400812E-9, 0., -1.718902114185E-11]

const semiDiamLP357775 = 6.33/2
const thickP357775 = 2.9
const bflP357775 = 1.862919035852 + .25/1.5+ .3
const osurf2C_P357775 = 0.0
const osurf1C_P357775 = 3.497173380949716859E-01
const osurf1ϵ_P357775 = conicToϵ( -1.035706996427) 
const osurf1ASP_P357775 = [0.0027459145825400001, 0., 3.9169060086760001e-05, 0., -2.0239651316750001e-08, 0., -1.867865008283e-07]
const glassP357775="D-LAK6M"



const thickSM1L03 = 8.4
const thickSM1L05 = 13.5
const thickSM1L10 = 26.2
const thickSM1L15 = 38.9
const thickSM1L20 = 51.6
const thickSM1L25 = 64.3
const thickSM1L30 = 77.0
const thickSM1L35 = 89.7
const thickSM1L40 = 102.4

#thicknesses are positive

#LED characteristics

const a375e=[1.0,1.0,0.9968000000000001,0.9872,0.9776,0.9648,0.9488,0.9328, 0.9136000000000001,0.9039999999999999,0.8880000000000001,0.8623999999999999,0.8464000000000002,0.8208,0.7984,0.7728,0.7375999999999999,0.696,0.6511999999999999,0.5967999999999999,0.5423999999999999,0.4847999999999999,0.4336,0.36960000000000004,0.30879999999999996,0.24480000000000002,0.1935999999999999,0.13599999999999987,0.08799999999999997,0.03999999999999989,0.024000000000000042,0.0015999999999999322]
const a375a=[0.,0.028302636518827057,0.08962501564295212,0.1509473947670772,0.2169868799776732,0.27830925910179827,0.35378295648533675,0.4103882295229909,0.4717106086471159,0.5235987755982986,0.5754869425494815,0.6509606399330204,0.6981317007977316,0.7688882920947991,0.8207764590459818,0.8679475199106934,0.919835686861876,0.9717238538130593,1.0188949146777708,1.0707830816289534,1.1085199303207227,1.1604080972719053,1.1981449459636746,1.2453160068283864,1.2830528555201557, 1.3396581285578097,1.3821120833360498,1.4245660381142906,1.4670199928925307,1.500039735497829,1.5189081598437133,1.5330594781031273]

const a455e=[1.,0.9895104895104895,0.9790209790209788,0.965034965034965,0.9580419580419581,0.9335664335664335,0.9125874125874126,0.8846153846153845,0.8601398601398601,0.8111888111888111,0.7517482517482518,0.6783216783216782,0.576923076923077,0.47202797202797203,0.38461538461538464,0.27972027972027974,0.2027972027972027,0.13986013986013987,0.09440559440559433,0.07692307692307684,0.06293706293706282,0.048951048951048994,0.03846153846153842,0.027972027972027858]
const a455a=[0.0,0.05416539057913449,0.13154451997789787,0.1779719976171557,0.2243994752564141,0.28114417014884024,0.3275716477880986,0.36884051680077196,0.4152679944400303,0.4616954720792881,0.5132815583451301,0.5648676446109724,0.631929556756568,0.693832860275578,0.7505775551680048,0.8176394673136004,0.8692255535794426,1.0085079864972162,1.1787420711744963,1.2870728523327641,1.3592933731049441,1.4263552852505386,1.4985758060227177,1.5707963267948966]

const a505e=[1.,0.9933993399339934,0.9867986798679867,0.9768976897689767,0.9570957095709569,0.9438943894389439,0.9174917491749175,0.8811881188118811,0.8481848184818482,0.8151815181518152,0.7755775577557755,0.7458745874587459,0.7095709570957096,0.6633663366336634,0.6105610561056105,0.5577557755775577,0.49174917491749176,0.4158415841584158,0.3366336633663366,0.2541254125412541,0.16171617161716179,0.10561056105610572,0.10231023102310233,0.059405940594059466,0.013201320132013215]
const a505a=[0.0,0.06088357855794187,0.16803867681991916,0.2264869122355433,0.2946765202204381,0.36773681448996787,0.44566779504413334,0.5333401481675698,0.6015297561524644,0.6745900504219944,0.7427796584068889,0.7963572075378782,0.8499347566688669,0.9181243646537613,0.9960553452079267,1.0642449531928213,1.1373052474623515,1.224977600585788,1.298037894855318,1.361356816555577,1.429546424540472,1.487994659956096,1.4977360325253666,1.531830836517814,1.5707963267948966]

const a617e=[1.,0.9928315412186381,0.9749103942652331,0.9498207885304659,0.9175627240143369,0.8817204301075269,0.8458781362007168,0.8064516129032258,0.7777777777777778,0.7383512544802866,0.6953405017921145,0.6523297491039427,0.5663082437275985,0.4910394265232973,0.40860215053763443,0.3369175627240143,0.26881720430107525,0.21505376344085994,0.18279569892473094,0.16129032258064507,0.125448028673835,0.09677419354838704,0.07526881720430076,0.05376344086021489,0.05017921146953384,0.03225806451612901]
const a617a=[0.,0.05271128613405693,0.12123595810833114,0.17394724424238805,0.23192965898985077,0.2846409451239077,0.3373522312579646,0.39006351739202155,0.432232546299267,0.4744015752065124,0.5165706041137579,0.5481973757941924,0.6009086619282493,0.6430776908354948,0.6799755911293344,0.7168734914231741,0.7537713917170143,0.7959404206242598,0.9119052501191851,1.017327822387299,1.1438349091090352,1.2756131244441773,1.3915779539391027,1.4811871403669994,1.5075427834340283,1.5707963267948966]

const a730e=[1.,0.9857651245551602,0.9537366548042704,0.9252669039145907,0.896797153024911,0.8434163701067616,0.7900355871886122,0.7437722419928825,0.6868327402135233,0.6441281138790035,0.5943060498220639,0.5266903914590749,0.4448398576512454,0.3451957295373663,0.26334519572953724,0.1957295373665482,0.15658362989323837,0.11387900355871866,0.0818505338078291,0.04982206405693953,0.03202846975088957,0.017793594306049956]
const a730a=[0.,0.06829549246934283,0.1365909849386863,0.19963297798731114,0.25742147161521667,0.3414774623467167,0.4202799536574975,0.4938289455475596,0.5673779374376217,0.6251664310655273,0.6777014252727148,0.7407434183213389,0.7985319119492452,0.8668274044185887,0.9246158980464942,1.013925388198713,1.1347558748752433,1.266093360393211,1.3816703476490226,1.4867403360633973,1.549782329112022,1.5760498262156153]

const a940e=[1.,0.9966329966329965,0.9865319865319866,0.9696969696969696,0.936026936026936,0.888888888888889,0.8451178451178452,0.7878787878787877,0.7239057239057239,0.6464646464646465,0.5622895622895623,0.47811447811447816,0.37710437710437716,0.2861952861952861,0.21212121212121224,0.14814814814814817,0.10774410774410775,0.07407407407407408,0.05723905723905734,0.04377104377104387,0.023569023569023666]
const a940a=[0.,0.04960409453036531,0.15432384965002494,0.2314857744750374,0.3141592653589795,0.4023443223018509,0.46848311500900436,0.5401334737750872,0.6062722664822411,0.6889457573661831,0.7495729840144069,0.8212233427804903,0.9094083997233617,0.9920818906073032,1.0692438154323158,1.1519173063162573,1.2456139293180588,1.3448221183787887,1.4219840432038016,1.4770997037930957,1.5707963267948966]


function lens_ACL12708U(base, dir, wl)
    [
        refractSphere("ACL12708i", base, dir, refIndexDefault,wl,
            -osurf2C_ACL12708U, semiDiamACL12708U, "testcoat"),
        refractEvenAsphere("ACL12708o", base .+thickACL12708U*dir, dir, wl, refIndexDefault,
             -osurf1C_ACL12708U, osurf1ϵ_ACL12708U, .-osurf1ASP_ACL12708U, semiDiamACL12708U, "testcoat")
    ]
end

const semiDiamF220APC = 2.75
const thickF220APC = 5.031595527853
#const bflF220APC = 4.
const osurf1C_F220APC = 0.0
const osurf2C_F220APC = -1.555661881368497800E-001
const osurf2ϵ_F220APC = conicToϵ(-7.31283599511E-1 )
const osurf2ASP_F220APC = [0.0, 0. ,0. , -8.924167336705E-5,0.,-4.384364140114E-7]
const glass_F220APC="D-ZK3"


function lens_TLF220APC(base, dir,  λ; order = "forward", lensname = "TL_F220APC")
    riD_ZK3 = getRefractiveIndexFunc(dirBaseRefractiveIndex, "glass/cdgm/D-ZK3.yml")

    lensEASinglet(base, dir, 0.0, 0.0, [0.], osurf2C_F220APC, osurf2ϵ_F220APC, osurf2ASP_F220APC, thickF220APC, λ, riD_ZK3, semiDiamF220APC; order = order, lensname = lensname)
end

# From TL Zemax file: 357775-405-Zemax(ZMX).zmx 

const semiDiamTL357775_405 = 2.4
const thickTL357775_405 = 2.8983318458129999
export semiDiamTL357775_405, thickTL357775_405


"""
    lens_TLF357775_405(base, dir, λ; order = "forward", lensname = "TL_357775_405")
    
    Creates a lens with the characteristics of the Thorlabs TL357775-405 lens.
    located at base in the direction dir with wavelength λ.

    Option parameters are order ("forward" or "reverse") and lensname (default "TL_357775_405").

"""
function lens_TLF357775_405(base, dir,  λ; order = "forward", lensname = "TL_357775_405")

    osurf1C_TL357775_405 = -3.497173380949716859E-01
    osurf1ϵ_TL357775_405 = conicToϵ(-1.035706996427)
    osurf1ASP_TL357775_405 = [-0.0027459145825400001, -3.9169060086760001e-05, 2.0239651316750001e-08, 1.867865008283e-07]
    riD_LAK6 = getRefractiveIndexFunc(dirBaseRefractiveIndex, "glass/cdgm/D-LAK6.yml")

    lensEASinglet(base, dir, 0.0, 0.0, [0.], osurf1C_TL357775_405, osurf1ϵ_TL357775_405, osurf1ASP_TL357775_405, thickTL357775_405, λ, riD_LAK6, semiDiamTL357775_405; order = order, lensname = lensname)
end

"""
    lens_TLAC254_060(base, dir, λ; order = "forward", lensname = "TL_AC254-060")
    
    Creates a lens with the characteristics of the Thorlabs AC254-060
    located at base in the direction dir with wavelength λ.

    Option parameters are order ("forward" or "reverse") and lensname (default "TL_254_060").

"""
function lens_TLAC254_060(base, dir,  λ; order = "forward", lensname = "TL_AC254-060")

    surf1C = 2.398656752218759900E-002
    surf2C = -3.863987635239570000E-002    
    surf3C = -4.334633723450400300E-003
    surf1t = 8.0
    surf2t = 2.5
    semiDiam = 25.4/2
    glass1 = "glass/hoya/BAF11.yml"
    glass2 = "glass/hoya/E-FD10.yml"

    ri_glass1 = getRefractiveIndexFunc(dirBaseRefractiveIndex, glass1)
    ri_glass2 = getRefractiveIndexFunc(dirBaseRefractiveIndex, glass2)
    riLens1 = ri_glass1(λ)
    riLens2 = ri_glass2(λ)
    if (order == "forward" )
        base1 = base + surf1t .* dir
        base2 = base1 + surf2t .* dir

        lens = [
        refractSphere("$(lensname)_1", base, dir, refIndexDefault,riLens1,
            surf1C, semiDiam, "testcoat"),
        refractSphere("$(lensname)_2", base1, dir, riLens1,riLens2,
            surf2C, semiDiam, "testcoat"),
        refractSphere("$(lensname)_3", base2, dir, riLens2, refIndexDefault,
            surf3C, semiDiam, "testcoat")
        ]
    elseif (order == "reverse" )
        base1 = base + surf2t .* dir
        base2 = base1 + surf1t .* dir

        lens = [
        refractSphere("$(lensname)_3", base, dir, refIndexDefault, riLens2,
            -surf3C, semiDiam, "testcoat"),
        refractSphere("$(lensname)_2", base1, dir, riLens2, riLens1,
            -surf2C, semiDiam, "testcoat"),
        refractSphere("$(lensname)_1", base2, dir, riLens1, refIndexDefault,
            -surf1C, semiDiam, "testcoat")
        ]
    else
        error("lens_TLAC254_060: order must be 'forward' or 'reverse'")
    end
    lens
end

