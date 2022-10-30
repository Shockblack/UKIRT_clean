##################################################################
### Code adapted from perl script provided by D. Nataf in 2019 ###
###    Used originally in Nataf et al. 2011, 2013a, 2013b?     ###
##################################################################

import numpy as np

def RC_MCMC(data,init_EWRC=1.2, init_B=0.50, init_K_RC=14, init_sigma_RC=0.2,iterations=10000,print_output=False):
    """
    EXPECTS FILTERED DATA
    data : Nx2 array, first column K-band magnitudes, second column H minus K (HMK) colors
    init_EWRC : initial guess of equivavlent width N_RC/A, seems to of order 1
    init_B : initial guess for exponent B, should be order 0.5?
    init_K_RC : initial guess for K-band magnitude of RC, added default of 13.5
    init_sigma_RC : initial guess for dispersion in RC K-band, default of 0.2
    iterations : number of MCMC iterations
    print_output : bool flag to print summary

    returns:
    List of determined values [A,Aerr,B,Berr,M_RC,M_RCerr,sigma_RC,sigma_RCerr,N_RC,N_RCerr]    
    """
    # input
    # K-band magnitude of the clump
    clump = init_K_RC 
    #bins for histogram (peak finding?)
    binnumber = 100
    # MCMC iterations
    Iterations = iterations

    #STARTING VALUES, reassignment to existing variables
    M0 = clump
    SIGMA = init_sigma_RC
    B = init_B
    EWRC = init_EWRC
    #don't matter but leaving anyway.
    #Just used a place holders for variables that are set to zero or multiplied by zero in case of single RC 
    SIGMA2 = 0.17
    M02 = clump + 0.25
    clumpfraction = 0.5
    #for two clumps, doesn't matter?
    bumpseparation = 0 #MagDiff
    
    #good magnitudes? previously called "value", changed to K_rcstars
    K_rcstars = data[:,0]
    #good colors? changed to HMK_rcstars
    HMK_rcstars = data[:,1]
    # a status for peak finding?
    status = np.zeros(len(K_rcstars))

    
    #sort K_rcstars by magnitude. Added a "list" suffix so copying for now
    K_rcstars_list = np.sort(K_rcstars)
    xstart = clump- 1.5
    xend = clump + 1.5
    binwidth = (xend-xstart)/binnumber
    totalpop = 0
    bincenter = []
    binpopulation = []
    for i in np.arange(binnumber):
        center = xstart + binwidth*(i-0.5)
        bincenter.append(center)
        binpopulation.append(0)
        for j in range(len(K_rcstars_list)):
            diff = abs(K_rcstars_list[j]-center)
            if ((diff<=0.5*binwidth)):
                binpopulation[i] += 1
                totalpop +=1
                status[j] +=1

    #
    addedstars = 0
    for i in range(len(status)):
        k1 = i
        if status[k1]!=1:
            mindiff = 10000
            mindiff2 = 200000
            for j in range(len(bincenter)):
                k2 = j
                linediff = np.abs(K_rcstars_list[k1] - bincenter[k2])
                if (linediff<mindiff):
                    mindiff = linediff
                    locationfit1 = k2
                elif linediff<mindiff2:
                    mindiff2 = linediff
                    locationfit2 = k2

            if status[i]==0:
                binpopulation[locationfit1] += 0.5
                # binpopulation[locationfit2] += 0.5
                addedstars += 1
            elif status[i]==2:
                binpopulation[locationfit1] -= 0.5
                # binpopulation[locationfit2] -= 0.5
                addedstars -=1
            else:
                print("ERROR: Exiting")
                quit()
    

    # more code
    AGBfrac = 0.028
    AGBDelta = -1.07
    AGBSigma = 0.18
    #assuming one RGB
    bumpfraction = 0.201
    MagDiff = 0.737
    DeltaSIGMAb = 0.00



    #jump sizes for MCMC
    jumpfactor = 0.9/(1+len(K_rcstars)/7000.)
    Bj = 0.04*jumpfactor
    SIGMAj = 0.04*jumpfactor
    M0j = 0.04*jumpfactor
    DiffM0j = 0.04*jumpfactor
    EWRCj = 0.12*jumpfactor
    clumpfractionj = 0.07*jumpfactor
    #bumpfractionj = 0.02
    #bumpseparationj = 0.02
    zloc1 = 1
    zloc1 = 0 ### zloc1 = 1 means use maximum likelihood analysis, 0 means bin data
    chi2min = 1000000
    tempchi2min = 1000000
    chi2 = chi2min
    chi2minsingleclump = 1000000
    acceptance = 0

    #assuming one clump, removing some logic
    singleIterations = 16 * Iterations / 8.0

    chances = 4
    chances = 4
    shifter = np.max([int((Iterations-singleIterations)/chances),3000])
    varnumber = 4-1
    temp = 0
    k3 = 0
    singleclump = 1 ### 1 means do a singleclump fit

    chi2list = []
    Blist = []
    Alist = []
    DeltaM0list = []
    SIGMAlist = []
    EWRClist= []
    fractionlist= []
    M0list= []
    M02list= []
    RCtotallist=[]
    chi2list2 = []
    Blist2 = []
    Alist2 = []
    DeltaM0list2 = []
    SIGMAlist2 = []
    EWRClist2= []
    fractionlist2= []
    M0list2= []
    M02list2= []
    RCtotallist2=[]

    for i in range(Iterations):
        k1 = i
        k3 += 1
        k4 = k1 - singleIterations
        if ( ((k4 % shifter) == 0) and (k4 > 1)):
                chi2 = 10000000
                clumpfraction = 0.3 + 0.6*np.random.rand()
                SIGMA = 0.1+0.3*np.random.rand()
                #SIGMA2 = 0.1+0.3*np.random.rand()
                B = 0.300 + 0.700*np.random.rand()
                #printf ("M0 M02\n")
                DeltaM0 = M02 - M0
                M0 = np.max([M0 - 0.5*np.random.rand(),(xstart+0.1)])
                M02 = M0 + DeltaM0
                #printf ("M0 M02\n\n")
                EWRC = 0.7+0.7*(np.random.rand())
                #printf ("Modifications!\n")

        if ((k4 >= 0) and (singleclump == 1)):
            #printf ("Shifting to double!\n")
            #printf ("singleclump Iterations singleclump\n")
            singleclump = 0
            clumpfraction = 0.5
                
        varselect = 0
        Bprop = B
        SIGMAprop = SIGMA
        M0prop = M0
        #SIGMAprop2 = SIGMA2
        M0prop2 = M02
        EWRCprop = EWRC
        clumpfractionprop = clumpfraction
        bumpseparationprop = bumpseparation
        bumpfractionprop = bumpfraction
                
        if (k1 > 0):
            varchooser = np.random.rand()
            #I think there is some missing code here, no prior referecenes to NRC, NRCj, NRCpass
            #XXX May be relic, as NRC and A were replaced with equivalent width EWRC = NRC/A
            """
            NRCpass = 1
            while (NRCpass == 1):
                NRCpass = 0
                DeltaNRC = np.random.rand(2*NRCj) - NRCj
                if (NRC+DeltaNRC < 0):
                    NRCpass = 1
                else: 
                    NRCprop = NRC + DeltaNRC
            """


            SIGMApass = 1
            while (SIGMApass == 1):
                SIGMApass = 0
                DeltaSigma =  (2*SIGMAj)*np.random.rand() - SIGMAj
                if ( (SIGMA+DeltaSigma > 2.)  or (SIGMA+DeltaSigma < -0.01)  ):
                    SIGMApass = 1
                    #print "fail!\n"
                else: 
                    SIGMAprop = SIGMA + DeltaSigma
            
            M0pass  = 0
            while M0pass == 0 :
                M0pass = 1
                DeltaM0 = 2*M0j*np.random.rand() - M0j
                if ((M0 + DeltaM0) <= xstart) :
                    M0pass = 0
                elif  ((M02 + DeltaM0) >= xend) :
                    M0pass = 0
                else :
                    M0prop = M0 + DeltaM0
                    M0prop2 = M02 + DeltaM0
            
            # Two copies of SIGMApass? Idenitical? Removing this one for now
            """        
            SIGMApass = 1
            while (SIGMApass == 1):
                SIGMApass = 0
                DeltaSigma =  (2*SIGMAj)*np.random.rand() - SIGMAj
                if ( (SIGMA2+DeltaSigma > 1.0)  or (SIGMA2+DeltaSigma < 0.0)  ):
                    SIGMApass = 1
                else:
                    SIGMAprop2 = SIGMA2 + DeltaSigma
            """
    
        
            DeltaPasser = 0
            while (DeltaPasser == 0):
                DeltaPasser = 1
                DeltaDiffM0 = 2*DiffM0j*np.random.rand() - DiffM0j
                if ((M0prop2 - M0prop - DeltaDiffM0) <= 0):
                    DeltaPasser = 0
                elif ((M0prop + DeltaDiffM0/2) <= xstart) :
                    DeltaPasser = 0
                elif ((M0prop2 - DeltaDiffM0/2) >= xend) :
                    DeltaPasser = 0
                elif ( np.abs(M0prop2 - M0prop - DeltaDiffM0 - 0.5) >= 0.50  ):
                    #value1 = M0prop2 - M0prop - DeltaDiffM0
                    #value2 = M0predicted
                    #printf ("value1 value2\n")
                    DeltaPasser = 0
                else:
                    M0prop = M0prop + DeltaDiffM0/2
                    M0prop2 = M0prop2 - DeltaDiffM0/2

            clumppasser = 0

            # k9 is some weird diagnositc thing, does not show up again
            #k9 = 9
            while (clumppasser == 0):
                #k9 += 1
                clumppasser = 1
                Deltaclumpfraction = 2*clumpfractionj*np.random.rand() - clumpfractionj
                #printf clumpfractionprop . "\n"
                if ( (clumpfraction+Deltaclumpfraction) > 1.00):
                    clumppasser = 0
                elif ( ((clumpfraction+Deltaclumpfraction) < 0.00) and (singleclump == 0)):
                    clumppasser = 0
                else :
                    clumpfractionprop = clumpfraction +  Deltaclumpfraction
            
        
            EWpasser = 1
            while (EWpasser == 1):
                EWpasser = 0
                DeltaEWRC = 2*EWRCj*np.random.rand()-EWRCj
                if ((EWRC +DeltaEWRC) <= 0):
                    EWpasser = 1
                else :
                    EWRCprop = EWRC +DeltaEWRC
            
        
            Bpasser = 1
            while (Bpasser == 1):
                Bpasser = 0
                DeltaB = 2*Bj*np.random.rand()-Bj
                #printf ("hello!\n")
                if ((B +DeltaB) < 0.100):
                    Bpasser = 1
                elif ((B+DeltaB) > 1.500):
                    Bpasser = 1
                else: 
                    Bprop = B +DeltaB


        SIGMAprop2 = SIGMAprop
        if ((SIGMAprop*SIGMAprop+DeltaSIGMAb) >0 ):
            SIGMApropB = np.sqrt(SIGMAprop*SIGMAprop + DeltaSIGMAb)
        else :
            SIGMApropB = SIGMAprop
            
        if ((SIGMAprop2*SIGMAprop2+DeltaSIGMAb) >0 ):
            SIGMAprop2B = np.sqrt(SIGMAprop2*SIGMAprop2 + DeltaSIGMAb)
        else :
            SIGMAprop2B = SIGMAprop2
    
        #if (M0prop2 < M0prop):
        #    print ("M0prop M0prop2  DeltaDiffM0 \n")
        if (singleclump == 1):
            clumpfractionprop = 0


        ##########################
        ### Function Integrator###
        ##########################
        LocFunctionSet = []
        Integral = 0
        modifier2 = clumpfractionprop/(1-clumpfractionprop)
        if (zloc1 == 0):
            #modifier2 = clumpfractionprop
            AGBint = 0
            term1 = 0
            for i in range(binnumber-1):
                k2 = i
                x1 = bincenter[k2]-M0prop
                #printf ("hello!\n")
                ###### front population
                LocFunction = np.exp(Bprop*(x1))
                LocFunction += (0.39894228*EWRCprop/SIGMAprop)*np.exp(-0.5*((x1/SIGMAprop)**2))
                x1B = x1 - bumpseparation
                LocFunction += (0.39894228*EWRCprop*bumpfraction/SIGMApropB)*np.exp(-0.5*((x1B/SIGMApropB)**2))
                X1C = x1 - AGBDelta
                LocFunction += (0.39894228*EWRCprop*AGBfrac/SIGMAprop)*np.exp(-0.5*((X1C/SIGMAprop)**2))
                ####### back population
                x2 = bincenter[k2]-M0prop2
                LocFunction += modifier2*np.exp(Bprop*(x2))
                LocFunction += (modifier2*(0.39894228)*EWRCprop/SIGMAprop2)*np.exp(-0.5*((x2/SIGMAprop)**2))
                x2B = x2 - bumpseparation
                LocFunction += (modifier2*(0.39894228)*EWRCprop*bumpfraction/SIGMAprop2B)*np.exp(-0.5*((x2B/SIGMAprop2B)**2))
                X2C = x2 - AGBDelta
                LocFunction += (modifier2*(0.39894228)*EWRCprop*AGBfrac/SIGMAprop2)*np.exp(-0.5*((X2C/SIGMAprop2)**2))
                #push @LocFunctionSet, LocFunction
                LocFunctionSet.append(LocFunction)
                Integral += LocFunction*binwidth
            Aprop = (len(K_rcstars)+1)/Integral
            A1prop = Aprop
            A2prop = A1prop*modifier2
            NRCprop = A1prop*EWRCprop
            NRCprop2 = A2prop*EWRCprop
            Aprop = A1prop*(1+modifier2*np.exp(Bprop*(M0prop-M0prop2)))
    

        ##############################
        ##### Parameter Evaluation ###
        ##############################
        if (zloc1 == 0):
            chi2prop = 0
            for i in range((binnumber)-1):
                k2 = i
                Expected = LocFunctionSet[k2]*binwidth*A1prop
                #printf ("k2 LocFunctionSet[k2] binwidth A1prop Aprop\n")
                #printf ("Integral:\tIntegral\n")
                diffchi = ((Expected-binpopulation[k2])**2)/Expected
                chi2prop += diffchi
            #NEEDED?    
            #undef @LocFunctionSet
        else :
            lprop = 0
            for i in range(len((K_rcstars))):
                k2 = i
                x = K_rcstars[k2]-M0prop
                FN = Aprop*np.exp(Bprop*(x))
                #### bright clump
                Term1 = (0.39894228/SIGMAprop)*NRCprop
                Exponent1 = (-1)*(x)*(x)/(2*SIGMAprop*SIGMAprop)
                FN += Term1*np.exp(Exponent1)
                ## bright bump
                xB = K_rcstars[k2]-(M0prop+bumpseparation)
                Term1B = (0.39894228/SIGMApropB)*NRCprop*bumpfraction
                Exponent1B = (-1)*(xB)*(xB)/(2*SIGMApropB*SIGMApropB)
                FN += Term1B*np.exp(Exponent1B)
                ## bright AGBbump
                xC = K_rcstars[k2]-(M0prop+AGBDelta)
                Term1C = (0.39894228/SIGMAprop)*NRCprop*AGBfrac
                Exponent1C = (-1)*(xC)*(xC)/(2*SIGMAprop*SIGMAprop)
                FN += Term1C*np.exp(Exponent1C)
                
                #### faint clump
                x2 = K_rcstars[k2] - M0prop2
                Term2 = (0.39894228/SIGMAprop2)*NRCprop2
                Exponent2 = (-1)*(x2)*(x2)/(2*SIGMAprop2*SIGMAprop2)
                FN += Term2*np.exp(Exponent2)
                ## faint bump
                x2B = K_rcstars[k2]-(M0prop2+bumpseparation)
                Term2B = (0.39894228/SIGMAprop2B)*NRCprop2*bumpfraction
                Exponent2B = (-1)*(x2B)*(x2B)/(2*SIGMAprop2B*SIGMAprop2B)
                FN += Term2B*np.exp(Exponent2B)
                ## faint AGBbump
                x2C = K_rcstars[k2]-(M0prop2+AGBDelta)
                Term2C = (0.39894228/SIGMAprop2)*NRCprop2*AGBfrac
                Exponent2C = (-1)*(x2C)*(x2C)/(2*SIGMAprop2*SIGMAprop2)
                FN += Term2C*np.exp(Exponent2C)
                if (FN>0):
                    lprop += np.log(FN)
                else :
                    lprop -= 1000000
            
        
            chi2prop = -2*lprop
    

        ##########################
        ##### Priors #############
        ##########################
        ### I've removed priors, other than the requirement that the two RCs have equal dispersion Gaussians
        ########################
        ### Markov #############
        ########################
        Improvement = (0.5)*(chi2-chi2prop)
        Factor = np.exp(Improvement)
        modular = (k1 % shifter)
        passer = np.random.rand()
        if ((chi2prop < chi2minsingleclump) and (singleclump == 1)):
            chi2minsingleclump = chi2prop
            Aoptsingleclump =Aprop
            Boptsingleclump = Bprop
            SIGMAoptsingleclump = SIGMAprop
            NRCoptsingleclump = NRCprop
            M0optsingleclump = M0prop
            # faint clump
            #SIGMAopt2 = SIGMAprop2
            #NRCopt2 = NRCprop2
            #M0opt2 = M0prop2
            #EWRCopt = EWRCprop
            #clumpfractionopt = clumpfractionprop
            
        if (chi2prop < chi2min):
            chi2min = chi2prop
            Aopt =Aprop
            Bopt = Bprop
            SIGMAopt = SIGMAprop
            NRCopt = NRCprop
            M0opt = M0prop
            # faint clump
            SIGMAopt2 = SIGMAprop2
            NRCopt2 = NRCprop2
            M0opt2 = M0prop2
            EWRCopt = EWRCprop
            clumpfractionopt = clumpfractionprop
    

        if (Factor>1):
            A = Aprop
            B = Bprop
            NRC = NRCprop
            SIGMA = SIGMAprop
            M0 = M0prop
            # faint clump
            SIGMA2 = SIGMAprop2
            NRC2 = NRCprop2
            M02 = M0prop2
            EWRC = EWRCprop
            bumpseparation = bumpseparationprop
            bumpfraction = bumpfractionprop
            chi2 = chi2prop
            clumpfraction = clumpfractionprop
            acceptance += 1
            #printf ("Here we are!!!\n")
        else :
            if (Factor>passer):
                modulus = (k1 % shifter)
                A = Aprop
                B = Bprop
                NRC = NRCprop
                SIGMA = SIGMAprop
                M0 = M0prop
                # faint clump
                SIGMA2 = SIGMAprop2
                NRC2 = NRCprop2
                M02 = M0prop2
                EWRC = EWRCprop
                chi2 = chi2prop
                clumpfraction = clumpfractionprop
                acceptance += 1
                #printf ("Here we are!!!\n")
        
        #THIS NEEDS TO BE ON THE END LEVEL OF THE ITERATIONS LOOP
        #need to define empty lists before loop starts
        #printf ("M0prop2\n")
        chi2list.append(chi2prop)
        Blist.append(Bprop)
        Alist.append(Aprop)
        DeltaM0list.append((M0prop2-M0prop))
        SIGMAlist.append(SIGMAprop)
        EWRClist.append(EWRCprop)
        fractionlist.append(clumpfractionprop)
        M0list.append(M0prop)
        M02list.append(M0prop2)
        RCtotallist.append((NRCprop+NRCprop2))
        

    for i in range(len(chi2list)):
        k1 = i
        if (chi2list[k1] <= (chi2min+1)):
            Blist2.append(Blist[k1])
            Alist2.append(Alist[k1])
            DeltaM0list2.append(DeltaM0list[k1])
            SIGMAlist2.append(SIGMAlist[k1])
            EWRClist2.append(EWRClist[k1])
            fractionlist2.append(fractionlist[k1])
            M0list2.append(M0list[k1])
            M02list2.append(M02list[k1])
            RCtotallist2.append(RCtotallist[k1])
    #pdb.set_trace()    
    #Sort thess
    Alist3 = np.sort(Alist2)
    Blist3 = np.sort(Blist2)
    DeltaM0list3 = np.sort(DeltaM0list2)
    SIGMAlist3 = np.sort(SIGMAlist2)
    EWRClist3 = np.sort(EWRClist2)
    fractionlist3 = np.sort(fractionlist2)
    M0list3 = np.sort(M0list2)
    M02list3 = np.sort(M02list2)
    RCtotallist3 = np.sort(RCtotallist2)

    Aerror = 0.5*(Alist3[-1] - Alist3[0])
    Berror = 0.5*(Blist3[-1] - Blist3[0])
    DeltaM0error = 0.5*(DeltaM0list3[-1] - DeltaM0list3[0])
    SIGMAerror = 0.5*(SIGMAlist3[-1] - SIGMAlist3[0])
    EWRCerror = 0.5*(EWRClist3[-1] - EWRClist3[0])
    fractionerror =  0.5*(fractionlist3[-1] - fractionlist3[0])
    M01error = 0.5*(M0list3[-1] - M0list3[0])
    #print ("M02ist3[-1] M02list3[0]\n")
    M02error = 0.5*(M02list3[-1] - M02list3[0])
    RCtotalerror = 0.5*(RCtotallist3[-1] - RCtotallist3[0])
    del chi2list
    del Blist
    del Alist
    del DeltaM0list
    del SIGMAlist
    del EWRClist
    del fractionlist
    del Blist2
    del Alist2
    del M0list
    del M02list
    del DeltaM0list2
    del SIGMAlist2
    del EWRClist2
    del fractionlist2
    del M0list2
    del M02list2
    del Blist3
    del Alist3
    del DeltaM0list3
    del SIGMAlist3
    del EWRClist3
    del fractionlist3
    del M0list3
    del M02list3
    
    DeltaM0opt = M0opt2 - M0opt
    #print SIGMAopt . "\n";
    #exit;
    #open "modelclump", ">modelclump.txt"
    chi2prop = 0
    inputpop = 0
    bumpseparationopt = MagDiff
    bumpfractionopt = bumpfraction
    SIGMAoptB = np.sqrt(SIGMAopt*SIGMAopt + DeltaSIGMAb)
    SIGMAopt2B = np.sqrt(SIGMAopt2*SIGMAopt2 + DeltaSIGMAb)
    totres = 0
    for i in range(binnumber-1):
        k1 = i
        x = bincenter[k1]-M0opt
        #FN = Aopt*((bincenter[k1]/M0opt)**Bopt)
        FN = Aopt*np.exp(Bopt*(x))
        ### bright clump
        ExpectedRG = FN*binwidth
        Term1 = (0.39894228/SIGMAopt)*NRCopt
        Exponent1 = (-1)*(x)*(x)/(2*SIGMAopt*SIGMAopt)
        FN += Term1*np.exp(Exponent1)
        ## bright bump
        xB = bincenter[k1]-(M0opt+bumpseparationopt)
        Term1B = (0.39894228/SIGMAoptB)*NRCopt*bumpfractionopt
        Exponent1B = (-1)*(xB)*(xB)/(2*SIGMAoptB*SIGMAoptB)
        FN += Term1B*np.exp(Exponent1B)
        ExpectedRGBB = ExpectedRG +  Term1B*np.exp(Exponent1B)*binwidth
        ## bright AGBbump
        xC = bincenter[k1]-(M0opt+AGBDelta)
        Term1C = (0.39894228/SIGMAopt)*NRCopt*AGBfrac
        Exponent1C = (-1)*(xC)*(xC)/(2*SIGMAopt*SIGMAopt)
        FN += Term1C*np.exp(Exponent1C)
        ### faint clump
        x2 = x - (M0opt2-M0opt)
        Term2 = (0.39894228/SIGMAopt2)*NRCopt2
        Exponent2 = (-1)*(x2)*(x2)/(2*SIGMAopt2*SIGMAopt2)
        FN += Term2*np.exp(Exponent2)
        
        totres += (Expected-binpopulation[k1])
        #
        chi2prop += diffchi
        inputpop += Expected
        #SAVE OUTPUT
        #print modelclump ("%3.3f %3.3f %g %3.3f %3.3f\n", x,Expected,binpopulation[k1],ExpectedRG, ExpectedRGBB)
    
    #close "modelclump"
    percentage2 = 100*acceptance/Iterations
    chi2report = chi2prop

    Kclump = (M0opt*NRCopt + M0opt2*NRCopt2)/(NRCopt+NRCopt2)
        
    TotalRC = NRCopt + NRCopt2
    DeltaM0opt = M0opt2 - M0opt
    TotalRC = NRCopt + NRCopt2
    deltachi2double = chi2minsingleclump - chi2min

    #if (RCnumber == 1):
    if print_output:
        print ("M0 is:\t\t\t%2.3f +\- %1.3f"%(M0opt, M01error))
        print ("Total RC is:\t\t%4.1f +\- %4.1f"%(TotalRC, RCtotalerror))
        print ("SIGMAopt:\t\t%1.3f +\- %1.3f"%(SIGMAopt, SIGMAerror))
        print ("Aopt:\t\t\t%1.3f +\- %1.3f"%(Aopt, Aerror))
        print ("Bopt:\t\t\t%1.3f +\- %1.3f"%(Bopt, Berror))
        print ("Chi Square for binnumber bins, 5 parameters and 1 constraint:\t%3.3f"%(chi2min))
        M0opt2 = M0opt

    output_dict = {'A':Aopt, 'Aerror':Aerror, 'B':Bopt, 'Berror':Berror, 'MRC':M0opt, 'MRCerror':M01error, 'SIGMA':SIGMAopt, 'SIGMAerror':SIGMAerror, 'NRC':NRCopt, 'NRCerror':RCtotalerror}
    #returns a list of A,Aerr,B,Berr,M_RC,M_RCerr,sigma_RC,sigma_RCerr,N_RC,N_RCerr
    return output_dict



if __name__=="__main__":
    import numpy as np
    import pdb
    #pdb.set_trace()
    data = np.genfromtxt("../data/cmd_example.csv",delimiter=',')[1:]
    data = data[data[:,1]>1]
    data = data[np.logical_and(data[:,0]>12, data[:,0]<16)]
    data2 = RC_MCMC(data)
    print(data2)