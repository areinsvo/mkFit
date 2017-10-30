import os.path, glob, sys
import ROOT
import array
import math

def run():
    arch   = sys.argv[1] # SNB, KNC, KNL
    sample = sys.argv[2] 

    g = ROOT.TFile('benchmark_'+arch+'_'+sample+'.root','recreate')

    for build in ['BH','STD','CE']:
        print arch,sample,build

        # Vectorization data points
        vuvals = ['1','2','4','8']
        nth = '1'

        if arch == 'KNC' or arch == 'KNL' :
            vuvals.append('16')
            vuvals.append('16int')
        else : 
            vuvals.append('8int')

        # call the make plots function
        makeplots(arch,sample,build,vuvals,nth,'VU')

        # Parallelization datapoints
        if arch == 'KNC' or arch == 'KNL' :
            nvu = '16int'
            if arch == 'KNC' :
                thvals = ['1','2','4','8','15','30','60','90','120','150','180','210','240']
            else : # KNL
                thvals = ['1','2','4','8','16','32','64','96','128','160','192','224','256']
        else : # SNB
            nvu = '8int'
            thvals = ['1','2','4','6','8','12','16','20','24']
    
        # call the make plots function
        makeplots(arch,sample,build,thvals,nvu,'TH')

    g.Write()
    g.Close()

def makeplots(arch,sample,build,vals,nC,text):

    # position in logs
    if   build is 'BH'  : pos = 8  
    elif build is 'STD' : pos = 11  
    elif build is 'CE'  : pos = 14 
    else :
        print build,'is not a valid test! Exiting...'
        exit 

    # time    
    print arch,sample,build,text

    # define tgraphs vs absolute time and speedup
    g_time    = ROOT.TGraphErrors(len(vals)-1)
    g_speedup = ROOT.TGraphErrors(len(vals)-1)

    # make separate plot for intrinsics measurement
    if text is 'VU' :
        g_time_int    = ROOT.TGraphErrors(1)
        g_speedup_int = ROOT.TGraphErrors(1)

    point = 0
    for val in vals :
        if    val is '16int': xval = 16.0
        elif  val is '8int' : xval = 8.0
        else                : xval = float(val)

        # array of time values
        yvals = array.array('d');

        # always skip the first event
        firstFound = False

        # open the correct log file, store times into temp file
        if   text is 'VU' : os.system('grep Matriplex log_'+arch+'_'+sample+'_'+build+'_NVU'+val+'_NTH'+nC +'.txt >& log_'+arch+'_'+sample+'_'+build+'_'+text+'.txt')
        elif text is 'TH' : os.system('grep Matriplex log_'+arch+'_'+sample+'_'+build+'_NVU'+nC +'_NTH'+val+'.txt >& log_'+arch+'_'+sample+'_'+build+'_'+text+'.txt')
        else :
            print 'VU or TH are the only options for extra text! Exiting...'
            exit

        # open temp file, store event times into yvals
        with open('log_'+arch+'_'+sample+'_'+build+'_'+text+'.txt') as f :
            for line in f :
                if 'Matriplex' not in line : continue
                if 'Total' in line : continue
                if not firstFound :
                    firstFound = True
                    continue
                lsplit = line.split()
                yvals.append(float(lsplit[pos]))

        # Compute mean and uncertainty on mean from yvals
        sum = 0.;
        for yval in range(0,len(yvals)):
            sum = sum + yvals[yval]
        mean = sum/len(yvals)
        emean = 0.;
        for yval in range(0,len(yvals)):
            emean = emean + ((yvals[yval] - mean) * (yvals[yval] - mean))
        emean = math.sqrt(emean / (len(yvals) - 1))
        emean = emean/math.sqrt(len(yvals))

        # Printout value for good measure
        print val,mean,'+/-',emean

        # store intrinsics val into separate plot
        if 'int' not in val :
            g_time.SetPoint(point,xval,mean)
            g_time.SetPointError(point,0,emean)
            point = point+1
        else :
            g_time_int.SetPoint(0,xval,mean)
            g_time_int.SetPointError(0,0,emean)

    # always write out the standard plot
    g_time.Write('g_'+build+'_'+text+'_time')

    # write out separate intrinsics plot
    if text is 'VU' :
        g_time_int.Write('g_'+build+'_'+text+'_time_int')

    # Speedup calculation
    xval0 = array.array('d',[0])
    yval0 = array.array('d',[0])
    yerr0 = array.array('d',[0])

    # Get first point to divide by
    g_time.GetPoint(0,xval0,yval0)
    yerr0.append(g_time.GetErrorY(0))

    point = 0
    for val in vals :
        # set up inputs
        xval = array.array('d',[0])
        yval = array.array('d',[0])
        yerr = array.array('d',[0])

        # get standard plots from standard plot
        if 'int' not in val :
            g_time.GetPoint(point,xval,yval)
            yerr.append(g_time.GetErrorY(point))
        else :
            g_time_int.GetPoint(0,xval,yval)
            yerr.append(g_time_int.GetErrorY(0))

        speedup  = 0.
        espeedup = 0.
        if yval[0] > 0. : 
            speedup  = yval0[0]/yval[0]
            espeedup = speedup * math.sqrt(math.pow(yerr0[0]/yval0[0],2) + math.pow(yerr[0]/yval[0],2))

        # store in the correct plot
        if 'int' not in val :
            g_speedup.SetPoint(point,xval[0],speedup)
            g_speedup.SetPointError(point,0,espeedup)
            point = point+1
        else :
            g_speedup_int.SetPoint(0,xval[0],speedup)
            g_speedup_int.SetPointError(0,0,espeedup)

    # always write out the standard plot
    g_speedup.Write('g_'+build+'_'+text+'_speedup')

    # write out separate intrinsics plot
    if text is 'VU' :
        g_speedup_int.Write('g_'+build+'_'+text+'_speedup_int')

    # all done
    return

if __name__ == "__main__":
    run()
