#!!usr/bin/env python
import pdb,sys,os
import warnings
warnings.filterwarnings("ignore")
from File import *
from BioUtils import BioList
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import networkx
import Bio.Phylo
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor,DistanceMatrix
from scipy.spatial.distance import cityblock
import random
from BioUtils import BioList
import pyBigWig as pybw
import argparse
from functools import reduce

class Group:
    def __init__(self,cells,types=None,features=None):
        self.cells=cells
        self.types=types
        self.features=features
        if self.types!=None:
            self._typestats=sorted([[item,self.types.count(item),round(self.types.count(item)*1.0/len(self.types),2)] for item in set(self.types)],key=lambda x:x[1],reverse=True)

def List2String(X):
    out=[]
    for i in X:
            out.append(",".join([str(item) for item in i]))
    out=";".join(out)
    return out


# get the mutation matrix, each row represents all mutations found in the cell
def buildMutationMatrix(dcell2snp,keptMutations):
    matrix=[]
    for cell in dcell2snp:
        mcell=dcell2snp[cell]
        row=[]
        for j in keptMutations:
            if j in mcell:
                    row.append(1)
            else:
                    row.append(0)
        row=[cell]+row
        matrix.append(row)
    return matrix

# build groups based on clustering
def buildGroups(matrixPreds,matrixCells,matrixValues,dsra2type):
    groups=[]
    for i in set(matrixPreds):
        ci=[]
        di=[]
        for j in range(len(matrixPreds)):
            if i==matrixPreds[j]:
                ci.append(matrixCells[j])
                di.append(matrixValues[j])
        if len(dsra2type)>0:
            ti=[dsra2type[item] for item in ci]
            gi=Group(ci,ti,di)
        else:
            gi=Group(ci,None,di)
        groups.append(gi)
    return groups

# build the distance for building phylogenetic tree
def buildPhyloDM(groups):
    names=[]
    ct=0
    pcut=0.2
    #pcut=0.0
    for i in groups:
        if i.types!=None:
            di=[[item,i.types.count(item)] for item in set(i.types) \
                    if i.types.count(item)*1.0/len(i.types)>pcut]
            di=sorted(di,key=lambda x: x[1],reverse=True)
            strdi=str(ct)+"|"+List2String(di)
            #pdb.set_trace()
        else:
            strdi=str(ct)
        names.append(strdi)
        ct+=1
    matrix=[]
    for i in range(len(groups)):
        irow=[]
        for j in range(i+1):
            mij=distanceij(groups[i],groups[j]) if i!=j else 0
            irow.append(mij)
        matrix.append(irow)
    dm=DistanceMatrix(names,matrix)
    return dm

# build the phylogenetic tree based on calculated distance
def buildPhylogenicTree(dm):
    TreeConstructor=DistanceTreeConstructor()
    ptree=TreeConstructor.nj(dm)
    return ptree

# distance between two groups
def distanceij(gi,gj):
    gif=gi.features
    gjf=gj.features
    ds=0
    ct=0
    for p in gif:
        for q in gjf:
            dpg=cityblock(p,q)
            ds+=dpg
            ct+=1
    ds=ds*1.0/ct
    return ds


def drawTree(net,outstr,outputdir):
    pos=networkx.spring_layout(net,seed=0)
    #Xlabels={net.nodes.keys()[k]:net.nodes.keys()[k].name.split("|")[0] for k in range(len(net.nodes))}
    Xlabels={}
    for k in range(len(net.nodes)):
        nkeys=list(net.nodes.keys())
        Xlabels[nkeys[k]]=nkeys[k].name.split("|")[0]
        
    XlabelTexts=[item.name.split("|") for item in Xlabels]
    Xnames=[item[0] for item in XlabelTexts]
    XlabelTexts=[item[-1] if len(item)>1 else "" for item in XlabelTexts]
    networkx.draw(net,pos,labels=Xlabels,with_labels=True)
   
    datgraph=[]
    for i in range(len(pos.values())):
        [x,y]=list(pos.values())[i]
        datgraph.append([Xnames[i],(x,y)])
        plt.text(x-0.003*len(XlabelTexts[i]),y+0.05,s=XlabelTexts[i],fontsize=5)
    plt.subplots_adjust(right=0.5)
    plt.savefig("%s/%s.png"%(outputdir,outstr),bbox_inches="tight",dpi=600)
    BioList(datgraph).ex2File("%s/%s.dat"%(outputdir,outstr))

# check whether new feature signature list is better than the old one, if yes -> update the feature list
def evalChange(xold,eNew,dfsnp,icells,ccells):
    mold=getMutationX(xold,icells,ccells,dfsnp)
    mnewList=[]
    # try 2 potential operations: + or |
    for o in ["|","+"]:
            xnew= [o+eNew]+xold if o=="|" else  xold+[o+eNew]
            mnew=getMutationX(xnew,icells,ccells,dfsnp)
            mnewList.append([xnew,mnew])

    maxNew=max(mnewList,key=lambda x: x[1]-mold)
    if maxNew[1]>mold:
            return maxNew[0]
    return xold

# evalulate the updated mutation list
def evalUpdateMut(updatedMutations):
    raiseException("Not implemented")


# use greedy algorithm to find the best mutation makers for given cells (compared with a group of controll cells)
def findMutationMarker(icells,ccells,dfsnp,pcut):
    iMutations={}
    cMutations={}
    for i in dfsnp:
            imut=dfsnp[i]
            iMutations[i]=[item for item in imut if item in icells]
            cMutations[i]=[item for item in imut if item in ccells]

    rankedMutationList=[]
    for i in iMutations:
            mi=len(iMutations[i])
            mc=0 if i not in cMutations else len(cMutations[i])
            rankedMutationList.append([i,mi*1.0/len(icells)-mc*1.0/len(ccells)])
    rankedMutationList=sorted(rankedMutationList,key=lambda x: x[1],reverse=True)

    # greedy algorithm  to find the mutation marker
    feMut=[]  # feature mutation list
    for i in rankedMutationList:
            xi=evalChange(feMut,i[0],dfsnp,icells,ccells)
            feMut=xi

    ffmut=[]
    mutscore=getMutationX(feMut,icells,ccells,dfsnp)
    for i in feMut:
        ffmut.append(i)
        if getMutationX(ffmut,icells,ccells,dfsnp)/mutscore>pcut:
            break

    return ffmut


# get the score of given mutation, the score is showing whether the mutation is found significantly more in icells vs ccells
def getMutationX(x,icells,ccells,dfsnp):
    xmut=[]
    for i in x:
            iop=i[0]
            ikey=i[1:]
            mikey=dfsnp[ikey]
            xmut.append([iop,ikey])

    xcells=[]
    for i in xmut:
            if i[0]=="|":
                    xcells+=dfsnp[i[1]]
            else:
                    xcells=[item for item in xcells if item in dfsnp[i[1]]]

    xcells=list(set(xcells))
    xicells=[item for item in icells if item in xcells]
    xccells=[item for item in ccells if item in xcells]

    score=len(xicells)*1.0/len(icells)-len(xccells)*1.0/len(ccells)
    return score


# get cells using mutation marker
def getCellsByMK(x,dfsnp):
    xmut=[]
    for i in x:
            iop=i[0]
            ikey=i[1:]
            mikey=dfsnp[ikey]
            xmut.append([iop,ikey])

    xcells=[]
    for i in xmut:
            if i[0]=="|":
                    xcells+=dfsnp[i[1]]
            else:
                    xcells=[item for item in xcells if item in dfsnp[i[1]]]

    xcells=list(set(xcells))
    return xcells


# get the dictionary : cell->snp
# filter the snp with low reads coverage (potential non-expressing snps)
def getcell2snp(dsnp,cutl,cuth,cells,dbw):
    umuts=sorted([item for item in dsnp if isUsefulMutation(item,dsnp,cutl,cuth,cells)])
    # cell-> snp
    dcell2snp={}
    for i in umuts:
        icells=dsnp[i]
        for j in icells:
            if j not in dcell2snp:
                dcell2snp[j]=[i]
            else:
                dcell2snp[j].append(i)
   
    # filtered snp (filter the snp with low reads converage)
    keptMutations=[] 
    dsnp2reads={} # snps->reads
    for i in range(len(umuts)):
        if len(dbw)>0:
            [iflag,ireads]=isExpressing(umuts[i],cells,dbw)
            dsnp2reads[umuts[i]]=ireads
        else:
            iflag=1
        #----------------
        if iflag:
            keptMutations.append(umuts[i])
        print(i) 
    #pdb.set_trace()
    dfsnp={item:dsnp[item] for item in keptMutations}
    return [dfsnp,dcell2snp,umuts,keptMutations,dsnp2reads]

# get the signature mutation list for each of the groups
def getGroupSignatureMutation(groups,allcells,dfsnp,pcut):
    selmuts=[]
    for g in groups:
        icells=g.cells
        ccells=[item for item in allcells if item not in icells]
        mk=findMutationMarker(icells,ccells,dfsnp,pcut)
        g.mk=mk
        for i in mk:
            imk=i[1:]
            if imk not in selmuts:
                selmuts.append(imk)
    return selmuts

# get all mutation list (not necessarily the signature ones) for each of the groups
def getGroupAllMutation(groups,selmuts,dfsnp):
    cut=0.8
    dg={}
    for m in selmuts:
        mcells=dfsnp[m]
        for i in range(len(groups)):
            g=groups[i]
            gcells=g.cells
            mgcells=[item for item in mcells if item in gcells]
            fc=len(mgcells)*1.0/len(gcells)
            if fc>cut:
                if i not in dg:
                    dg[i]=[m]
                else:
                    dg[i].append(m)
    for i in dg:
        dg[i].sort()
    return dg

def getMutationDistribution(x,dfsnp,groups):
    xcells=dfsnp[x]
    xg=[]
    for g in groups:
        xgcells=[item for item in xcells if item in g.cells]
        xg.append(len(xgcells)*1.0/len(g.cells))
    return xg

def getSNPReadsDistribution(x,cells,dsnp2reads,groups):
    #TODO implement the function
    xreads=dsnp2reads[x]

    xg=[]
    for g in groups:
        gcellindexList=[cells.index(item) for item in g.cells]
        xgreads=[xreads[item] for item in gcellindexList]
        xg.append(sum(xgreads)/len(xgreads))
    return xg

# remove non-expressing SNPs
def isExpressing(i,icells,dbw,w=50,rcut=8,pcut=0.8):
    [ichr,ipos]=i.split(",")
    istart=max(0,int(ipos)-w/2)
    iend=int(ipos)+w/2
    
    rcells=[]
    for icell in icells:
        ri=dbw[icell].stats(ichr,int(istart),int(iend))[0]
        rcells.append(ri)

    rfcells=[item for item in rcells if item>rcut]

    flag=len(rfcells)*1.0/len(icells)>pcut
    return [flag,rcells]

# remove unuseful mutations
def isUsefulMutation(i,dsnp,cutl,cuth,allcells):
    icells=dsnp[i]
    res=len(icells)*1.0/len(allcells)
    if res>cuth:
            return False
    if res<cutl:
            return False
    return True

# perform clustering
def performClustering(dcell2snp,keptMutations,dsra2type,k=7):
    k=7 if k==None else k
    matrix=buildMutationMatrix(dcell2snp,keptMutations)
    matrixCells=[item[0] for item in matrix]
    matrixValues=[[float(j) for j in item[1:]] for item in matrix]
    kk=KMeans(n_clusters=k,random_state=0)
    matrixPreds=kk.fit_predict(matrixValues)
    groups=buildGroups(matrixPreds,matrixCells,matrixValues,dsra2type)
    ss=silhouette_score(matrixValues,matrixPreds,random_state=0)  # silhouette score
    return [groups,ss]

# update mutation list by using the keptMutations and randomly choosing other mutations
def updateMutList(dcell2snp,keptMutations,dsra2type,dfsnp,k):
    unchosenMut=[item for item in dfsnp if item not in keptMutations]
    rs=[]
    for m in unchosenMut:
        updatedMuts=keptMutations+[m]
        [ug,us]=performClustering(dcell2snp,updatedMuts,dsra2type,k)
        rs.append([us,updatedMuts])
    rs=sorted(rs,reverse=True)
    updatedMuts=rs[0][1]
    return updatedMuts

def exportGroupCells(groups,outputdir,dsra2type):
    out=[]
    for i in range(len(groups)):
        icells=groups[i].cells
        for j in icells:
            ki=[j,"cluster:%s"%(i)] if j not in dsra2type else [j,"cluster:%s"%(i),dsra2type[j]]
            out.append(ki)
    BioList(out).ex2File("%s/GroupCells.txt"%(outputdir))


def exportSNP(dfsnp,groups,keptMutations,outputdir):
    #mutmat[getMutationDistribution(item,dfsnp,groups) for item in keptMutations]
    mutmat=[]
    for i in keptMutations:
        icells=dfsnp[i]
        imut=[]
        for g in groups:
            gcells=g.cells
            gmut=[1 if item in icells else 0 for item in gcells]
            imut+=gmut
        mutmat.append(imut)
    #-------------------------------------------
    fig,axs=plt.subplots(1,len(groups),sharey=True,gridspec_kw={'width_ratios':[len(item.cells) for item in groups]})
    st=0
    for i in range(len(groups)):
        ed=st+len(groups[i].cells)
        axs[i].imshow([item[st:ed] for item in mutmat],cmap="hot",aspect="auto",interpolation="nearest")
        st=ed

    for ax,l in zip(axs,range(len(groups))):
        ax.set_xticklabels([])
        ax.set_xlabel(l)   
    
    fig.add_subplot(111,frameon=False)
    plt.tick_params(labelcolor='none',top='off',bottom='off',left='off',right='off')
    plt.grid(False)
    plt.ylabel("SNPs")
    plt.xlabel("Clusters")
    
    plt.savefig("%s/SNP_matrix.pdf"%(outputdir),dpi=600)
    
    
    mutmatout=[]

    headers=['SNP']+reduce(lambda x,y:x+y,[item.cells for item in groups])
    mutmatout.append(headers)
    for i in range(len(keptMutations)):
        mutmat[i]=[keptMutations[i]]+mutmat[i]
        mutmatout.append(mutmat[i])

    BioList(mutmatout).ex2File("%s/SNP_matrix.tsv"%(outputdir))
    

def bestK(dcell2snp,keptMutations,dsra2type):
    matrix=buildMutationMatrix(dcell2snp,keptMutations)
    matrixCells=[item[0] for item in matrix]
    matrixValues=[[float(j) for j in item[1:]] for item in matrix]
    ssl=[]
    K=range(2,11)
    for k in K: 
        kk=KMeans(n_clusters=k,random_state=0)
        matrixPreds=kk.fit_predict(matrixValues)
        groups=buildGroups(matrixPreds,matrixCells,matrixValues,dsra2type)
        ss=silhouette_score(matrixValues,matrixPreds,random_state=0)  # silhouette score
        ssl.append(ss)
    ssl=[item/sum(ssl) for item in ssl]
    pks=detPeak(K,ssl,0.05,type='max')
    pks=[item for item in pks if item>2]
    k=pks[0] if len(pks)>0 else None
    return k

def detPeak(K,A,deltaP,type='max'):
    LM=[]
    LX=[]
    if len(A)==0:
        return []
    minA=min(0,min(A))
    for i in range(len(A)-1):
        delta=deltaP* (A[i]-minA)
        if i==0:
            ad=A[i+1]-A[i]
            if ad > delta:
                LM.append(i)
            if ad<=1*delta:
                LX.append(i)
        else:
            pd=A[i]-A[i-1]
            ad=A[i+1]-A[i]
            if (pd>delta) and (ad <-0.25* delta):
                LX.append(i)
            if (pd<-1* delta) and (ad > 0.25*delta):
                LM.append(i)
    if type=='min':
        LM=sorted(LM,key=lambda x: A[x])
        return [K[item] for item in LM]
    LX=sorted(LX,key=lambda x: A[x], reverse=True)
    

    return [K[item] for item in LX]

# main program starts here
def main():
    # parse the arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--ivcf",help="Required,directory with all input .vcf files. This specifies the directory of SNP files (.vcf) for the cells (one .vcf file for each cell). These .vcf files can be obtained using the \
            provided bam2vcf script or other RNA-seq variant calling pipelines preferred by the users.",required=True)
    parser.add_argument("-b","--ibw",help="Optional,directory with all input bigwig (.bw) files with the information about the number of aligned reads at each genomic position. These bigwig files are used to filter the SNPs, \
            which are redundant to expression information.",nargs="?")
    parser.add_argument("-k","--kcluster",help="Optional, number of clusters, Integer. If not specified, the program will choose the k with best silhouette score.", default='auto')
    parser.add_argument("-l","--cell_label", help="Optional, labels for the cells. This is used only to annotate the cells with known information, not used for building the model.",nargs="?")
    parser.add_argument("-o","--output",help="Required,output directory",required=True)
    parser.add_argument("--cutl",help="Optional, lower bound cutoff to remove potential false positive SNPs, default=0.1",default=0.1)
    parser.add_argument("--cuth",help="Optional, upper bound cutoff to remove baseline SNPs, which are common in most cells, default=0.8",default=0.8)
    parser.add_argument("--greedycut",help="Optional, the stopping cutoff for the greedy search of candidate SNPs, default=0.05 (less than 0.05 score improvement). A larger cutoff means less strict SNP candidate search", default=0.05)
    parser.add_argument("--cutc",help="Optional, convergence cutoff, a smaller cutoff represents a stricter convergence criterion,default=0.001",default=0.001)
    parser.add_argument("--maxiter",help="Optional, the maximal number of iterations allowed",default=10)
    
    args=parser.parse_args()
	
    vcfolder=args.ivcf
    bwfolder=args.ibw
    sralabelfile=args.cell_label
    outputdir=args.output
    try:
        cutl=float(args.cutl)
        cuth=float(args.cuth)
        loopcut=float(args.cutc)
        Loop=float(args.maxiter)
        gcut=1-float(args.greedycut)
    except:
        print("please check your input! Only float numbers are allowed for cutoffs")
        sys.exit(0)
    k=args.kcluster
    
    #--------------------------------------------
    if os.path.exists(outputdir)==False:
        os.mkdir(outputdir)

    # L: the result folder of all vcf files
    L=os.listdir(vcfolder)
    cells=[item.split(".")[0] for item in L]

    dsnp={}
    ct=0
    for i in L:
        isnp=TabFile("%s/%s"%(vcfolder,i)).read("\t")
        isnp=[item for item in isnp if len(item)>1]
        for j in isnp:
            key=j[0]+","+j[1]
            if key not in dsnp:
                dsnp[key]=[i.split(".")[0]]
            else:
                dsnp[key].append(i.split(".")[0])
        ct+=1
        print("cell:%s"%(ct))

    # read in bw files
    dbw={}
    if bwfolder!=None:
        bwFileList=os.listdir(bwfolder)
        bwList=[]
        for i in bwFileList:
            bi=pybw.open('%s/%s'%(bwfolder,i))
            si=i.split(".")[0]
            dbw[si]=bi
    
    # the labels for each cell
    dsra2type={}

    if sralabelfile!=None:
        sra2labels=TabFile(sralabelfile).read("\t")

        for i in sra2labels:
            dsra2type[i[0]]=i[1]

    
    #============================================================================
    [dfsnp,dcell2snp,umuts,keptMutations,dsnp2reads]=getcell2snp(dsnp,cutl,cuth,cells,dbw)
   
    #---------------------------------------------
    if k=='auto':
        k=bestK(dcell2snp,keptMutations,dsra2type)
    else:
        try:
            k=int(k)
        except:
            print("please check your -k parameter")
            sys.exit(0)
    #--------------------------------------------

    #pdb.set_trace()

    Loop=10
    ssList=[]
     
    #FIXME: update the main function
    print('clustering...')
    
    # initial clustering...
    [groups,ss]=performClustering(dcell2snp,umuts,dsra2type,k)
    ssList.append(ss)

    # find informative snps based on current groups
    keptMutations=getGroupSignatureMutation(groups,cells,dfsnp,gcut)
    [groups,ss]=performClustering(dcell2snp,keptMutations,dsra2type,k)
    ssList.append(ss)
    
    for l in range(Loop):
        print("loop:%s"%(l))
        updatedMuts=updateMutList(dcell2snp,keptMutations,dsra2type,dfsnp,k)
        [ugroups,uss]=performClustering(dcell2snp,updatedMuts,dsra2type,k)
        if (uss>ss):
            keptMutations=updatedMuts
            [groups,ss]=performClustering(dcell2snp,keptMutations,dsra2type,k)

        if abs(ss-ssList[-1])<loopcut:
            ssList.append(ss) 
            break
        ssList.append(ss)
    #pdb.set_trace()
    dm=buildPhyloDM(groups)
    ptree=buildPhylogenicTree(dm)
    net=Bio.Phylo.to_networkx(ptree)
    drawTree(net,"Trajectory",outputdir)
    
    # exporting snps and snp reads
    exportSNP(dfsnp,groups,keptMutations,outputdir)
    exportGroupCells(groups,outputdir,dsra2type)

if __name__=="__main__":
    main()
