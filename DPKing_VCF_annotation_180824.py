import datetime
import requests
import sys
import json 
#from pprint import pprint #used to view json data while developing


dt=str(datetime.datetime.now())
dts = ''.join(filter(str.isdigit, dt))
dtstr = (dts[2:8]+"_"+dts[8:12])
print(dtstr)

f = open("Challenge_data (1).vcf", "r")
r = open("VCF file annotations_"+dtstr+".txt", "w")

print("Annotating...")

#The lists below are used to accommodate the occurence of
#multiple alternate alleles at a given locus (i.e., VCF file row)
M_ALT     = []
M_TYPE    = []
M_ALTDPTH = []
x=0 #counter for lists above

z=0 #VCF row counter

tabs = []
tabindex = -1

#The sets below are used to generate non-redundant lists from the ExAC database 
EXACvartype_set    = set()
allEXACvartype_set = set()
EXACgene_set       = set()
VARrank = [] #for ranking severity of ExAC variants (details below)

version = str(sys.argv[0][:-3])

r.write("Script: "+version[version.find("DPK"):]+"\n\n")
#annotation data column headers:
r.write("{}\t{:4}\t{:10}\t{:7}\t{:7}\t{:7}\t{:7}\t{:>6}\t{:>6}\t{:7}\t{:23}\t{:15}\t{}\n"
        .format("Var#","Chr","Pos","VCFvar","Depth","R-reads","V-reads","%Var","%Ref","ExACfrq","ExACvar","ExACgene","Chr-Pos-R-A ID"))

for row in f:
    if row[0]=="#":
        continue 

    z+=1
    if z%25 == 0:
        print("row",z)

    #if z > 30:
        #break
    
    while True: 
        tabindex = row.find("\t",tabindex+1)
        if tabindex == -1:
            break
        tabs.append(tabindex)
        
    #VCF row parsing for relevant details:
    CHR    =  (row[0:tabs[0]])
    POS    =  (row[tabs[0]+1:tabs[1]])
    REF    =  (row[tabs[2]+1:tabs[3]])
    preALT =  (row[tabs[3]+1:tabs[4]])
    M_ALT  =   preALT.split(",") #to capture multiple alternate alleles, when present

    TYPEindex   =   row.find(";TYPE=")
    preTYPE     =  (row[TYPEindex+6:(row.find("\t",TYPEindex+1))])
    M_TYPE      =   preTYPE.split(",") #to capture multiple alternate allele types, when present
    DPTHindx    =   row.find(";DP=")
    DPTH        =  (row[DPTHindx+4:(row.find(";",DPTHindx+1))])
    REFDPTHindx =   row.find(";RO=")
    REFDPTH     =  (row[REFDPTHindx+4:(row.find(";",REFDPTHindx+1))])
    ALTDPTHindx =   row.find(";AO=")
    ALTDPTHstr  =  (row[ALTDPTHindx+4:(row.find(";",ALTDPTHindx+1))])
    M_ALTDPTH   =   ALTDPTHstr.split(",") #to capture multiple alternate allele read depths, when present
    
    '''
    #the lines below may be used to sum all alternate reads associated with a given reference locus.
    ALTDPTH = 0    
    for item in M_ALTDPTH:
        ALTDPTH+=int(item)
    '''
    
    NUMALTindx = row.find(";NUMALT=")  #to identify the number of alternate alleles for each reference locus
    NUMALT     = int(row[NUMALTindx+8:(row.find(";",NUMALTindx+1))])

    x=0   
    while x < NUMALT:

        EXACid       =  CHR+"-"+POS+"-"+REF+"-"+M_ALT[x] 
        EXACdatapull =  requests.get("http://exac.hms.harvard.edu/rest/variant/"+EXACid)
        EXACdata     =  EXACdatapull.json()
        
        if "vep_annotations" in EXACdata["variant"]:
            for item in EXACdata["variant"]["vep_annotations"]:
                EXACvartype = item["major_consequence"] 
                EXACvartype_set.add(EXACvartype)
                #allEXACvartype_set.add(EXACvartype) #this set was used while developing the program
                                                     #to identify all ExAC variant types in the VCF file
                                                     #in order to generate the VARrank list below.
                EXACgene = item["SYMBOL"]
                if item["SYMBOL_SOURCE"] == "HGNC":  
                    EXACgene_set.add(EXACgene)
        else:
            EXACgene_set.add("N/A")
            
        #Rank order of severity of variant types is adapted from Daniel MacArthur's lab -
        #https://github.com/macarthur-lab/seqr/blob/master/xbrowse/annotation/vep_annotations.py
        VARrank = ["splice_acceptor_variant",
                    "splice_donor_variant",
                    "stop_gained",
                    "initiator_codon_variant",
                    "stop_lost",
                    "missense_variant",
                    "splice_region_variant",
                    "synonymous_variant",
                    "stop_retained_variant",
                    "3_prime_UTR_variant",
                    "5_prime_UTR_variant",
                    "intron_variant",
                    "non_coding_transcript_exon_variant"]

        for item in VARrank:
            if item in EXACvartype_set:
                EXACvartype = item
                break
            else:
                EXACvartype = "N/A"
                        
        if "allele_freq" in EXACdata["variant"]:
            EXACfreq = "%7.5f" % float(EXACdata["variant"]["allele_freq"]) 
        else:
            EXACfreq = "N/A"

        #writing to annotation file. Generally, I anticipate taking advantage of excel, etc., to examine
        #an output file such as this one, but I have also applied some formatting here, for improved readability.
        #(Note that some content in the 'ExACgene_set' column overflows its designated width.)
            
        if NUMALT == 1:
            r.write(str(z)+"\t")
        else:
            r.write(str(z)+"."+str(x+1)+"\t")
            
        r.write("{:4}\t{:10}\t{:7}\t{:7}\t{:7}\t{:7}\t".format( CHR, POS, M_TYPE[x], DPTH, REFDPTH, M_ALTDPTH[x] ))
        r.write("{:>6.1%} \t{:>6.1%} \t".format( int(M_ALTDPTH[x])/int(DPTH), int(REFDPTH)/int(DPTH) ))
        r.write("{:7}\t{:23}\t{:15}\t{}\t".format( EXACfreq, EXACvartype, ",".join(EXACgene_set), EXACid ))
        r.write("\n")

        x+=1
        tabs.clear()
        EXACvartype_set.clear()
        EXACgene_set.clear()

print("Done!")
print("VCF rows annotated:",z-1)    

f.close()
r.close()

