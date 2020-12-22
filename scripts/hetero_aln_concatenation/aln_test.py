import os,sys
import subprocess
from string import Template,ascii_lowercase
import requests
import time

target_nf=[
    128,  # 64 from gremlin
    129, # 128 from baker casp12
]

id_cut=[
    99, # from metapsicov and deepcontact
    90, # from gremlin
]

cov_cut=[
    50, # from metapsicov 2.0.3 and deepcontact. for final alignment output
    60, # from metapsicov 1.04. for Nf calculation
    75, # from gremlin. for qhmmbuild searching
]

HHLIB = '/storage/htc/bdm/zhiye/DNCON4_db_tools/tools/deepmsa/hhsuite2/'
bin_dict=dict(
    #### upstream hhsuite executables ####
    hhblits =os.path.join(HHLIB,"bin/hhblits"),
    hhfilter      =os.path.join(HHLIB,"bin/hhfilter"),
    #### MSAParser ####
    rmRedundantSeq=os.path.join(HHLIB,"bin/rmRedundantSeq"),
    calNf     =os.path.join(HHLIB,"bin/calNf"),
)

hhblits_template=Template(bin_dict["hhblits"]+ \
" -i $infile -diff inf -d $db -cpu $ncpu -oa3m $outprefix.a3m "+ \
" -id $id_cut -cov $cov_cut -o $outprefix.log -n 3 -diff inf; "+ \
"grep -v '^>' $outprefix.a3m|sed 's/[a-z]//g' > $outprefix.aln")

hhfilter_template=Template(bin_dict["hhfilter"] + \
" -i $prefix.a3m -o $prefix.$cov_cut.a3m -id $id_cut -cov $cov_cut; "   + \
"grep -v '^>' $prefix.$cov_cut.a3m|sed 's/[a-z]//g' > $prefix.$cov_cut.aln")

alnfilter_template=Template(bin_dict["rmRedundantSeq"]+ \
" $id_cut $cov_cut $prefix.aln > $prefix.$cov_cut.aln")

calNf_template=Template(bin_dict["calNf"]+" $infile 0.8 0 $target_nf")

def getNf(prefix):
    ''' return Nf on -cov 60 filtered MSA. input file is prefix.a3m.
    output files are prefix.60.a3m and prefix.60.aln
    '''
    cov   =cov_cut[1]
    infile="%s.%d.aln"%(prefix,cov)

    #### filter input alignment ####
    if not os.path.isfile(infile):
        if os.path.isfile(prefix+".a3m"):
            cmd=hhfilter_template.substitute(
                prefix =prefix,
                id_cut =id_cut[0],
                cov_cut=cov,
            )
        else:
            cmd=alnfilter_template.substitute(
                prefix=prefix,
                id_cut=id_cut[0],
                cov_cut=cov,
            )
        sys.stdout.write(cmd+'\n')
        os.system(cmd)

    #### calculate Nf ####
    cmd=calNf_template.substitute(
        infile   =infile,
        target_nf=target_nf[-1],
    )
    stdout,stderr=subprocess.Popen(cmd,
        shell=True,stdout=subprocess.PIPE).communicate()
    return float(stdout)

def run_hhblits(query_fasta,db,ncpu,hhblits_prefix):
    ''' run hhblits with -cov 50 and return Nf for -cov 60'''
    cmd=hhblits_template.substitute(
        infile   =query_fasta,
        db   =db,
        ncpu     =ncpu,
        outprefix=hhblits_prefix,# outputs are $outprefix.a3m 
                                 # $outprefix.log $outprefix.aln
        id_cut   =id_cut[0],     # 99 in metapsicov
        cov_cut  =cov_cut[0],   # 50 in metapsicov 2.0.3
    )
    sys.stdout.write(cmd+'\n')
    os.system(cmd)
    return getNf(hhblits_prefix)

def get_stringID(proteinID):
    string_api_url = "https://string-db.org/api"
    output_format = "tsv-no-header"
    method = "get_string_ids"

    my_genes = [proteinID]

    request_url = "/".join([string_api_url, output_format, method])

    params = {
        "identifiers" : "%0d".join(my_genes), # your protein
    }
    response = requests.post(request_url, data=params)

    for line in response.text.strip().split("\n"):
        l = line.strip().split("\t")
        stringID = l[1]
        break # query one seq a time
    return stringID 

def get_interaction_partners(stringID, limit=10):
    string_api_url = "https://string-db.org/api"
    output_format = "tsv-no-header"
    method = "interaction_partners"

    my_genes = [stringID]

    request_url = "/".join([string_api_url, output_format, method])

    params = {
        "identifiers" : "%0d".join(my_genes), # your protein
        "limit" : limit,
    }
    response = requests.post(request_url, data=params)
    partners = []
    for line in response.text.strip().split("\n"):
        l = line.strip().split("\t")
        query_ensp = l[0]
        query_name = l[2]
        partner_ensp = l[1]
        partner_name = l[3]
        combined_score = l[5]
        ## print
        partners.append("\t".join([partner_ensp, combined_score]))
        # print("\t".join([partner_ensp, combined_score]))
    return partners

if __name__ == '__main__':
    query_fasta1 = '/storage/htc/bdm/zhiye/DeepComplex2/test/1A02F.fasta'
    query_fasta2 = '/storage/htc/bdm/zhiye/DeepComplex2/test/1A02J.fasta'

    db='/storage/htc/bdm/tools/multicom_db_tools/databases/uniprot20/uniprot20_2016_02/uniprot20_2016_02'
    ncpu=4
    print('run hhblits for 1A02F')
    # run_hhblits(query_fasta1,db,ncpu,'1A02F')
    print('run hhblits for 1A02J')
    # run_hhblits(query_fasta2,db,ncpu,'1A02J')

    query_name1 = '1A02F'
    query_name2 = '1A02J'

    # os.system('grep \'>\' %s.a3m > %s.seq'%(query_name1, query_name1))
    os.system('awk -F \'[||]\' \'{print $2}\' %s.a3m > %s.tmp'%(query_name1, query_name1))
    os.system('grep -v \'^$\' %s.tmp > %s.seq'%(query_name1, query_name1))
    
    time.sleep(t)
    #  ## python -m pip install requests

    # string_api_url = "https://string-db.org/api"
    # output_format = "tsv-no-header"
    # method = "interaction_partners"

    # my_genes = ["9606.ENSP00000000233", "9606.ENSP00000000412",
    #             "9606.ENSP00000000442", "9606.ENSP00000001008"]

    # request_url = "/".join([string_api_url, output_format, method])

    # params = {
    #     "identifiers" : "%0d".join(my_genes), # your protein
    #     "limit" : 5,
    # }
    # response = requests.post(request_url, data=params)

    # for line in response.text.strip().split("\n"):
    #     l = line.strip().split("\t")
    #     query_ensp = l[0]
    #     query_name = l[2]
    #     partner_ensp = l[1]
    #     partner_name = l[3]
    #     combined_score = l[5]
    #     ## print
    #     print("\t".join([query_ensp, query_name, partner_name, combined_score]))