#!/bin/env python3
# Author: dfguan9
# Contact: dfguan9@gmail.com
# Insitute: Zoology Institute, CAS
# Date: 2020-08-14
# Function: parse an assembly summary json file to obtain the assembly information and print to a file in table delimited format 

import sys, json, argparse, requests, time
import re

def load_tax_tree(treepath):
    tax_tree = {}
    with open(treepath) as fl:
        for line in fl:
            llst = line.strip().split('\t|\t') 
            if llst[3] != '':
                tax_tree[llst[0]] = [llst[6], llst[5],llst[4]]
    return tax_tree 

def get_sN50_others(asmid):
    asmurl = "https://www.ncbi.nlm.nih.gov/assembly/{}".format(asmid)
    time.sleep(5)
    a = ["NA", "NA", "NA", "NA", "NA"]
    try:
        res = requests.get(asmurl)
    except: 
        return a
        # print (res.text)
    m = re.search("<td>Scaffold N50</td><td class=\"align_r\">(.+?)</td>", res.text)
    if m is not None:
        a[0] = m.group(1)
    m = re.search("<dt>Submitter: </dt><dd>(.+?)</dd>", res.text)
    if m is not None:
        a[1] = m.group(1)
    m = re.search("<dt>Sequencing technology: </dt><dd>(.+?)</dd>", res.text)
    if m is not None:
        a[2] = m.group(1)
    m = re.search("<dt>BioProject: </dt><dd><a href=(.+?)>(.+?)</a></dd><dt>", res.text)
    if m is not None:
        a[3] = m.group(2)
    m = re.search("<dt>BioSample: </dt><dd><a href=(.+?)>(.+?)</a></dd>", res.text)
    if m is not None:
        a[4] = m.group(2)
    # print (a)
    print ("Extracted detailed assembly information from {}".format(asmurl), file=sys.stderr)
    return a

def get_ptname(taxid):
    taxurl = "https://api.ncbi.nlm.nih.gov/datasets/v1alpha/genome/taxon/{}/tree".format(taxid)
    res = requests.get(taxurl)
    res_dict = json.loads(res.text)
    rv = []
    # print (res_dict)
    if res_dict is not None:
        if 'rank' in res_dict:
            rv.append(res_dict['rank'])
        else:
            rv.append('NA')
        if 'sci_name' in res_dict:
            rv.append(res_dict['sci_name'])
        else:
            rv.append('NA')
        if 'parent_tax_id' in res_dict:
            rv.append(res_dict['parent_tax_id'])
        else:
            rv.append('NA')
    return rv # [rank, sci_name, ptaxid]


def parse_asminfo2(url, tax_tree, nrec):

    res = requests.get(url)
    if res is not None:
        print ("Get taxonomy list, start analyzing...", file=sys.stderr)
    # with open(res.json) as f:
        res_dict = json.loads(res.text)
        asm_sumy = None
        if 'datasets' in res_dict:
            asm_sumy = res_dict['datasets']
        else:
            print ("Error: No assembly is available, please check species name or taxonomy ID")
            return 1

        # {'assembly_accession': 'GCA_000005465.1', 'display_name': 'BGIAF', 'org': {'tax_id': '9606', 'sci_name': 'Homo sapiens', 'common_name': 'human', 'sex': 'male', 'rank': 'SPECIES', 'parent_tax_id': '9605', 'assembly_counts': {'node': 128, 'subtree': 128}, 'key': '9606', 'title': 'human'}, 'chromosomes': ['Un'], 'assembly_level': 'Scaffold', 'submission_date': '2010-06-21', 'contig_n50': 887, 'estimated_size': '652893711', 'seq_length': '2676008911'}
        # print header
        hdr = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}\t{15}".format("Class", "Order", "Family","Accession ID",  "Latin name", "Common name", "Assembly level", "Submission date", "Contig N50", "Genome size", "Annotated", "Scaffold N50", "Submitter", "Sequencing technology", "Project ID", "Sample ID")
        print (hdr)
        n_records = len(asm_sumy)
        onep_step = int(0.01 * n_records)
        if onep_step == 0:
            onep_step += 1
        if nrec == -1:
            nrec = n_records
        elif nrec == 0:
            print ("No records are required, exiting")
            return 0
        
        print ("There are {0} records in total, downloading {1} of them".format(n_records, nrec), file=sys.stderr)
        counter = 0
        for ele in asm_sumy:
            if counter >= nrec:
                break   
            attr = []
            if "org" in ele and "tax_id" in ele["org"]:
                tid = ele["org"]["tax_id"]
                if tid in tax_tree:
                    attr.extend(tax_tree[tid])
                else:
                    attr = ["NA", "NA", "NA"]
            else:
                attr = ["NA", "NA", "NA"]
            # get taxonomy ranks 
            # if "org" in ele:
                # if "rank" in ele["org"]:
                    # rak = ele["org"]["rank"]
                # else:
                    # rak = "NA"
                # if "parent_tax_id" in ele["org"]:
                    # pid = ele["org"]["parent_tax_id"] 
                # else:
                    # pid = "NA"
            # z = 1
            # [rak, sname, pid] = get_ptname(pid)
            # print (pid)
            # while rak != "PHYLUM" and z < 10:
                # if rak in ["FAMILY", "ORDER", "CLASS"]:
                    # attr.append(sname)
                # [rak, sname, pid] = get_ptname(pid)
                # z += 1
            # get assembly information
            for kew in ["assembly_accession", "org sci_name", "org common_name", "assembly_level", "submission_date", "contig_n50", "seq_length"]:
                kewlst = kew.split(" ")
                if len(kewlst) > 1:
                    attr.append(str(ele[kewlst[0]][kewlst[1]]) if kewlst[0] in ele and kewlst[1] in ele[kewlst[0]] else "NA") 
                else:
                    attr.append(str(ele[kewlst[0]]) if kewlst[0] in ele else "NA")
            attr.append("Y" if "annotation_metadata" in ele else "N")
            # print ("getSnt5")
            # print (get_sN50_others(attr[3])) 
            # print ("fins")
            attr.extend(get_sN50_others(attr[3]))
            # most the last element to the first
            # attr[0] = attr[-1]
            print ("\t".join(attr))
            counter += 1
            if counter % onep_step == 0:
                print ("Finished processing {0} records".format(counter, counter/n_records*100), file=sys.stderr)
        return 0
    else:
        print ("No resource found")
        return 1
def parse_asminfo(jfl, outf):

    jfl=sys.argv[1]
    outf=sys.argv[2]

    fout = open(outf, "w")

    with open(jfl) as f:
        asm_sumy = None
        if "datasets" in asm_sumy: 
            asm_sumy = json.load(f)['datasets']
        if asm_sumy is None:
            print ("Error: No assembly is available, please check species name or taxonomy ID")
            return 1

        # {'assembly_accession': 'GCA_000005465.1', 'display_name': 'BGIAF', 'org': {'tax_id': '9606', 'sci_name': 'Homo sapiens', 'common_name': 'human', 'sex': 'male', 'rank': 'SPECIES', 'parent_tax_id': '9605', 'assembly_counts': {'node': 128, 'subtree': 128}, 'key': '9606', 'title': 'human'}, 'chromosomes': ['Un'], 'assembly_level': 'Scaffold', 'submission_date': '2010-06-21', 'contig_n50': 887, 'estimated_size': '652893711', 'seq_length': '2676008911'}
        # print header
        hdr = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}".format("Accession ID", "Latin name", "Common name", "assembly level", "Submission date", "Contig N50", "Genome size", "Annotated")
        print (hdr, file=fout)
        for ele in asm_sumy:
            attr = []
            for kew in ["assembly_accession", "org sci_name", "org common_name", "assembly_level", "submission_date", "contig_n50", "seq_length"]:
                kewlst = kew.split(" ")
                if len(kewlst) > 1:
                    attr.append(str(ele[kewlst[0]][kewlst[1]]) if kewlst[0] in ele and kewlst[1] in ele[kewlst[0]] else "NA") 
                else:
                    attr.append(str(ele[kewlst[0]]) if kewlst[0] in ele else "NA")
            attr.append("Y" if "annotation_metadata" in ele else "N")
            print ("\t".join(attr), file=fout)
    fout.close()

def get_url(use_taxid, val):
    exturl = ""
    if use_taxid:
        exturl = "/assembly_descriptors/taxid/{}".format(val)
    else:
        exturl = "/assembly_descriptors/organism/{}".format(val)

    baseurl = "https://api.ncbi.nlm.nih.gov/datasets/v1alpha"
    sufurl = "?returned_content=COMPLETE"
    return "{0}{1}{2}".format(baseurl, exturl, sufurl)
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Parse assembly summary information')
    # g = parser.add_mutually_exclusive_group(required=True)
    parser.add_argument('-s',  action="store", dest = "species", help="species name")
    parser.add_argument('-t',  action="store", dest = "taxid", help="taxonomy ID")
    parser.add_argument('-p',  action="store", dest = "treepath", help="path of rankedlineage.dmp", default="rankedlineage.dmp")
    parser.add_argument('-n',  action="store", dest = "nrec", help="number of records, -1 for all", default="-1")
    args = parser.parse_args() 
    if args.species is None and args.taxid is None:
        print ("Require a species name or taxonomy ID")
        parser.print_help()
        exit(1)
    elif args.species is not None and args.taxid is not None:
        print ("Both species name and taxonomy ID are given, only use the taxonomy ID {}".format(args.taxid))
        use_taxid = 1
        url = get_url(use_taxid, args.taxid) 
    elif args.taxid is not None:
        use_taxid = 1
        url = get_url(use_taxid, args.taxid) 
    else:
        use_taxid = 0
        url = get_url(use_taxid, args.species) 
    tax_tree = load_tax_tree(args.treepath)
    nrec = args.nrec
    parse_asminfo2(url, tax_tree, nrec) 
